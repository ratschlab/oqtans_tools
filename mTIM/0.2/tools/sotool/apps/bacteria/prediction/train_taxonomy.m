function train_taxonomy(base_dir, fname, params_c, params_b, max_runs)
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

addpath('..');
addpath('../model');
addpath('rproc');
addpath('utils2');

addpath('../../hmsvm');
addpath('../../hmsvm/native');

addpath('../../../src');
addpath('../../../src/linesearch');
addpath('../../../src/losses');
addpath('../../../src/solver/prloqo');

model_name = 'bacteria';

load([base_dir fname]);

num_runs = min(size(runs,1),max_runs); 

whos

% PARALLEL TRAINING
PARALLEL = 0;
% number of training sessions
cnt = 1;
% maximum of parallel jobs allowed
MAX_JOBS = 30;

if (PARALLEL),
  fprintf('Parallel training is enabled.\n');
end

% parameter combinations to be used for independent training
% check if mtl_b contains 1.0 as parameter (=independet training) 
if (~any(params_b==1.0)),
    fprintf('Adding mtl_b=1.0 as a possible parameter (=independent training).\n');
    params_b = [params_b, 1.0];
end
params_mtl = build_params_vector({params_c, params_b});
params_root = build_params_vector({params_c, 0});

% convert the uint8 arrays into double again
signal = double(signal);
label = double(label);

% setup new exm id intervals 
new_exm_id_intervals = [example_id', exm_id_intervals'];

% store final fscores test values for each taxonomy branch
% 1 :union 
% 2 :independent
% 3-:single organism
num_single_orgs = size(org_idx_names,2); 
fscores = zeros(size(taxonomy,1),2+num_single_orgs);
fscores_std = zeros(size(taxonomy,1),2+num_single_orgs); 

% mark finished nodes (for parallel training)
% -1:not started yet
%  0:already started
% +1:finished
is_node_finished = -ones(1,size(taxonomy,1));

% for all taxonomy nodes
for t=1:size(taxonomy,1),
    fprintf('Train node %s (depends on node %i).\n',tax_idx_names{t}, taxonomy(t,end));

    % create a directory with the node name
    out_dir = ([base_dir tax_idx_names{t}]);
    mkdir(out_dir);
   
    % load organisms binary vector
    org_mask = taxonomy(t,1:end-1);
    org_inds = find(org_mask==1);
    num_orgs = sum(org_mask==1);

    % if there is a parent node load 
    % the mtl parameter set
    params = params_root;
    % convert from python 0..T-1 to matlab 1..T
    root_idx = taxonomy(t,end) + 1;
    if (root_idx>0),
       % load mtl parameter set
        params = params_mtl;
    end

    % going to solve node t
    is_node_finished(t) = 0;

    % store the results for the current node
    results = {};
    evals = {};
    fscore_vald = zeros(num_runs,size(params,1));
    fscore_test = zeros(num_runs,num_single_orgs,size(params,1));

    % if the next job depends on an unfinished parent node
    % wait until parent node is finished
    if (PARALLEL & root_idx>0 & is_node_finished(root_idx)<1),
      fprintf('Waiting...');
      jobinfo = rproc_wait(jobinfo, 20, 1, -1);
      cnt = 1;
      fprintf('Done!\n');
    end

    % for all repetitions
    for r=1:num_runs,

        % load the best parent classifier for the current run
        % if there is a root node, load solution w_star 
        PAR.mtl.mtl_enable = 0; 
        PAR.mtl.mtl_b = 1.0;
        PAR.mtl.mtl_wstar = [];
        if (root_idx>0),
            % assert that the training of the parent node is already finished
            assert(is_node_finished(root_idx)==1);

            % get the parent directory
            parent_dir = [base_dir tax_idx_names{root_idx}];
            % load the results to get best parameter
            foo = load([parent_dir '/final.mat']);
            opt_param_ind = foo.opt_param_ind;
            % load the parent result w_star of optimal parameter and current run r 
            %foo = load([parent_dir '/rep' num2str(r) '_param' num2str(opt_param_ind) '/sosvm_final.mat']); % output directory
 
            % load the solution vector
            PAR.mtl.mtl_enable = 1; 
            PAR.mtl.mtl_wstar = foo.result{r,opt_param_ind}.w;
        end

        % setup training/validation and test data
        % load slice
        slice = reshape(runs(r,org_mask==1,:), num_orgs, num_train+num_vald+num_test);
        % extract ids 
        PAR.num_train_exm = num_orgs * num_train;
        PAR.train_exms = reshape(slice(:,1:num_train),1,num_orgs*num_train); 
        PAR.vald_exms = reshape(slice(:,num_train+1:num_train+num_vald),1,num_orgs*num_vald)  
        PAR.test_exms = reshape(slice(:,num_train+num_vald+1:end),1,num_orgs*num_test) 
        
        % for all parameters
        for p=1:size(params,1),
            fprintf('Run(%i/%i) params(%i/%i) #Orgs(%i).\n',r,size(runs,1),p,size(params,1),num_orgs);

            % train so-svm
    
            % store current repetition/parameter index
            PAR.n = r;
            PAR.i = p;
            PAR.check_acc = 0;
            
            % constant parameters
            PAR.out_dir = [out_dir '/rep' num2str(r) '_param' num2str(p) '/']; % output directory
            PAR.num_plif_nodes = 64;                         % number of supporting points
            PAR.num_features = 1;
            PAR.model_config = model_config;

            % and current parameter
            PAR.C_small = params(p,1);
            % set coupling and smoothing to zero:
            % while coupling can be (and is) disabled by the model
            % smoothing can't be disabled therefore it is important
            % to set it to 0.0
            PAR.C_coupling = 0.0;
            PAR.C_smooth = 0.0;
            % parameter for convex combination for the additional
            % 2-norm distance regularizer:
            % mtl_b * w'Qw + (1-mtl_b) * ||w-mtl_wstar||^2
            PAR.mtl.mtl_b = params(p,2);
    

            % set the permitted feature ranges such that state transitions
            % are enforced
            PAR.perm_feature_ranges = set_permitted_feature_ranges();

            % init PAR vector with default values
            PAR = init_par(PAR);

            % model is part of the argument list
            state_model = make_model(PAR);
            [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR);

            % set model related data
            model.state_model = state_model;
            model.transition_scores = transition_scores;
            model.score_plifs = score_plifs;
            model.LABELS = get_label_set();

            % ..as well as examples, labels and holdout
            trainset.num_examples = PAR.num_train_exm;
            trainset.train_exm_ids = PAR.train_exms;
            trainset.holdout_exm_ids = PAR.vald_exms;
           
            trainset.exm_id_intervals = new_exm_id_intervals;
            trainset.signal = signal;
            trainset.label = label;
            trainset.state_label = state_label;

            % start training
            if (~PARALLEL),
              % sequential training
              time = clock();
              [information] = train_primal_sosvm(PAR, model, trainset);
              fprintf('Complete training took %3.2f sec\n',etime(clock(), time));
              %store the final result
              %result{r,p} = information{end};
            else
              % train parallel and collect the results
              % later on
              PAR.model=model;
              PAR.trainset=trainset;
              opts.force_matlab = 1;
              jobinfo(cnt) = rproc('train_taxonomy_atom', PAR, 6000, opts, 100000);
              cnt = cnt+1;

              % if job counter exceeds MAX_JOBS then wait until they are finished
              if (cnt>MAX_JOBS),
                fprintf('Waiting...');
                jobinfo = rproc_wait(jobinfo, 20, 1, -1);
                cnt = 1;
                fprintf('Done!\n');
              end
            end
            
       end
    end


    if (PARALLEL),
        % wait for job finished
        fprintf('Waiting for jobs to finish...');
        jobinfo = rproc_wait(jobinfo, 20, 1, -1);
        fprintf('Done!\n');
    end

    % collect the results
    fprintf('Collecting results...');
    for n = 1:num_runs,
        for i=1:size(params,1),
          curr_out_dir = [out_dir '/rep' num2str(n) '_param' num2str(i) '/']; % output directory
          % information about the current directory
          cdir = dir(curr_out_dir);

          info = {};
          for j = 1:length(cdir),
            if (~cdir(j).isdir && ~isempty(strfind(cdir(j).name,'final'))),
              fprintf('Load %s.\n',[curr_out_dir cdir(j).name]);
              foo = load([curr_out_dir cdir(j).name]);
              PAR = foo.PAR;
              information = foo.information;
            end
          end
          result{n,i} = information{end};
          
          % validate 
          res_vald = evaluate(PAR, PAR.vald_exms, information{end}.model, new_exm_id_intervals, signal, label);
          evals{n,i}.vald = res_vald;
          fscore_vald(n,i) = res_vald.f_score;

          % test with respect to the single organism performance
          for o=org_inds,
              inds = find(org_id(PAR.test_exms)==(o-1));
              res_test = evaluate(PAR, PAR.test_exms(inds), information{end}.model, new_exm_id_intervals, signal, label);
              fscore_test(n,o,i) = res_test.f_score;
          end
          evals{n,i}.test = reshape(fscore_test(n,:,i),1,num_single_orgs);
        end
    end
    fprintf('Done!\n');


    % select the best classifier based on the validation data 
    mean_vald = mean(fscore_vald,1);
    mean_test = reshape(mean(fscore_test,1),num_single_orgs,size(params,1));
    [foo ind] = max(mean_vald);    

    % if fscore_* has only one row, than 
    % std should be zero
    std_vald = zeros(1,length(mean_vald));
    std_test = zeros(num_single_orgs,size(params,1));
    if (size(mean_vald,1)>1),
        std_vald  = std(fscore_vald);

        for o=1:num_single_orgs,
            std_test(o,:) = std(reshape(fscore_test(:,o,:),num_runs,size(params,1)));
        end
    end

    % union performance
    fscores(t,1) = mean(mean_test(org_inds,ind));
    fscores_std(t,1) = mean(std_test(org_inds,ind));

    % single performance
    fscores(t,3:end) = mean_test(:,ind);
    fscores_std(t,3:end) = std_test(:,ind);

    % independent
    if (root_idx>0),
        % also get independent-performance (where mtl_b = 1.0) 
        % find all params where params(:,end) = 1.0
        ninds = find(params(:,end)==1.0);
        [foo nind] = max(mean_vald(ninds));
        nind = ninds(nind);
        assert(~isempty(ninds));

        fscores(t,2) = mean(mean_test(org_inds,nind));
        fscores_std(t,2) = mean(std_test(org_inds,nind));
    end

    % store results
    opt_param_ind = ind;
    save([out_dir '/final.mat'], 'result','fscore_vald','fscore_test','evals','opt_param_ind','fscores','fscores_std');
    
    % going to solve node t
    is_node_finished(t) = +1;
    
    % mark the current node as finished
    is_node_finished(t) = 1;
end


% print results
fprintf('\n\n==================F-SCORE RESULTS================\n');
fprintf('Performance: Union, Independent, Org1, Org2, ...\n');
for t=1:size(taxonomy,1),
    name = tax_idx_names{t};

    fprintf('f-score=');
    for i=1:size(fscores,2),
        score = '. . . ';
        score_std = ' . .';
        if (fscores(t,i)>0.0), 
            score = sprintf('%1.4f',fscores(t,i)); 
            score_std = sprintf('%1.2f',fscores_std(t,i));
        end

        delim = ', ';
        if (i<=2), delim='| '; end

        fprintf('%s+%s%s',score,score_std,delim);
    end
    fprintf(' for %s\n',name);
end








function res = evaluate(PAR, ids, model, exm_id_intervals, signal, label)
% Evaluate the example with id 'ids'.

% get valid labels
LABELS = get_label_set();

tp=0; tn=0; fp=0; fn=0;
for j=1:length(ids),
  % get the index of the current example id
  idx = find(exm_id_intervals(:,1) == ids(j));
  % start and stop position in the sequence
  % and observation sequence and corresponding label sequence
  start=exm_id_intervals(idx,2);
  stop=exm_id_intervals(idx,3);
  obs_seq = signal(start:stop);
  label_seq = label(start:stop);    

  % prediction
  [pred_path] = decode_viterbi(obs_seq, model.transition_scores, model.score_plifs, PAR);
  
  % evaluate the nukleotide exonic performance
  [f1,f2,f3,f4] = eval_performance(label_seq, pred_path.label_seq, LABELS.exonic);
  tp=tp+f1; tn=tn+f2; fp=fp+f3; fn=fn+f4;
end
% get sensitivity, specificity, f_score, precision
res = eval_confusion_matrix(tp,tn,fp,fn);

% eof

