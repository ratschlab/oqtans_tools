function train_orgs_parallel(base_dir, in_dirs, out_dir, num_train, num_vald, params_c, ...
  params_b, mtl_filename)
% Load all the organisms located in 'base_dir/in_dirs'. Choose
% for each organism 'num_train' training data and 'num_val' validation data.
% Last half of each organism's data is used for testing purposes.
% Use params_c and params_b for model selection.
%
% Attention: the first half of the data is for training/validation
% purposen while the last half is used for testing.
% Therefore, never ever train/validate on the last half of the data.
%
% IN
%   params_c     : [1 x n_c] parameters for C_small 
%   params_b     : [1 x n_b] parameters for MTL b
%   mtl_filename : (optional) if file exist use MTL learning
%
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

addpath('model');
addpath('rproc');
addpath('utils2');

addpath('../hmsvm');
addpath('../hmsvm/native');

addpath('../../src');
addpath('../../src/linesearch');
addpath('../../src/losses');
addpath('../../src/solver/prloqo');

model_name = 'bacteria';

out_dir = [base_dir '/' out_dir];
if (~exist(out_dir,'dir')),
  fprintf('Create new directory %s.\n',out_dir);
  mkdir(out_dir);
end

% parameter combinations to be used for independent training
parameters = build_params_vector({params_c, params_b});

mtl_PAR = [];

% if multi-task is enabled
if (exist('mtl_filename','var') && exist(mtl_filename,'file')),
  % load the final wstar 
  fprintf('Load MTL file %s ...',mtl_filename);
  load(mtl_filename);
  fprintf('finished.\n');
  mtl_PAR = PAR;
  PAR = [];
  % a w-vector should exists -> this becomes our
  % mtl_wstar
  assert(exist('w','var')>=1);
  PAR.mtl.mtl_enable = 1; 
  PAR.mtl.mtl_wstar = w;
  PAR.mtl.filename = mtl_filename;
else
  PAR.mtl.mtl_enable = 0; 
  PAR.mtl.mtl_b = 0;
  PAR.mtl.mtl_wstar = [];
  parameters = build_params_vector({params_c, []});
end

PAR.parameters = parameters;
PAR.num_repetitions = 5;
PAR.num_features = 1;
PAR.model_config = model_config;

% final result vector
result = cell(PAR.num_repetitions,size(parameters,1));

% number of training sessions
cnt = 1;

% maximum of parallel jobs allowed
MAX_JOBS = 30;

% repeat sampling
for n = 1:PAR.num_repetitions,

  [PAR, label, signal, state_label, exm_id_intervals] = ...
      load_multiple_organisms(PAR, base_dir, in_dirs, num_train, num_vald);
  data_filename = [out_dir sprintf('/data_%i.mat',n)];
  fprintf('Saving data set as %s.\n',data_filename);
  save(data_filename,'label','signal','state_label','exm_id_intervals','PAR');    
    
  for i=1:size(parameters,1),
    % store current repetition/parameter index
    PAR.n = n;
    PAR.i = i;
    PAR.check_acc = 0;
    
    % constant parameters
    PAR.out_dir = [out_dir '/rep' num2str(n) '_param' num2str(i) '/']; % output directory
    PAR.num_plif_nodes = 64;                         % number of supporting points

    % and current parameter
    PAR.C_small = parameters(i,1);
    % set coupling and smoothing to zero:
    % while coupling can be (and is) disabled by the model
    % smoothing can't be disabled therefore it is important
    % to set it to 0.0
    PAR.C_coupling = 0.0;
    PAR.C_smooth = 0.0;
    % parameter for convex combination for the additional
    % 2-norm distance regularizer:
    % mtl_b * w'Qw + (1-mtl_b) * ||w-mtl_wstar||^2
    if (PAR.mtl.mtl_enable),
      PAR.mtl.mtl_b = parameters(i,2);
    end
    
    % set the permitted feature ranges such that state transitions
    % are enforced
    PAR.perm_feature_ranges = set_permitted_feature_ranges();

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
    trainset.exm_id_intervals = exm_id_intervals;
    trainset.signal = signal;
    trainset.label = label;
    trainset.state_label = state_label;

    % start training
    disp(PAR)
    %time = clock();
    %[information] = train_primal_sosvm(PAR, model, trainset);
    %fprintf('Complete training took %3.2f sec\n',etime(clock(), time));
    
    PAR.model=model;
    PAR.trainset=trainset;
    opts.force_matlab = 1;
    jobinfo(cnt) = rproc('train_orgs_parallel_atom', PAR, 6000, opts, 100000);
    cnt = cnt+1;

		% if job counter exceeds MAX_JOBS then wait until they are finished
		if (cnt>MAX_JOBS),
			fprintf('Waiting...');
			jobinfo = rproc_wait(jobinfo, 20, 1, -1);
			cnt = 1;
			fprintf('Done!\n');
		end

    %store the final result
    %result{n,i} = information{end};
  end
end

% wait for job finished
fprintf('Waiting...');
jobinfo = rproc_wait(jobinfo, 20, 1, -1);
fprintf('Done!\n');


% collect the results
fprintf('Collecting results...');
for n = 1:PAR.num_repetitions,
  for i=1:size(parameters,1),
    curr_out_dir = [out_dir '/rep' num2str(n) '_param' num2str(i) '/']; % output directory
    % information about the current directory
    cdir = dir(curr_out_dir);
    
    info = {};
    for j = 1:length(cdir),
      if (~cdir(j).isdir && ~isempty(strfind(cdir(j).name,'final'))),
        fprintf('Load %s.\n',[curr_out_dir cdir(j).name]);
        load([curr_out_dir cdir(j).name]);
      end
    end

    information{end}.PAR = PAR;
    information{end}.isConverged = isConverged;
    
    result{n,i} = information{end};
  end
end
fprintf('Done!\n');

% store the result file
fprintf('Saving final results to %s/final.mat.\n\n',out_dir);
save('-v7.3', [out_dir '/final.mat'],'result','PAR','mtl_PAR');
% eof
