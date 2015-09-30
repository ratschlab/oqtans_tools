function save_model_dbg(PAR, info, fout)
% Saves the trained model as a human-readable file
% for further inspection.
%
% info :
% fout : (optional) output directory
%
% written by Nico Goernitz, 2012

% find directory
foutput = '';
if (exist('fout','var')), foutput = fout; end

% write to file
fh = fopen([foutput 'dbg_model.txt'],'w');

model = info{end}.model;
config = model_config();

% check if feature description is available otherwise
% suppression would not work
if (~isfield(config,'func_get_feature_set')),
    fprintf('No feature description found in model.\n');
    return;
end

% convert to function handle and get the feature description
get_feature_desc = str2func(config.func_get_feature_set);
[FEATS, FEAT_MAP] = get_feature_desc();
% if features have been suppressed 'PAR.suppressed_feats_map' contains
% the real index in the FEATS array
if isfield(PAR,'suppressed_feats_map'),
    fprintf(fh,'Suppressed_feats_map found in PAR.\n');
    map = PAR.suppressed_feats_map;
    FEATS = FEATS(map(:,2)>0);
else
    fprintf(fh,'Warning! No suppressed_feats_map found in PAR. The model is an older version.\n');
    fprintf(fh,'If features have been suppressed during training the following feature\n');
    fprintf(fh,'list per state might be not correct.\n\n');
end


% single entries
state_model = model.state_model;
plifs = model.score_plifs;

% check for non-decodable state transitions in the training set
if (isfield(info{end},'dbg')),
    fprintf('Check for non-decodable state transitions in the training set.\n');
    dbg = info{end}.dbg;
    % the current state
    true_seq = dbg.train_state_seq;

    % remove index and ids of the chunks that cannot be decoded
    pos = [];
    % change in the label is marked by 1
    grad = true_seq(2:end)~=true_seq(1:end-1);
    inds = find(grad==1);

    % ind and the corresponding successor state differ
    % check all pairs 
    for j=1:length(inds),
        si = true_seq(inds(j));
        sj = true_seq(inds(j)+1);
        ind = find([state_model.id]==si);
        succs = state_model(ind).successors;
        if (~any(succs==sj)), pos=[pos, inds(j)]; end
    end
    fprintf('There are %i transitions in the training sequence of which %i are non-decodable.\n', ...
        length(inds), length(pos));

end

em_dim_count = 0;
for s=1:length(state_model),

    % current state
    state = state_model(s);
    
    % check for existing dbg field which enables
    % us to check the accuracy for the specific state 
    acc_txt = '';
    if (isfield(info{end},'dbg')),
        dbg = info{end}.dbg;
        % the current state
        true_seq = dbg.train_state_seq==state.id;
        pred_seq = dbg.train_pred_state_seq==state.id;
        
        acc = sum(true_seq & pred_seq) / sum(true_seq);
        acc_abs = sum(true_seq & pred_seq);
        acc_txt = sprintf('true:%int pred:%int acc:%i(=%2.2f)',...
            sum(true_seq),sum(pred_seq),acc_abs,acc);
    end


    % first line: state_name: 
    fprintf(fh,'%s(%s):\n',state.name,acc_txt);

    % plot feature information
    learn_inds = find(state.learn_scores~=0);
    learn_inds = sort(learn_inds,'ascend'); % for monot_scores the sorting is important
    learn_names = FEATS(learn_inds);
    for i=1:length(learn_names),
        monot_sign = 'x'; % no monotonicity defined
        if (state.monot_scores(i)==+1), monot_sign='/'; end;
        if (state.monot_scores(i)==-1), monot_sign='\'; end;

        fprintf(fh,'  -%s(%s): (',learn_names{i},monot_sign);
        scores = plifs(learn_inds(i),s).scores;
        em_dim_count = em_dim_count + length(scores);
        for j=1:length(scores)-1;
            fprintf(fh,'%2.2f,',scores(j));
        end
        if (~isempty(scores)),
            fprintf(fh,'%2.2f',scores(end));
        end
        fprintf(fh,')\n');
    end


    fprintf(fh,'\n');
end

fprintf('Number of emission dimensions: %i\n',em_dim_count);
fprintf('Number of states: %i\n',length(state_model));
fprintf('Problem dimensionality: %i\n',length(info{end}.w));
