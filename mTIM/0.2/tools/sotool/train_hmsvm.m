function [information] = train_hmsvm(PAR)
% information = train_hmsvm(PAR)
%
% Wrapper for HM-SVMs and HMMs
%
% PAR -- a struct to configure the HM-SVM (for specification see
%   setup_hmsvm_training.m)
% returns a struct recording the training progress
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008
% modified by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2010

% add paths
this_name = mfilename();
this_path = mfilename('fullpath');
idx = find(this_path=='/', 1, 'last');
src_dir = this_path(1:idx);
fprintf('Hmsvm-source directory is %s.\n',src_dir);

% init sotool paths
init_paths();

% initialize and check PAR struct
PAR = init_par(PAR);

% prepare empty return value
progress = [];

% include user-specified include paths
if isfield(PAR, 'include_paths'),
  for i=1:length(PAR.include_paths),
    addpath(PAR.include_paths{i});
  end
end

% init state model
PAR.model_config = model_config();
disp(PAR);
disp(PAR.data_file);

% load data and select training examples
load(PAR.data_file, 'label', 'signal', 'exm_id');

if ~exist('signal', 'var'),
  fn_data = PAR.data_file;
  if isequal(fn_data(end-3:end), '.mat'),
    fn_data = fn_data(1:end-4);
  end
  fn_data = [fn_data '_signal'];
  signal = load_struct(fn_data, 'signal');
end

PAR.num_features = size(signal,1);

if ~exist('exm_id', 'var'),
 load(PAR.data_file, 'exm_id_intervals');
 assert(exist('exm_id_intervals', 'var') ~= 0);
else
  unq_exm_id = unique(exm_id);
  exm_id_intervals = zeros(length(unq_exm_id),3);
  for i=1:length(unq_exm_id),
    idx = find(exm_id==unq_exm_id(i));
    exm_id_intervals(i,:) = [unq_exm_id(i), idx(1), idx(end)];
  end
  clear exm_id
end
state_label = nan(size(label));

if isfield(PAR, 'train_exms'),
  % randomize order of potential training example before subselection
  train_exm_ids = PAR.train_exms;
  train_exm_ids = train_exm_ids(randperm(length(train_exm_ids)));
  train_exm_ids = train_exm_ids(1:PAR.num_train_exm);
  fprintf('\nusing %i sequences for training.\n', ...
          length(train_exm_ids));
  % for performance checks use sequences from validation set if given
  if isfield(PAR, 'vald_exms'),
    holdout_exm_ids = PAR.vald_exms;
  else
    holdout_exm_ids = [];
    fprintf('skipping performance estimation.\n\n');
    keyboard
  end
else
  % if training examples are not specified use all loaded sequences
  warning('No training set specified, treating whole data as training set!');
  assert(~isfield(PAR, 'vald_exms'));
  assert(~isfield(PAR, 'test_exms'));
  
  % randomize order of potential training example before subselection
  unq_exm_ids = unique(exm_id_intervals(:,1)');
  train_exm_ids = unq_exm_ids;
  train_exm_ids = train_exm_ids(randperm(length(train_exm_ids)));
  train_exm_ids = train_exm_ids(1:PAR.num_train_exm);
  fprintf('\nusing %i sequences for training.\n', ...
          length(train_exm_ids));
  % from the remainder take sequences for performance checks
  holdout_exm_ids = setdiff(unq_exm_ids, train_exm_ids);
  holdout_exm_ids = holdout_exm_ids(randperm(length(holdout_exm_ids)));
  assert(isempty(intersect(train_exm_ids, holdout_exm_ids)));
  fprintf('using %i sequences for performance estimation.\n\n', ...
          length(holdout_exm_ids));
end
% choose random subset for validation if there are too many
% validation examples
if length(holdout_exm_ids) > PAR.max_num_vald_exms,
  holdout_exm_ids = holdout_exm_ids(randperm(length(holdout_exm_ids)));
  holdout_exm_ids = holdout_exm_ids(1:PAR.max_num_vald_exms);
end
assert(isempty(intersect(train_exm_ids, holdout_exm_ids)));
fprintf('using %i sequences for performance estimation.\n\n', ...
        length(holdout_exm_ids));


%%%%% assemble model and score function structs,
LABELS = eval(sprintf('%s();', PAR.model_config.func_get_label_set));
state_model = eval(sprintf('%s(PAR);', PAR.model_config.func_make_model));
[score_plifs transition_scores] = eval(sprintf('%s(signal, label, state_model, PAR);', ...
                                               PAR.model_config.func_init_parameters));

assert(~any(isnan([score_plifs.limits])));
assert(~any(isnan([score_plifs.scores])));
assert(~any(isnan(transition_scores)));

%%%%% determine the true state sequence for each example from its label sequence
if isfield(PAR, 'label_noise_prop') && PAR.label_noise_prop > 0,
  noise_cnt = 0;
end
  warning('Let the state sequence decoding include problems.');
for i=1:size(exm_id_intervals,1),
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);

  if isfield(PAR, 'label_noise_prop') && PAR.label_noise_prop > 0,
    [label(idx) cnt] = add_label_noise(label(idx), PAR);
    noise_cnt = noise_cnt + cnt;
  end  
  true_label_seq = label(idx);
  obs_seq = signal(:,idx);
  [true_state_seq, problems] = eval(sprintf('%s(true_label_seq, state_model, obs_seq, PAR);', ...
                                            PAR.model_config.func_labels_to_states));
  % assert(isempty(problems));
  state_label(idx) = true_state_seq;
end

%{
% check whether there are roughly equally many positions per expression level
STATES = eval(sprintf('%s(PAR);', PAR.model_config.func_get_state_set));
num_p = zeros(2,PAR.num_levels);
for l=1:PAR.num_levels
  exw_state = getfield(STATES, sprintf('EXW_%02i', l));
  num_p(1,l) = sum(state_label == exw_state);

  exc_state = getfield(STATES, sprintf('EXC_%02i', l));
  num_p(2,l) = sum(state_label == exc_state);
end
%}

if isfield(PAR, 'label_noise_prop') && PAR.label_noise_prop > 0,
  fprintf('  converted %i segments (label noise level: %2.1f%%)\n', noise_cnt, ...
          100*PAR.label_noise_prop);
end


% PRE-PROCESS training data
% if switch features is on
if (isfield(PAR,'remove_low_cover_blocks') && PAR.remove_low_cover_blocks),
    % filter out MIN_COVER_BLOCKS>CFG.gene_states_low_cover_cutoff
    config = model_config();
    % check if feature description is available otherwise
    % suppression would not work
    if (~isfield(config,'func_get_feature_set')),
        error('No feature description found in model.\n');
    end
    % convert to function handle and get the feature description
    get_feature_desc = str2func(config.func_get_feature_set);
    [FEATS, FEAT_MAP] = get_feature_desc();
    LOW_COVER_FEAT = strmatch('low_cover_blocks', FEATS, 'exact');

    % retain these positions
    fprintf('Remove low-cover-blocks that have a value below %2.1f.\n',PAR.gene_states_low_cover_cutoff);
    retain_inds = find(signal(LOW_COVER_FEAT,:)<=PAR.gene_states_low_cover_cutoff);
    rm_inds = setdiff(1:size(signal,2),retain_inds);

    fprintf('%i positions marked for deletion.\n',length(rm_inds));
    fprintf('%i positions remain.\n',length(retain_inds));

    % new observation vector
    fprintf('Re-arranging observation matrix...');
    signal = signal(:,retain_inds);
    label = label(retain_inds);
    state_label = state_label(retain_inds);
    fprintf('Done!\n');

    fprintf('Build new example intervals vector...');
    new_exm_id_intervals = [];
    start = 1;
    stop = -1;
    for i=1:size(exm_id_intervals,1);
        id    = exm_id_intervals(i,1);
        o_start = exm_id_intervals(i,2);
        o_stop  = exm_id_intervals(i,3);

        inds = retain_inds(find(retain_inds>=o_start & retain_inds<=o_stop));
        stop = start+length(inds)-1;
        % new intervals
        new_exm_id_intervals = [new_exm_id_intervals; [id,start,stop]];
        if (length(inds)>0), start = stop+1; end
    end
    exm_id_intervals = new_exm_id_intervals;
    fprintf('Done!\n');

    % remove all zero length intervals
    lens = abs(exm_id_intervals(:,3)-exm_id_intervals(:,2)+1);
    fprintf('All chunks have a length>0 = %i\n',all(lens>0));
    if (any(lens==0)),
        ids = exm_id_intervals(find(lens>0),1);
        fprintf('Retain %i of %i chunks.\n',length(ids),length(lens));
        train_exm_ids = intersect(train_exm_ids, ids);
        holdout_exm_ids = intersect(holdout_exm_ids, ids);
    end

    % remove all remaining sequences that cannot be decoded 
    % by the current state model

    % remove index and ids of the chunks that cannot be decoded
    fprintf('Check for non-decodable state transitions in the training set.\n');
    rm_inds = [];
    rm_ids = [];
    % change in the label is marked by 1
    for i=1:size(exm_id_intervals,1),
        id = exm_id_intervals(i,1);
        o_start = exm_id_intervals(i,2);
        o_stop  = exm_id_intervals(i,3);
      
        s_l = state_label(o_start:o_stop);
        assert(~any(isnan(s_l))); % sanity check
        grad = s_l(2:end)~=s_l(1:end-1);
        inds = find(grad==1);

        % ind and the corresponding successor state differ
        % check all pairs 
        for j=1:length(inds),
            si = s_l(inds(j));
            sj = s_l(inds(j)+1);
            ind = find([state_model.id]==si);
            succs = state_model(ind).successors;
            if (~any(succs==sj)),
                rm_inds = [rm_inds, i];
                rm_ids = [rm_ids, id];
                continue;
            end
        end
    end
    fprintf('There are %i remaining examples that have non-decodable state transitions.\n', ...
        length(rm_inds));
    train_exm_ids = setdiff(train_exm_ids, rm_ids);
    holdout_exm_ids = setdiff(holdout_exm_ids, rm_ids);
    fprintf('%i training examples and %i validation examples remain.\n', ...
        length(train_exm_ids),length(holdout_exm_ids));
end


% filter dimensions (if suppress features method is present in model)
if (exist('suppress_features','file')),
    PAR.dbg_signal_before_suppress = check_signal(PAR, signal);
    [PAR, signal, score_plifs, state_model] = suppress_features(PAR, signal, score_plifs, state_model);
    PAR.dbg_signal_after_suppress = check_signal(PAR, signal);
    fprintf('Histograms of each signal before and after suppressing features.\n');
    PAR.dbg_signal_before_suppress
    PAR.dbg_signal_after_suppress
else
    warning('Function `suppress_features` is not defined in model.');
end
 


% display training sequence size and further information
idx = find(ismember(exm_id_intervals(:,1), train_exm_ids));
L = sum(exm_id_intervals(idx,3) - exm_id_intervals(idx,2) + length(idx));
lens = abs(exm_id_intervals(idx,3)-exm_id_intervals(idx,2)+1);
fprintf('\nTraining on sequences with a total length of %i nt.\n', L);
fprintf('Shortest example has a size of %i nt.\n', min(lens));
fprintf('Longest example has a size of %i nt.\n', max(lens));
fprintf('Median example has a size of %i nt.\n', median(lens));
fprintf('Mean example size is %i nt.\n', mean(lens));

if (isfield(LABELS,'ambiguous')),
    ambs = sum(label==LABELS.ambiguous);
    fprintf('%i(=%2.2f%%) ambigiuous labels present.\n',ambs,100.0*ambs/L);
else
    fprintf('No ambiguous labels present.\n');
end



% choose correct classifier according
% to PAR.classifier 
classifier = @train_primal_sosvm;
if (isfield(PAR,'classifier'))
  classifier = PAR.classifier;
else
  warning('PAR.classifier is not defined. Chose default instead.');
end
fprintf('Chosen classifier: %s \n',func2str(classifier));

% model is part of the argument list
model.state_model = state_model;
model.transition_scores = transition_scores;
model.score_plifs = score_plifs;
model.LABELS = LABELS;

% ..as well as examples, labels and holdout
trainset.num_examples = length(train_exm_ids);
trainset.train_exm_ids = train_exm_ids
trainset.holdout_exm_ids = holdout_exm_ids;
trainset.exm_id_intervals = exm_id_intervals;
trainset.signal = signal;
trainset.label = label;
trainset.state_label = state_label;

% train classifier
[information] = classifier(PAR, model, trainset);

% add information (cumulants) about the signal sequences
% which is neccessary for 'zero-shot transfer learning'
[sig_means,sig_stds,sig_valid] = get_signal_cumulants(signal);
PAR.transfer.sig_means = sig_means;
PAR.transfer.sig_stds = sig_stds;
PAR.transfer.sig_valid = sig_valid;

% for debugging save model in extra text file
save_model_dbg(PAR, information, PAR.out_dir);




function bins = check_signal(PAR, signal)
dims = size(signal,1);
levels = 5;
bins = zeros(dims,levels);
for i=1:dims,
    smin = min(signal(i,signal(i,:)>-inf));
    smax = max(signal(i,signal(i,:)<+inf));

    lvl = linspace(smin,smax,levels);

    foo = 0;
    for j=1:levels,
        bins(i,j) = length(find(signal(i,:)<=lvl(j)));
        foo = foo+bins(i,j);
    end
    bins(i,:) = bins(i,:)./foo;
end
