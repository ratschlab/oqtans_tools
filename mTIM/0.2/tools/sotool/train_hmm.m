function [information] = train_hmm(PAR, model, trainset)
% Trains a baseline HMM.
%
% PAR -- a struct to configure the HMM (for specification see
%   setup_hmsvm_training.m)
% returns a struct summarizing training accuracy
%
% see train_hmsvm.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008
% adapted by Nico Goernitz, Berlin Institute of Technology, 2012



% model is part of the argument list
state_model = model.state_model;
transition_scores = model.transition_scores;
score_plifs = model.score_plifs;
LABELS = model.LABELS;

% ..as well as examples, labels and holdout
train_exm_ids = trainset.train_exm_ids;
holdout_exm_ids = trainset.holdout_exm_ids;
exm_id_intervals = trainset.exm_id_intervals;
signal = trainset.signal;
label = trainset.label;
state_label = trainset.state_label;



%%%%% start training
% a struct keeping track of training progress
progress = [];
% accuracy on training examples
trn_acc = zeros(1,length(train_exm_ids));
% accuracy on holdout validation examples
val_acc = zeros(1,length(holdout_exm_ids));
% record elapsed time
t_start = clock();

tic
% compute mapping (features, states) -> position in one-dimensional vector
PAR.num_trans_score = length(transition_scores);
next_score_start = length(transition_scores)+1;
res_map = zeros(PAR.num_features, length(state_model));
cnt = 1;
for i=1:length(state_model), % for all states
  idx = find(state_model(i).learn_scores);
  for j=1:length(idx),
    row_idx = state_model(i).feature_scores(j,1);
    col_idx = state_model(i).feature_scores(j,2);
    if res_map(row_idx, col_idx) == 0,
      res_map(row_idx, col_idx) = next_score_start;
      score_starts(cnt) = next_score_start;
      next_score_start = next_score_start + PAR.num_plif_nodes;
      cnt = cnt + 1;
    end
  end
end
PAR.num_param = next_score_start - 1;
assert(PAR.num_features == size(score_plifs,1));
assert(length(state_model) == size(score_plifs,2));
assert(PAR.num_features == size(res_map,1));
assert(length(state_model) == size(res_map,2));


%%% do one round of Viterbi decoding for all training examples to count
%%% transitions and feature values
W = zeros(1, PAR.num_param);
for i=1:length(train_exm_ids),
  idx = find(exm_id_intervals(:,1)==train_exm_ids(i));
  idx = exm_id_intervals(idx,2):exm_id_intervals(idx,3);
  obs_seq = signal(:,idx);
  true_label_seq = label(idx);
  true_state_seq = state_label(idx);
    
  [tmp true_path] = decode_viterbi(obs_seq, transition_scores, score_plifs, ...
                                   PAR, true_label_seq, true_state_seq);
  W = W + weights_to_vector(true_path.transition_weights, ...
                            true_path.plif_weights, state_model, ...
                            res_map, PAR);
end
fprintf('Decoded %i training sequences\n\n', length(train_exm_ids));
fprintf('Decoding took %3.2f sec\n\n', toc);


%%% add pseudo counts
if ~isfield(PAR, 'hmm_pseudo_cnt');
  pseudo_cnt = 1;
else
  pseudo_cnt = PAR.hmm_pseudo_cnt;  
end
W(W<pseudo_cnt) = pseudo_cnt;
% this does not work as well
%W = W + pseudo_cnt;

%%% compute transition scores as transition frequencies
for i=1:length(state_model),
  idx = state_model(i).trans_scores;
  % LOG-TRANSFORM scores as Viterbi assumes ADDITIVE scores
  if (idx>0)
    transition_scores(idx) = log(W(idx) ./ sum(W(idx)));
  end
end

%%% compute feature scoring function values as frequencies of occurrence in
%%% true paths (accumulated in W)
for i=1:size(res_map,1), % for all features
  for j=1:size(res_map,2), % for all states
    if res_map(i,j) ~= 0,
      idx = res_map(i,j):res_map(i,j)+PAR.num_plif_nodes-1; 
      % LOG-TRANSFORM scores as Viterbi assumes ADDITIVE scores
      score_plifs(i,j).scores = log(W(idx) ./ sum(W(idx)));
    end
  end
end

progress.el_time = etime(clock(), t_start);


%%% check prediction accuracy on training examples
for j=1:length(train_exm_ids),
  trn_idx = find(exm_id_intervals(:,1)==train_exm_ids(j));
  trn_idx = exm_id_intervals(trn_idx,2):exm_id_intervals(trn_idx,3);
  trn_obs_seq = signal(:,trn_idx);
  trn_true_label_seq = label(trn_idx);
  trn_true_state_seq = state_label(trn_idx);
  trn_pred_path = decode_viterbi(trn_obs_seq, transition_scores, score_plifs, PAR);
  trn_pred_label_seq = trn_pred_path.label_seq;
  trn_acc(j) = mean(trn_true_label_seq(1,:)==trn_pred_label_seq(1,:));
  
  if PAR.verbose>=3 && j<=25,
    view_label_seqs(gcf, trn_obs_seq, trn_true_label_seq, trn_pred_label_seq);
    title(gca, ['Training example ' num2str(train_exm_ids(j))]);
    fprintf('Training example %i\n', train_exm_ids(j));
    fprintf('  Example accuracy: %3.2f%%\n', 100*trn_acc(j));
    pause
  end
end
fprintf('  HMM training accuracy:              %2.2f%%\n', ...
        100*mean(trn_acc));
progress.trn_acc = trn_acc';


%%% check prediction accuracy on holdout examples
if ~isempty(holdout_exm_ids),
  for j=1:length(holdout_exm_ids),
    val_idx = find(exm_id_intervals(:,1)==holdout_exm_ids(j));
    val_idx = exm_id_intervals(val_idx,2):exm_id_intervals(val_idx,3);
    val_obs_seq = signal(:,val_idx);
    val_pred_path = decode_viterbi(val_obs_seq, transition_scores, score_plifs, PAR);
    val_true_label_seq = label(val_idx);
    val_true_state_seq = state_label(val_idx);
    val_pred_label_seq = val_pred_path.label_seq;
    val_acc(j) = mean(val_true_label_seq(1,:)==val_pred_label_seq(1,:));
    
    if PAR.verbose>=3 && j<=25,
      view_label_seqs(gcf, val_obs_seq, val_true_label_seq, val_pred_label_seq);
      title(gca, ['Hold-out example ' num2str(holdout_exm_ids(j))]);
      fprintf('Hold-out example %i\n', holdout_exm_ids(j));
      fprintf('  Example accuracy: %3.2f%%\n', 100*val_acc(j));
      pause
    end
  end
  fprintf('  HMM validation accuracy:            %2.2f%%\n\n', ...
          100*mean(val_acc));
  progress.val_acc = val_acc';
end

if PAR.verbose>=3,
  plot_progress(progress, fh1);
  pause(1);
end    

if PAR.verbose>=3,
  eval(sprintf('%s(state_model, score_plifs, PAR, transition_scores);', ...
               PAR.model_config.func_view_model));
end  

% save all into model 
model.state_model = state_model;
model.transition_scores = transition_scores;
model.score_plifs = score_plifs;

information{1}.w = W'; 
information{1}.model = model;
% plot(W)

% save results
fprintf('Saving result...\n\n\n');
% fname = sprintf('hmm_training_minl%i', PAR.hmm_min_level);
fname = sprintf('sosvm_final.mat');
save([PAR.out_dir fname], 'information','PAR', 'state_model', 'score_plifs', 'transition_scores', ...
     'trn_acc', 'val_acc', 'train_exm_ids', 'holdout_exm_ids', 'progress');

% eof
