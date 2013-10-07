function [accs, dbg] = check_accuracy(PAR, model, trainset, w)
% check_accuracy(trainset)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

% get transition scores and plif scores from solution vector w
[transition_scores, score_plifs] = res_to_scores(w, model.state_model, model.res_map, model.score_plifs, PAR);

% enable the calculation of the empirical loss 
PAR.empirical_loss =1;
PAR.force_predict_only = 1;

for i=1:length(PAR.include_paths),
  addpath(PAR.include_paths{i});
end
  
t_start_ac = clock();
state_seq = [];
pred_state_seq = [];
train_loss = [];
for j=1:length(trainset.train_exm_ids),
  trn_idx = find(trainset.exm_id_intervals(:,1)==trainset.train_exm_ids(j));
  trn_idx = trainset.exm_id_intervals(trn_idx,2):trainset.exm_id_intervals(trn_idx,3);
  trn_obs_seq = trainset.signal(:,trn_idx);
  trn_pred_path = decode_viterbi(trn_obs_seq, model.transition_scores, ...
     model.score_plifs, PAR, trainset.label(trn_idx), trainset.state_label(trn_idx));
  trn_true_label_seq = trainset.label(trn_idx);
  trn_pred_label_seq = trn_pred_path.label_seq;
  trn_acc(j) = mean(trn_true_label_seq(1,:)==trn_pred_label_seq(1,:));
  train_loss = [train_loss, sum(trn_pred_path.loss)];

  % glue all training state labels together
  state_seq = [state_seq, trainset.state_label(trn_idx)];
  % also predicted state labels
  pred_state_seq = [pred_state_seq, trn_pred_path.state_seq];
end
fprintf(['\n  LSL training accuracy:              %2.2f%%\n'], 100*mean(trn_acc));
accs(1)=mean(trn_acc);
dbg.train_emp_loss = sum(train_loss);
dbg.train_state_seq = state_seq;
dbg.train_pred_state_seq = pred_state_seq;


%%% check prediction accuracy on holdout examples
dbg.vald_emp_loss = 0;
dbg.vald_state_seq = [];
dbg.vald_pred_state_seq = [];

if ~isempty(trainset.holdout_exm_ids),
  state_seq = [];
  pred_state_seq = [];
  vald_loss = [];
  for j=1:length(trainset.holdout_exm_ids),
    val_idx = find(trainset.exm_id_intervals(:,1)==trainset.holdout_exm_ids(j));
    val_idx = trainset.exm_id_intervals(val_idx,2):trainset.exm_id_intervals(val_idx,3);
    val_obs_seq = trainset.signal(:,val_idx);
    val_pred_path = decode_viterbi(val_obs_seq, model.transition_scores, ...
        model.score_plifs, PAR, trainset.label(val_idx),trainset.state_label(val_idx));
    val_true_label_seq = trainset.label(val_idx);
    val_pred_label_seq = val_pred_path.label_seq;
    val_acc(j) = mean(val_true_label_seq(1,:)==val_pred_label_seq(1,:));
    vald_loss = [vald_loss, sum(val_pred_path.loss)];
  
    % glue all training state labels together
    state_seq = [state_seq, trainset.state_label(val_idx)];
    % also predicted state labels
    pred_state_seq = [pred_state_seq, val_pred_path.state_seq];
  end
  fprintf(['  LSL validation accuracy:            %2.2f%%\n\n'], ...
          100*mean(val_acc));
  accs(2) = mean(val_acc);
  dbg.vald_emp_loss = sum(vald_loss);
  dbg.vald_state_seq = state_seq;
  dbg.vald_pred_state_seq = pred_state_seq;

end
t_stop_ac = clock();
fprintf('Performance checks took %3.2f sec\n\n', etime(t_stop_ac, t_start_ac));

fprintf('Real empirical loss R(train/vald)=%f/%f\n', ...
     dbg.train_emp_loss,dbg.vald_emp_loss);

% enable the calculation of the empirical loss
PAR.empirical_loss =0;
PAR.force_predict_only = 0;

