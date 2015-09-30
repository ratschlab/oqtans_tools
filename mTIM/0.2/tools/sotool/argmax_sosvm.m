function [PAR, model, w_p, w_n, loss, trainAcc] = argmax_sosvm(PAR, model, w, Xi, Yi)


state_model = model.state_model;
res_map = model.res_map;
LABELS = model.LABELS;

% get transition scores and plif scores from solution vector w
[transition_scores, score_plifs] = res_to_scores(w, state_model, res_map, model.score_plifs, PAR);
model.score_plifs = score_plifs;
model.transition_scores = transition_scores;


% Viterbi decoding
[pred_path true_path pred_path_mmv] = decode_viterbi(...
  Xi.obs_seq, transition_scores, score_plifs, PAR, Yi.true_label_seq, Yi.true_state_seq);

if (PAR.extra_checks),
  wtest = weights_to_vector(pred_path.transition_weights, pred_path.plif_weights, ...
                        state_model, res_map, PAR);
  assert(abs(wtest*res(1:PAR.num_param) - pred_path.score) < PAR.epsilon);

  wtest = weights_to_vector(pred_path_mmv.transition_weights, ...
                        pred_path_mmv.plif_weights, state_model, ...
                        res_map, PAR);
  assert(abs(wtest*res(1:PAR.num_param) - pred_path_mmv.score) < PAR.epsilon);
end

w_p = weights_to_vector(true_path.transition_weights, ...
                        true_path.plif_weights, ...
                        state_model, res_map, PAR);
w_n = weights_to_vector(pred_path_mmv.transition_weights, ...
                        pred_path_mmv.plif_weights, ...
                        state_model, res_map, PAR);

%if (PAR.iter>2 && Xi.idx==1),
%    true_path.foo
%
%    w'*w_p'
%    keyboard
%end
%

if (isfield(LABELS, 'ambiguous')),
  eval_idx = true_path.label_seq ~= LABELS.ambiguous;
else
  eval_idx = logical(ones(size(true_path.label_seq)));      
end
trainAcc = mean(true_path.label_seq(eval_idx) ...
                == pred_path_mmv.label_seq(eval_idx));
              
loss = sum(pred_path_mmv.loss);


% check for weird constellation
if (norm(w_p-w_n)==0 && loss>PAR.epsilon),
    warning('Weird constellation: norm(weight_delta)==0 && loss>0.');
    save('weird.mat','-v7.3','pred_path_mmv','true_path','Yi','Xi');

    label_diff = sum(Yi.true_label_seq~=pred_path_mmv.label_seq);
    state_diff = sum(Yi.true_state_seq~=pred_path_mmv.state_seq);
 
    fprintf('There are %i label differences and %i state differences.\n',label_diff,state_diff);
end

% debug check
if (Xi.idx<10 && PAR.extra_checks),
  h = figure(1);
  subplot(3,3,Xi.idx);
  hold off;
  plot(pred_path_mmv.state_seq,'-r','linewidth',2);
  hold on;
  plot(pred_path.state_seq+5,'-b','linewidth',2)
  plot(true_path.state_seq+10,'-g','linewidth',2)
  hold off
  ylim([min(pred_path.state_seq)-1 max(true_path.state_seq)+11]);
  refresh(1);
end

              
              
