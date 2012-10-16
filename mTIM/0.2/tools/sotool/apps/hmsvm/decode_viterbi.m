function [pred_path true_path pred_path_mmv] = decode_viterbi(obs_seq, transition_scores, ...
                                                  score_plifs, PAR, true_label_seq, true_state_seq)
% [pred_path true_path pred_path_mmv] 
%    = decode_viterbi(obs_seq, transition_scores, score_plifs,
%    PAR, [true_label_seq], [true_state_seq])
%
% Calls the Viterbi algorithm implemented in the Shogun toolbox
% (http://www.shogun-toolbox.org/) to decode the best state sequence
% under the current parameters (pred_path) and if a true label and state
% sequence are given also the maximal margin violator (pred_path_mmv).
%
% obs_seq -- sequence of observations, i.e. the feature matrix
% transition_scores -- scores associated with allowed transitions between
%   states
% score_plifs -- a struct representation of feature scoring functions
%   (see also score_plif_struct.h / .cpp)
% PAR -- a struct to configure the HM-SVM (for specification see
%   setup_hmsvm_training.m and train_hmsvm.m)
% true_label_seq -- optional parameter indicating the true label sequence
%   if given, also true_state_seq has to be specified. In this case used
%   transitions and plif weights are computed for the true path,
%   pred_path is augmented with a loss and a struct pred_path_mmv is
%   returned corresponding to the maximal margin violator under the given
%   loss.
% true_state_seq -- true sequence of states, has to be specified iff
%   true_label_seq is given 
% returns a struct representing the decoded state sequence (pred_path),
%   and if true_label_seq is given, also a struct representing the true
%   state sequence and a struct for the max-margin violator
%
% see also compute_score_matrix.cpp
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008-2011
% updated by Nico Goernitz, TU Berlin, 2012:
% - moved the pred_path prediction into the 'else'-block practically
%   avoiding the double call for viterbi decoding during training

[state_model, A] = eval(sprintf('%s(PAR, transition_scores);', ...
                                PAR.model_config.func_make_model));

% Viterbi decoding to obtain best prediction WITHOUT loss

% compute score matrix for decoding
num_states = length(state_model);
STATES = eval(sprintf('%s(PAR);', PAR.model_config.func_get_state_set));

if isfield(PAR, 'perm_feature_ranges'),
  % if permitted feature ranges are specified those will be accounted for
  % during decoding to restrict the state sequence such that if a feature
  % value for a given state is outside this range, another state will be
  % decoded at this particular position
  score_matrix = compute_score_matrix(obs_seq, score_plifs, ...
                PAR.perm_feature_ranges{1}, PAR.perm_feature_ranges{2});
else
  score_matrix = compute_score_matrix(obs_seq, score_plifs);
end
p = -inf(1, num_states);
p(find([state_model.is_start])) = 0;
q = -inf(1, num_states);
q(find([state_model.is_stop]))  = 0;

% define returns
pred_path = [];
pred_path_mmv = [];


% if true_label_seq is given (for training examples),
% true_state_seq has to be specified as well.
% In this case used transitions and plif weights are computed
% for the true path, pred_path is augmented with a loss
% and a struct pred_path_mmv is returned corresponding 
% to the maximal margin violator under the given loss

if exist('true_label_seq', 'var'),
  assert(length(true_label_seq)==size(obs_seq,2));
  assert(all(size(true_state_seq)==size(true_label_seq)));
  
  % score, transition and plif weights for the true path 
  true_path.score = path_score(true_state_seq, score_matrix, A);
  true_path.state_seq = true_state_seq;
  true_path.label_seq = true_label_seq;

  [true_path.transition_weights, true_path.plif_weights] ...
      = path_weights(true_path.state_seq, obs_seq, score_plifs, length(state_model));
 
  % position-wise loss of the decoded state sequence
  loss = feval(PAR.model_config.func_calc_loss_matrix, true_path.state_seq, state_model, PAR);

  % ambiguously labeled parts should not receive losses
  LABEL = feval(PAR.model_config.func_get_label_set);
  if (isfield(LABEL,'ambiguous')),
      inds = find(true_path.label_seq == LABEL.ambiguous);
      % start with the most-left position
      inds = sort(inds,'ascend');
      % not optimal for large amounts of ambigiuous positions
      for i=1:length(inds),
        if (inds(i)>1), loss(:,inds(i))=loss(:,inds(i)-1);
        else loss(:,1)=zeros(size(loss,1),1); end
      end
  end

  % NOT NECCESSARY anymore 
  %pred_loss = zeros(size(pred_state_seq));
  %for i=1:size(pred_state_seq,2),
  %  pred_loss(i) = loss(pred_state_seq(i), i);
  %end
  %pred_path.loss = pred_loss;
  %[pred_path.transition_weights, pred_path.plif_weights] ...
  %    = path_weights(pred_state_seq, obs_seq, score_plifs, length(state_model));

  % Viterbi decoding to obtain best prediction WITH loss, 
  % i.e. the maximal margin violater (MMV)
  
  % add loss to score matrix
  score_matrix = score_matrix + loss;
  
  [pred_path_mmv.score, pred_state_seq] = best_path(p, q, A, score_matrix);
  pred_path_mmv.state_seq = pred_state_seq;
  pred_path_mmv.label_seq = eval(sprintf('%s(pred_state_seq, state_model);', ...
                                         PAR.model_config.func_states_to_labels));
  
  % position-wise loss of the decoded state sequence
  pred_loss = zeros(1, size(obs_seq,2));
  for i=1:size(obs_seq,2),
    pred_loss(i) = loss(pred_state_seq(i), i);
  end
  pred_path_mmv.loss = pred_loss;
  pred_path_mmv.score = pred_path_mmv.score - sum(pred_loss);
  
  [pred_path_mmv.transition_weights, pred_path_mmv.plif_weights] ...
      = path_weights(pred_state_seq, obs_seq, score_plifs, length(state_model));
else

  % if not training
  [pred_path.score, pred_state_seq] = best_path(p, q, A, score_matrix);
  pred_path.state_seq = pred_state_seq;
  pred_path.label_seq = eval(sprintf('%s(pred_state_seq, state_model);', ...
                                   PAR.model_config.func_states_to_labels));
  if PAR.extra_checks,
    if pred_path.score <= 20,  scale = 1; 
    else  scale = round(log10(pred_path.score)); end
    assert(abs(pred_path.score-path_score(pred_state_seq,score_matrix,A))<10^scale*PAR.epsilon);
  end

end
