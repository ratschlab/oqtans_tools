function [transition_scores, score_plifs] = res_to_scores(res, state_model, res_map, score_plifs, PAR)
% [transition_scores, score_plifs] 
%   = res_to_scores(res, state_model, res_map, score_plifs, PAR)
%
% Updates feature scoring functions from the parameter vector (typically
% obtained as a solution of the intermediate QP/LP).
%
% res -- solution vector of the training problem
% state_model -- graphical model specifying states and allowed
%   transitions between them
% res_map -- a mapping between components of the solution vector of the
%   training problem and the score_plif struct
% score_plifs -- a struct representation of feature scoring functions
%   (see also score_plif_struct.h / .cpp)
% PAR -- a struct to configure the HM-SVM (for specification see
%   setup_hmsvm_training.m and train_hmsvm.m)
% returns a vector of transition scores and a struct representation of
%   the feature scoring functions (score_plifs)
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008

assert(PAR.num_features == size(score_plifs,1));
assert(length(state_model) == size(score_plifs,2));
assert(PAR.num_features == size(res_map,1));
assert(length(state_model) == size(res_map,2));

transition_scores = res(1:PAR.num_trans_score);

% extract new feature scoring function values from the solution vector
% (this is only done if specified as learning parameters in make_model.m)
for i=1:size(res_map,1), % for all features
  for j=1:size(res_map,2), % for all states
    if res_map(i,j) ~= 0,
      idx = res_map(i,j):res_map(i,j)+PAR.num_plif_nodes-1; 
      score_plifs(i,j).scores = res(idx)';
    end
  end
end
