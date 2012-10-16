function [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)

% [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)
%
% Initializes the feature scoring PLiFs and transition scores.
%
% signal -- the feature matrix (sequence of observations) of size m x n
%   where m is equal to the number of features and n the combined length
%   of the training sequences
% label -- (combined) label sequence(s) of total length n
% state_model -- graphical model (see make_model.m)
% PAR -- a struct of parameters specified in setup_hmsvm_training.m and
%   train_hmsvm.m
% returns a struct representation of feature scoring functions
%   (score_plifs) and a vector of transition scores
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2007-2008

LABELS = get_label_set();

% init a score PLiF for each combination of features and states
num_features = size(signal, 1);
for f=1:num_features,
  s = signal(f,:);
  s = sort(s);

  % determine x values for supporting points of PLiFs
  limits = linspace(1, length(s), PAR.num_plif_nodes+1);
  limits = round((limits(1:end-1)+limits(2:end))/2);
  limits = s(limits);

  for s=1:length(state_model),
    score_plifs(f,s).limits = limits;
    score_plifs(f,s).scores = zeros(size(limits));
    score_plifs(f,s).dim = (f-1)*num_features + s;
  end
end

% init scores for transitions specified to have a score in make_model
for i=1:length(state_model),
  assert(length(state_model(i).successors) ...
         == length(state_model(i).trans_scores));
  idx = state_model(i).trans_scores;
  idx(idx==0) = [];
  transition_scores(idx) = 0;
end
transition_scores = transition_scores';

%view_model(state_model, score_plifs, transition_scores);
%keyboard

% eof