function [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)
% [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)
% initializes the feature scoring PLiFs and transition scores
  
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany

LABELS = get_label_set();

% init a score PLiF for each combination of features and states
num_features = size(signal, 1);
for f=1:num_features,
  s = signal(f,:);
  % filter out uninformative signals (i.e. splice signals for
  % hybridization positions and vice versa)
  s = s(s>0);
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