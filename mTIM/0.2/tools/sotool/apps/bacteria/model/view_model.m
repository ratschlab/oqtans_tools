function fhs = view_model(state_model, score_plifs, PAR, transition_scores)

% fhs = view_model(state_model, score_plifs, PAR, [transition_scores])
%
% Plots the score plifs and transition scores of the given model.
%
% state_model -- graphical model (see make_model.m)
% score_plifs -- a struct representation of feature scoring functions
% PAR -- a struct of parameters specified in setup_hmsvm_training.m and
%   train_hmsvm.m
% transition scores -- (optional) scores associated with allowed
%   transition between states (see make_model.m)
% returns a vector of file handles of the model figures
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2008

num_features = size(score_plifs,1);
for s=1:length(state_model),
  for t=1:num_features,
    scores(t,:) = score_plifs(t,s).scores;
    limits(t,:) = score_plifs(t,s).limits;
    feats(t,:) = repmat(t, size(score_plifs(t,s).limits));
  end
  figure
  plot3(limits', feats', scores', '.-');
  xlabel('signal');
  ylabel('feature');
  zlabel('score');
  title(state_model(s).name);
  grid on
  fhs(s) = gcf;
  colors = colormap;
  colors = colors(round(linspace(1,64,num_features)),:);
  ch = get(gca, 'Children');
  assert(length(ch) == num_features);
  for c=1:length(ch),
    set(ch(c), 'Color', colors(c,:));
  end
end

if exist('transition_scores', 'var'),
  figure
  A = zeros(length(state_model));
  for i=1:length(state_model),
    sc_idx = state_model(i).trans_scores;
    tr_idx = state_model(i).successors;
    tr_idx(sc_idx==0) = [];
    sc_idx(sc_idx==0) = [];
    A(i, tr_idx) = transition_scores(sc_idx);
  end
  imagesc(A);
end

% eof
