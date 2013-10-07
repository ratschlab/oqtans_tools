function fhs = view_model(state_model, score_plifs, PAR, transition_scores)

% fhs = view_model(state_model, score_plifs, PAR, [transition_scores])
%
% Plots the score plifs and transition scores of the given model.
%
% state_model -- graphical model (see make_model_mTIM.m)
% score_plifs -- a struct representation of feature scoring functions
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m 
% transition scores -- (optional) scores associated with allowed
%   transition between states (see make_model_mTIM.m)
% returns a vector of file handles of the model figures
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009

fhs = [];

FEATS = get_feature_set_mTIM();
STATES = get_state_set_mTIM(PAR);
fn = fieldnames(STATES);

% visualize intergenic state
ige_idx = strmatch('IGE', fn, 'exact');
assert(length(ige_idx)==1);
f_idx = find(state_model(ige_idx).learn_scores);
for i=1:length(f_idx),
  f = f_idx(i);
  s = ige_idx;
  fhs(end+1) = figure;
  scores = score_plifs(f,s).scores';
  limits = score_plifs(f,s).limits';
  plot(limits, scores, '.-', 'LineWidth', 1.5);
  xlabel('signal');
  ylabel('score');
  title([state_model(s).name ', feature ' num2str(f)]);
  grid on
end

% visualize exon/intron states in groups
submodel_names = {'EFW', 'EIW', 'ELW', 'EFC', 'EIC', 'ELC', 'IFW', 'IFC', 'IIW', 'IIC', 'ILW', 'ILC'};
for i=1:length(submodel_names),
  subm_idx = strmatch(submodel_names{i}, fn);
  assert(length(subm_idx) == PAR.num_levels);
  f_idx = find(state_model(subm_idx(1)).learn_scores);
  for j=1:length(f_idx),
    f = f_idx(j);
    for k=1:length(subm_idx),
      s = subm_idx(k);
      scores(:,k) = score_plifs(f,s).scores';
      limits(:,k) = score_plifs(f,s).limits';
      state_nos(:,k) = repmat(k,PAR.num_plif_nodes,1);
    end
    fhs(end+1) = figure;
    plot3(limits, state_nos, scores, '.-', 'LineWidth', 1.5);
    xlabel('signal');
    ylabel('state');
    zlabel('score');
    title([submodel_names{i} ', ' FEATS{f}]);
    grid on
    colors = color_grad(PAR.num_levels);
    colors = colors(end:-1:1,:);
    ch = get(gca, 'Children');
    assert(length(ch) == PAR.num_levels);
    for c=1:length(ch),
      set(ch(c), 'Color', colors(c,:));
    end
  end
end

% visualize transition scores if given
if exist('transition_scores', 'var'),
  fhs(end+1) = figure;
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
