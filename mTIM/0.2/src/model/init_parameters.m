function [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)

% [score_plifs transition_scores] = init_parameters(signal, label, state_model, PAR)
%
% Initializes the feature scoring PLiFs and transition scores.
%
% signal -- the feature matrix (sequence of observations) of size m x n
%   where m is equal to the number of features and n the combined length
%   of the training sequences
% label -- (combined) label sequence(s) of total length n
% state_model -- graphical model (see make_model_mTIM.m)
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
% returns a struct representation of feature scoring functions
%   (score_plifs) and a vector of transition scores
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008-2009

LABELS = get_label_set_mTIM();
STATES = get_state_set_mTIM(PAR);
[FEATS, FEAT_MAP]  = get_feature_set_mTIM();

% number of all features available
num_features = size(signal, 1);
% description an features should match
assert(length(FEATS)==num_features);

% init a score PLiF for each combination of features and states
fprintf('Number of (all) features: %i\n', num_features);

% PLiF supporting points for the read coverage features should be spaced
% such that roughly equally many data points lie between each pair of
% adjacent supporting points


%%% exon coverage feature (combined for both strands)
F = strmatch('exon_cover', FEATS, 'exact');
assert(F < 6);
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
exo_idx = label==LABELS.exon_W | label==LABELS.exon_C;
active_limits = find_limits(PAR.num_plif_nodes, signal(F,exo_idx));
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of exon PLiFs using the "active limits"
% that were computed on exons only
fn = fieldnames(STATES);
exon_states = strmatch('E', fn);
for i=1:length(exon_states),
  s = exon_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end
fprintf('Exon coverage limits:\n');
disp(general_limits);
disp(active_limits);

%%% exon difference feature
F = strmatch('exon_diff', FEATS, 'exact');
assert(F < 6);
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
exo_idx = label==LABELS.exon_W | label==LABELS.exon_C;
active_limits = find_limits(PAR.num_plif_nodes, signal(F,exo_idx));
%fprintf('Exon diff limits:\n');
%disp(general_limits);
%disp(active_limits);
% make PLiF limits symmetrical about zero
for i=1:floor(PAR.num_plif_nodes/2),
  mn = round((-general_limits(i) + general_limits(end-i+1)) / 2);
  general_limits(i) = -mn;
  general_limits(end-i+1) = mn;

  mn = round((-active_limits(i) + active_limits(end-i+1)) / 2);
  active_limits(i) = -mn;
  active_limits(end-i+1) = mn;
end
fprintf('Symmetrical exon diff limits:\n');
disp(general_limits);
disp(active_limits);
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of exon PLiFs using the "active limits"
% that were computed on exons only
fn = fieldnames(STATES);
exon_states = strmatch('E', fn);
for i=1:length(exon_states),
  s = exon_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end

%%% coverage gradient features
F = strmatch('cover_grad', FEATS, 'exact');
assert(F < 6);
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
exo_end_idx = find_blocks(label==LABELS.exon_W | label==LABELS.exon_C);
ridx = randperm(length(label));
exo_end_idx = sort([exo_end_idx(:)', ridx(1:10*length(exo_end_idx))]);
active_limits = find_limits(PAR.num_plif_nodes, signal(F,exo_end_idx));
%fprintf('Coverage gradient limits:\n');
%disp(general_limits);
%disp(active_limits);
% make PLiF limits symmetrical about zero
for i=1:floor(PAR.num_plif_nodes/2),
  mn = (-general_limits(i) + general_limits(end-i+1)) / 2;
  general_limits(i) = -mn;
  general_limits(end-i+1) = mn;

  mn = (-active_limits(i) + active_limits(end-i+1)) / 2;
  active_limits(i) = -mn;
  active_limits(end-i+1) = mn;
end
fprintf('Symmetrical coverage gradient limits:\n');
disp(general_limits);
disp(active_limits);
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of exon PLiFs using the "active limits"
% that were computed on exons only
fn = fieldnames(STATES);
exon_end_states = sort([strmatch('EF', fn), strmatch('EL', fn)]);
for i=1:length(exon_states),
  s = exon_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end

%%% intron coverage features
F = strmatch('intron_span', FEATS, 'exact');
assert(F < 6);
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
ino_idx = label==LABELS.intron_W | label==LABELS.intron_C;
active_limits = find_limits(PAR.num_plif_nodes, signal(F,ino_idx));
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of intron PLiFs using the "active limits"
% that were computed on introns only
fn = fieldnames(STATES);
intron_states = [strmatch('IF', fn); strmatch('II', fn); strmatch('IL', fn)];
for i=1:length(intron_states),
  s = intron_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end
fprintf('Intron span limits:\n');
disp(general_limits);
disp(active_limits);

%%% intron difference feature
F = strmatch('intron_diff', FEATS, 'exact');
assert(F < 6);
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
ino_idx = label==LABELS.intron_W | label==LABELS.intron_C;
active_limits = find_limits(PAR.num_plif_nodes, signal(F,ino_idx));
%fprintf('Intron diff limits:\n');
%disp(general_limits);
%disp(active_limits);
% make PLiF limits symmetrical about zero
for i=1:floor(PAR.num_plif_nodes/2),
  mn = round((-general_limits(i) + general_limits(end-i+1)) / 2);
  general_limits(i) = -mn;
  general_limits(end-i+1) = mn;

  mn = round((-active_limits(i) + active_limits(end-i+1)) / 2);
  active_limits(i) = -mn;
  active_limits(end-i+1) = mn;
end
fprintf('Symmetrical intron diff limits:\n');
disp(general_limits);
disp(active_limits);
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of intron PLiFs using the "active limits"
% that were computed on introns only
fn = fieldnames(STATES);
intron_states = [strmatch('IF', fn); strmatch('II', fn); strmatch('IL', fn)];
for i=1:length(intron_states),
  s = intron_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end

%%% intron start/end features (based on reads)
limits = find_limits(PAR.num_plif_nodes, [signal(6,:) signal(7,:)]);
% we can ignore signal(8,:) and signal(9,:) because these score
% distributions are identical to the intron start distribution used
for f=6:9,
  for s=1:length(state_model),
    score_plifs(f,s).limits = limits;
    score_plifs(f,s).scores = zeros(size(limits));
    score_plifs(f,s).dim = (f-1)*length(state_model) + s;
  end
end
fprintf('Intron read score limits:\n');
disp(limits);


%%% splice site prediction features 
% have PLiFs with uniformly spaced supporting points
limits = linspace(0, 1, PAR.num_plif_nodes);
for f=10:13,
  for s=1:length(state_model),
    score_plifs(f,s).limits = limits;
    score_plifs(f,s).scores = zeros(size(limits));
    score_plifs(f,s).dim = (f-1)*length(state_model) + s;
  end
end
fprintf('Splice site prediction limits:\n');
disp(limits);

%%% mate pair features (based on reads)
F = strmatch('pair_span', FEATS, 'exact');
assert(F > 13 & F < 17);
limits = find_limits(PAR.num_plif_nodes, signal(F,:));
for s=1:length(state_model),
  score_plifs(F,s).limits = limits;
  score_plifs(F,s).scores = zeros(size(limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
fprintf('Mate pair span limits:\n');
disp(limits);

F = strmatch('total_cover', FEATS, 'exact');
assert(F > 13 & F < 17);
limits = find_limits(PAR.num_plif_nodes, signal(F,:));
fprintf('Total coverage limits:\n');
disp(limits);
for s=1:length(state_model),
  score_plifs(F,s).limits = limits;
  score_plifs(F,s).scores = zeros(size(limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end

F = strmatch('low_cover_blocks', FEATS, 'exact');
assert(F > 13 & F < 17);
limits = find_limits(PAR.num_plif_nodes, signal(F,:));
fprintf('Low-coverage block length limits:\n');
disp(limits);
for s=1:length(state_model),
  score_plifs(F,s).limits = limits;
  score_plifs(F,s).scores = zeros(size(limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end

% TODO: ideally this should depend on whether IGE state is enforced in
% low-coverage blocks
%  if PAR.switch_features,
    % shift PLIF nodes to the range where learning can still take place
limits = linspace(0, PAR.gene_states_low_cover_cutoff, PAR.num_plif_nodes);
%limits = linspace(0, 10, PAR.num_plif_nodes)
score_plifs(F,1).limits = limits;
score_plifs(F,1).scores = zeros(size(limits));


%%% repeat feature
F = strmatch('repeats', FEATS, 'exact');
limits = find_limits(PAR.num_plif_nodes, signal(F,:));
for s=1:length(state_model),
    score_plifs(F,s).limits = limits;
    score_plifs(F,s).scores = zeros(size(limits));
    score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
fprintf('Repeat limits:\n');
disp(limits);

%% binned intron coverage features
for F=18:20,
general_limits = find_limits(PAR.num_plif_nodes, signal(F,:));
ino_idx = label==LABELS.intron_W | label==LABELS.intron_C;
active_limits = find_limits(PAR.num_plif_nodes, signal(F,ino_idx));
% use general limits for spacing of supporting points by default
for s=1:length(state_model),
  score_plifs(F,s).limits = general_limits;
  score_plifs(F,s).scores = zeros(size(general_limits));
  score_plifs(F,s).dim = (F-1)*length(state_model) + s;
end
% overwrite supporting points of intron PLiFs using the "active limits"
% that were computed on introns only
fn = fieldnames(STATES);
intron_states = [strmatch('IF', fn); strmatch('II', fn); strmatch('IL', fn)];
for i=1:length(intron_states),
  s = intron_states(i);
  score_plifs(F,s).limits = active_limits;
  score_plifs(F,s).scores = zeros(size(active_limits));
end
fprintf('Binned intron span %i limits:\n',F);
disp(general_limits);
disp(active_limits);
end

%%% repeat feature
assert(PAR.num_plif_nodes>4);
for F=21:22,
    limits = linspace(-0.5,2.5,PAR.num_plif_nodes);
    for s=1:length(state_model),
        score_plifs(F,s).limits = limits;
        score_plifs(F,s).scores = zeros(size(limits));
        score_plifs(F,s).dim = (F-1)*length(state_model) + s;
    end
    fprintf('Cufflinks %i limits:\n',F);
    disp(limits);
end


assert(size(score_plifs,1) == num_features);
assert(size(score_plifs,2) == length(state_model));
assert(isequal(sort([score_plifs.dim]), unique([score_plifs.dim])));

% init scores for transitions specified to have a score in make_model
for i=1:length(state_model),
  assert(length(state_model(i).successors) ...
         == length(state_model(i).trans_scores));
  idx = state_model(i).trans_scores;
  idx(idx==0) = [];
  transition_scores(idx) = 0;
end
transition_scores = transition_scores';

%view_model(state_model, score_plifs, PAR, transition_scores);
%keyboard


% eof
