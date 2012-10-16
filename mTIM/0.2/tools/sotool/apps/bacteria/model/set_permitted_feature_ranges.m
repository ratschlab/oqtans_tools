function feat_ranges = set_permitted_feature_ranges()

STATES = get_state_set();
num_states = length(fieldnames(STATES))

eps = 0.1;

num_feats = 1;
triplet_feat_idx = 1;

start_codon_feat_values = 0:4;
stop_codon_feat_values = 61:63;

min_feat_ranges = -inf*ones(num_states, num_feats);
max_feat_ranges = +inf*ones(num_states, num_feats);

% restrict startCodon state to start codon feature values
min_feat_ranges(STATES.startCodon, triplet_feat_idx) ...
    = min(start_codon_feat_values) - eps;
max_feat_ranges(STATES.startCodon, triplet_feat_idx) ...
    = max(start_codon_feat_values) + eps;

% restrict stopCodon state to stop codon feature values
min_feat_ranges(STATES.stopCodon, triplet_feat_idx) ...
    = min(stop_codon_feat_values) - eps;
max_feat_ranges(STATES.stopCodon, triplet_feat_idx) ...
    = max(stop_codon_feat_values) + eps;

% restrict exon1 state to exclude stop codon feature values; this is only
% necessary for the exonic1 state because it corresponds to in-frame codons
% (whereas exonic2 and exonic3 correspond to frame-shifted 3-grams)
max_feat_ranges(STATES.exonic1, triplet_feat_idx) ...
     = min(stop_codon_feat_values) - 1 + eps;

feat_ranges = {min_feat_ranges, max_feat_ranges};

% eof