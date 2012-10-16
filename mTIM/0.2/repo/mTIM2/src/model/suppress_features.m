function [PAR, signal, score_plifs, state_model] = suppress_features(PAR, signal, score_plifs, state_model)
% Turns off feature-groups.
% 
% If PAR.suppressed_feats_map ist found then assume that
%  1. model
%  2. permitted feature ranges
% already have been adjusted and check only 
% if the signal matrix needs adaptation.
%
% IN:
%   PAR         : parameter struct with possible feature-group parameters (e.g. PAR.use_pair_feats=0)
%   signal      : [features x n] observations
%   score_plifs : struct(features x npn) 
%   state_model : struct(#states) with fields (modified field are marked with *)
%                   -trans_scores
%                   -learn_scores*
%                   -feature_scores*
%                   -monot_scores*
%                   -score_coupling*
%                 if state_model==[] then return also [].                   
%
% ABBR:
%   npn = num_plif_nodes (as defined in PAR).
%
% ASSUMPTIONS:
%   1. There is no score_coupling between different feature groups.
%
%
%
% written by nico goernitz, TU Berlin, 2012

is_model_adapted = isfield(PAR,'suppressed_feats_map');
num_features = size(signal,1);
fprintf('Number of features in signal: %i\n', num_features);
fprintf('Number of score_plifs: %i\n', size(score_plifs,1));
fprintf('Is reduced model: %i\n',is_model_adapted);

% assume that all available features are present in the 
% signal and number score_plifs is the same (for newly
% initialized model as in training) or less (for a already
% trained model).
assert(num_features>=size(score_plifs,1));

if (num_features>size(score_plifs,1)),
    fprintf('Suppressing features only in the signal vector.\n');
end

config = model_config();

% check if feature description is available otherwise
% suppression would not work
if (~isfield(config,'func_get_feature_set')),
    fprintf('No feature description found in model.\n');
    exit;
end

% convert to function handle and get the feature description
get_feature_desc = str2func(config.func_get_feature_set);
[FEATS, FEAT_MAP] = get_feature_desc();

% features are by default active!
% reduce the effective number of features
fprintf('\nCheck feature-groups...\n');
remove_feats = [];
for i = 1:size(FEAT_MAP,1),
    par_name = FEAT_MAP{i,1};
    par_feats = FEAT_MAP{i,2};

    fprintf('Feature-group parameter <%s> has <%i> subfeatures.\n',par_name,length(par_feats));
    % if the feature group should not be active..
    if (isfield(PAR,par_name) && getfield(PAR,par_name)==0), 
        fprintf('  Turning OFF the following features:\n');
        for j = 1:length(par_feats),
            dim = par_feats(j);
            fprintf('    -%3i: %s\n',dim,FEATS{dim});    
            remove_feats = [remove_feats, dim];
        end
    end
end

% the remaining feature dimensions
dims = setdiff([1:num_features],remove_feats);

% setup dimension map
map = [[1:num_features];[1:num_features]]';
map(remove_feats,2) = 0;
map(map(:,2)>0,2) = [1:(num_features-length(remove_feats))];
fprintf('Dimension-map:\n');
disp(map)

% assume that the generated map is the same as in PAR
% (if the corresponding field exists).
if ~isfield(PAR,'suppressed_feats_map'),
    PAR.suppressed_feats_map = map;
else
    % check for equity
    assert(size(map,1)==size(PAR.suppressed_feats_map,1));
    assert(size(map,2)==size(PAR.suppressed_feats_map,2));
    assert(sum(sum(map-PAR.suppressed_feats_map,1))==0);
end


% compress signal
signal = signal(dims,:);
PAR.num_features = length(dims);


% change permitted feature ranges
if isfield(PAR, 'perm_feature_ranges') && ~is_model_adapted,
  % if permitted feature ranges are specified those will be accounted for
  % during decoding to restrict the state sequence such that if a feature
  % value for a given state is outside this range, another state will be
  % decoded at this particular position
  fprintf('Found permitted feature ranges map.\n');
  min_feature_ranges = PAR.perm_feature_ranges{1};
  max_feature_ranges = PAR.perm_feature_ranges{2};
  PAR.perm_feature_ranges = {min_feature_ranges(:,dims), max_feature_ranges(:,dims)};
end 

% change score_plifs
if (num_features==size(score_plifs,1)),
    fprintf('INFO: Old score_plifs adaptation (1).\n');
else
    fprintf('INFO: Old score_plifs adaptation (2).\n');
end
    
if (~is_model_adapted),
    fprintf('Adapt score plifs structure.\n');
    score_plifs = score_plifs(dims,:);
end

% change the model
if (~isempty(state_model) && ~is_model_adapted),
    % number of states
    states = length(state_model);
    
    for s = 1:states,

        % FEATURE SCORES
        % - which of the learn_score features is shared with another state
        % - k x 2 matrix; k=number of non-zero features in learn_scores
        % - (i,1) => feature index
        % - (i,2) => state index
        % - 1. remove all remove_features in (:,1) that have non-zero learn_scores
        % - 2. map remaining features (:,1) to new indices
        [foo, inds] = setdiff(state_model(s).feature_scores(:,1),remove_feats);
        state_model(s).feature_scores = state_model(s).feature_scores(inds,:);
        state_model(s).feature_scores(:,1) = map(state_model(s).feature_scores(:,1),2);

        % MONOT SCORES
        % - monotonicity of the features (0=not required, +1 and -1)
        % - k x 1 vector; k=number of non-zero features in learn_scores
        inds = find(state_model(s).learn_scores > 0);
        [foo, keepinds] = setdiff(inds,remove_feats);
        state_model(s).monot_scores = state_model(s).monot_scores(keepinds);

        % SCORE COUPLING
        % - coupling of features with the same feature belonging to another state
        % - k x 2 vector; k=number of non-zero features in learn_scores
        % - (i,1) => feature index
        % - (i,2) => state index
        % ATTENTION: no coupling yields a (0,0)-pair which is removed hereby!!!
        % - 1. remove all remove_features in (:,1) that have non-zero learn_scores
        % - 2. map remaining features (:,1) to new indices

        % remove all (0,0) pairs
        inds = find(state_model(s).score_coupling(:,1)>0);
        state_model(s).score_coupling = state_model(s).score_coupling(inds,:);

        [foo, inds] = setdiff(state_model(s).score_coupling(:,1),remove_feats);
        state_model(s).score_coupling = state_model(s).score_coupling(inds,:);
        state_model(s).score_coupling(:,1) = map(state_model(s).score_coupling(:,1),2);

        % LEARN SCORES
        % - 0/1-vector with length #features
        % - delete all dimensions in remove_feats
        state_model(s).learn_scores = state_model(s).learn_scores(dims);
     end
end




PAR.num_features = num_features-length(remove_feats);
fprintf('Number of (effective) features: %i\n', PAR.num_features);
assert(PAR.num_features==size(score_plifs,1));

% TODO: CHECK CHANGES !
warning('Suppress features is without sanity checks!');
