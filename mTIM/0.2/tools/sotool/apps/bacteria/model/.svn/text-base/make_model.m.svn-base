function [state_model, A, a_trans] = make_model(PAR, transition_scores)
% [state_model, A, a_trans] = make_model(PAR, [transition_scores])
%
% Combines user specification and automatic completion of the
% state-transition model.
% 
% PAR -- a struct of parameters specified in setup_hmsvm_training.m and
%   train_hmsvm.m
% transition_scores -- optional parameter; if specified, transition scores
%   of A and a_trans will be set accordingly (otherwise initialized to 0)
% returns the fully specified graphical model (state_model), the transition
%   matrix (A) and its sparse representation used for Viterbi decoding in
%   Shogun (a_trans)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2008

% define state names and corresponding state ids
STATES = get_state_set(PAR);

% initialize names and ids in state_model struct
fn = fieldnames(STATES);
for i=1:length(fn),
  state_model(i).id = i;
  state_model(i).name = fn{i};
  state_model(i).is_start = 0;
  state_model(i).is_stop  = 0;
end
state_model(STATES.start).is_start = 1;
state_model(STATES.stop).is_stop  = 1;


% associate a label with each state
LABELS = get_label_set();
state_model(STATES.start).label      = LABELS.intergenic;
state_model(STATES.stop).label       = LABELS.intergenic;
state_model(STATES.intergenic).label = LABELS.intergenic;
state_model(STATES.startCodon).label = LABELS.startCodon;
state_model(STATES.stopCodon).label  = LABELS.stopCodon;
state_model(STATES.exonic1).label    = LABELS.exonic;
state_model(STATES.exonic2).label    = LABELS.exonic;
state_model(STATES.exonic3).label    = LABELS.exonic;


% define allowed transitions:
% successors contains all ids of states reachable via an arc from this state
% trans_scores indicates whether a transition score will be learned (if 1)
% or whether transition score will be fixed to 0 (if 0).
% 1s are later replaced by an index into transition_scores.
state_model(STATES.stop).successors = [];
state_model(STATES.stop).trans_scores = [];
state_model(STATES.start).successors = [STATES.intergenic];
state_model(STATES.start).trans_scores = [0];

state_model(STATES.intergenic).successors = [STATES.intergenic STATES.startCodon STATES.stop];
state_model(STATES.intergenic).trans_scores = [1 1 0];

% start and stop codon
state_model(STATES.startCodon).successors = [STATES.exonic2];
state_model(STATES.startCodon).trans_scores = [1];
state_model(STATES.stopCodon).successors = [STATES.intergenic];
state_model(STATES.stopCodon).trans_scores = [1];

state_model(STATES.exonic1).successors = [STATES.exonic2];
state_model(STATES.exonic1).trans_scores = [1];
state_model(STATES.exonic2).successors = [STATES.exonic3];
state_model(STATES.exonic2).trans_scores = [1];
state_model(STATES.exonic3).successors = [STATES.exonic1 STATES.stopCodon];
state_model(STATES.exonic3).trans_scores = [1 1];



% initialize transition matrices and mapping to transition scores
if ~exist('transition_scores', 'var')
  transition_scores = zeros(1,length(state_model).^2);
end

A = -inf(length(state_model));
q = 1;
for i=1:length(state_model),
  assert(length(state_model(i).successors) == length(state_model(i).trans_scores));
  
  A(i, state_model(i).successors) = 0;
  
  idx_1 = find(state_model(i).trans_scores~=0);
  idx_2 = state_model(i).successors(idx_1);
  
  for j=1:length(idx_1),
    % assign score indices
    cslist_temp = state_model(i).trans_scores;
    cslist_temp(idx_1(j)) = q;
    state_model(i).trans_scores = cslist_temp;

    %state_model(i).trans_scores(idx_1(j)) = q;
    A(i, idx_2(j)) = transition_scores(q);
    q = q + 1;
  end
end


% convert transition matrix to shogun format (a_trans)
a_trans = zeros(3,sum(~isinf(A(:))));
k = 0;
for i=1:size(A,1),
  idx = find(~isinf(A(i,:)));
  val = A(i,idx);
  a_trans(1, k+1:k+length(idx)) = i-1;
  a_trans(2, k+1:k+length(idx)) = idx-1;
  a_trans(3, k+1:k+length(idx)) = val;
  k = k + length(idx);
end
a_trans = a_trans';
[tmp, idx] = sort(a_trans(:,2));
a_trans = a_trans(idx,:);


% specify whether feature scoring functions are learned
% expected is a 0/1 vector with nonzero entries for the features to be
% scored by functions included in the learning process
state_model(STATES.start).learn_scores        = zeros(PAR.num_features,1);
state_model(STATES.stop).learn_scores         = zeros(PAR.num_features,1);

state_model(STATES.intergenic).learn_scores   = ones(PAR.num_features,1);
state_model(STATES.startCodon).learn_scores   = ones(PAR.num_features,1);
state_model(STATES.stopCodon).learn_scores    = ones(PAR.num_features,1);
state_model(STATES.exonic1).learn_scores      = ones(PAR.num_features,1);
state_model(STATES.exonic2).learn_scores      = ones(PAR.num_features,1);
state_model(STATES.exonic3).learn_scores      = ones(PAR.num_features,1);


% specify whether scoring functions should be shared between several
% states as a matrix k x 2, where k is equal to the number of nonzeros
% in learn_scores of the same state
% first column is a feature index and second column indicates the state
% id  to which the scoring parameters correspond
state_model(STATES.start).feature_scores = zeros(sum(state_model(STATES.start).learn_scores),2);
state_model(STATES.stop).feature_scores = zeros(sum(state_model(STATES.stop).learn_scores),2);

state_model(STATES.intergenic).feature_scores = [[1:PAR.num_features]', STATES.intergenic*ones(PAR.num_features,1)];
state_model(STATES.startCodon).feature_scores = [[1:PAR.num_features]', STATES.startCodon*ones(PAR.num_features,1)];
state_model(STATES.stopCodon).feature_scores  = [[1:PAR.num_features]', STATES.stopCodon*ones(PAR.num_features,1)];
state_model(STATES.exonic1).feature_scores    = [[1:PAR.num_features]', STATES.exonic1*ones(PAR.num_features,1)];
state_model(STATES.exonic2).feature_scores    = [[1:PAR.num_features]', STATES.exonic2*ones(PAR.num_features,1)];
state_model(STATES.exonic3).feature_scores    = [[1:PAR.num_features]', STATES.exonic3*ones(PAR.num_features,1)];

for i=1:length(state_model),
  assert(size(state_model(i).feature_scores,1) == sum(state_model(i).learn_scores));
end

% specify monotonicity constraints for feature scoring functions
% as a vector of length k containing +1 (monotonically increasing
% scoring function), -1 (monotonically decreasing) and 0 (no
% monotonicity desired) entries where k is equal to the
% number of nonzeros in learn_scores of the same state.
% will not be considered when scoring functions are shared with
% another states
state_model(STATES.start).monot_scores = zeros(sum(state_model(STATES.start).learn_scores),1);
state_model(STATES.stop).monot_scores = zeros(sum(state_model(STATES.stop).learn_scores),1);

state_model(STATES.intergenic).monot_scores = zeros(sum(state_model(STATES.intergenic).learn_scores),1);
state_model(STATES.startCodon).monot_scores = zeros(sum(state_model(STATES.startCodon).learn_scores),1);
state_model(STATES.stopCodon).monot_scores  = zeros(sum(state_model(STATES.stopCodon).learn_scores),1);
state_model(STATES.exonic1).monot_scores    = zeros(sum(state_model(STATES.exonic1).learn_scores),1);
state_model(STATES.exonic2).monot_scores    = zeros(sum(state_model(STATES.exonic2).learn_scores),1);
state_model(STATES.exonic3).monot_scores    = zeros(sum(state_model(STATES.exonic3).learn_scores),1);

for i=1:length(state_model),
  assert(size(state_model(i).feature_scores,1) == sum(state_model(i).learn_scores));
end


% specify whether feature scoring functions will be coupled via
% regularization terms to those of other states as a k x 2 matrix where
% k is equal to the number of nonzeros in learn_scores of the same state. 
% first column is a feature index and second column indicates the state
% id  to which the scoring parameters correspond (both should be zero
% if no coupling is desired)
% AVOID TO COUPLE the same pair of states twice as (i,j) and (j,i).
% only feature scoring functions which are not shared between states
% can be coupled
state_model(STATES.start).score_coupling = zeros(sum(state_model(STATES.start).learn_scores),2);
state_model(STATES.stop).score_coupling = zeros(sum(state_model(STATES.stop).learn_scores),2);

state_model(STATES.intergenic).score_coupling = zeros(sum(state_model(STATES.intergenic).learn_scores),2);
state_model(STATES.startCodon).score_coupling = zeros(sum(state_model(STATES.startCodon).learn_scores),2);
state_model(STATES.stopCodon).score_coupling  = zeros(sum(state_model(STATES.stopCodon).learn_scores),2);
state_model(STATES.exonic1).score_coupling    = zeros(sum(state_model(STATES.exonic1).learn_scores),2);
state_model(STATES.exonic2).score_coupling    = zeros(sum(state_model(STATES.exonic2).learn_scores),2);
state_model(STATES.exonic3).score_coupling    = zeros(sum(state_model(STATES.exonic3).learn_scores),2);


for i=1:length(state_model),
  assert(all(size(state_model(i).score_coupling) == size(state_model(i).feature_scores)));
  for j=1:size(state_model(i).score_coupling,1),
    f_idx = state_model(i).score_coupling(j,1);
    s_idx = state_model(i).score_coupling(j,2);
    if f_idx ~= 0,
      assert(s_idx ~= 0);
      assert(state_model(s_idx).feature_scores(j,1) == f_idx);
      assert(state_model(s_idx).feature_scores(j,2) == s_idx);
    end
  end
end

% eof
