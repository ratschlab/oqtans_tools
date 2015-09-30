function state_model = specify_model(PAR)

% state_model = specify_model(PAR)
%
% User definition of the state-transition model.
%
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
% returns a preliminary version of the graphical model to be completed by
%   complete_model.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008

%fprintf('Model changed!\n');


%%% define state names and corresponding state ids
STATES = get_state_set_mTIM(PAR);

% define which states are start and stop states (for the purpose of decoding)
fn = fieldnames(STATES);
for i=1:length(fn),
  state_model(i).is_start = 0;
  state_model(i).is_stop  = 0;
end
state_model(STATES.IGE).is_start = 1;
state_model(STATES.IGE).is_stop  = 1;

% some state indices
% Watson strand exons (first, internal, last)
efw_idx = strmatch('EFW', fn);
assert(length(efw_idx) == PAR.num_levels);
eiw_idx = strmatch('EIW', fn);
assert(length(eiw_idx) == PAR.num_levels);
elw_idx = strmatch('ELW', fn);
assert(length(elw_idx) == PAR.num_levels);
% Crick strand exons (first, internal, last)
efc_idx = strmatch('EFC', fn);
assert(length(efc_idx) == PAR.num_levels);
eic_idx = strmatch('EIC', fn);
assert(length(eic_idx) == PAR.num_levels);
elc_idx = strmatch('ELC', fn);
assert(length(elc_idx) == PAR.num_levels);
% Watson strand introns (first, internal, last)
ifw_idx = strmatch('IFW', fn); 
assert(length(ifw_idx) == PAR.num_levels);
iiw_idx = strmatch('IIW', fn); 
assert(length(iiw_idx) == PAR.num_levels);
ilw_idx = strmatch('ILW', fn);
assert(length(ilw_idx) == PAR.num_levels);
% Crick strand introns (first, internal, last)
ifc_idx = strmatch('IFC', fn);
assert(length(ifc_idx) == PAR.num_levels);
iic_idx = strmatch('IIC', fn); 
assert(length(iic_idx) == PAR.num_levels);
ilc_idx = strmatch('ILC', fn);
assert(length(ilc_idx) == PAR.num_levels);
% Watson strand states
w_idx = sort([efw_idx; eiw_idx; elw_idx; ifw_idx; iiw_idx; ilw_idx]);
% Crick strand states
c_idx = sort([efc_idx; eic_idx; elc_idx; ifc_idx; iic_idx; ilc_idx]);

%%% associate a label with each state
LABELS = get_label_set_mTIM();
state_model(STATES.IGE).label   = LABELS.intergenic;
for i=1:PAR.num_levels,
  state_model(efw_idx(i)).label = LABELS.exon_W;
  state_model(eiw_idx(i)).label = LABELS.exon_W;
  state_model(elw_idx(i)).label = LABELS.exon_W;
  state_model(efc_idx(i)).label = LABELS.exon_C;
  state_model(eic_idx(i)).label = LABELS.exon_C;
  state_model(elc_idx(i)).label = LABELS.exon_C;
end
for i=1:PAR.num_levels,
  state_model(ifw_idx(i)).label = LABELS.intron_W;
  state_model(iiw_idx(i)).label = LABELS.intron_W;
  state_model(ilw_idx(i)).label = LABELS.intron_W;
  state_model(ifc_idx(i)).label = LABELS.intron_C;
  state_model(iic_idx(i)).label = LABELS.intron_C;
  state_model(ilc_idx(i)).label = LABELS.intron_C;
end
for s=1:length(state_model),
  assert(~isempty(state_model(s).label));
end
assert(length(fn) == length(state_model));

%%% define allowed transitions
% successors contains all ids of states reachable via an arc from this state
% trans_scores indicates whether a transition score will be learned (if 1)
% or whether transition score will be fixed to 0 (if 0).
% 1s are later replaced by an index into transition_scores.
state_model(STATES.IGE).successors       = [STATES.IGE efw_idx' efc_idx'];
state_model(STATES.IGE).trans_scores     = [1 ones(1,2*PAR.num_levels)];

for i=1:PAR.num_levels,
  state_model(efw_idx(i)).successors     = [eiw_idx(i)];
  state_model(efw_idx(i)).trans_scores   = [1];
  
  state_model(efc_idx(i)).successors     = [eic_idx(i)];
  state_model(efc_idx(i)).trans_scores   = [1];

  state_model(eiw_idx(i)).successors     = [eiw_idx(i) elw_idx(i) ifw_idx(i)];
  state_model(eiw_idx(i)).trans_scores   = [1 1 1];
  
  state_model(eic_idx(i)).successors     = [eic_idx(i) elc_idx(i) ifc_idx(i)];
  state_model(eic_idx(i)).trans_scores   = [1 1 1];

  state_model(elw_idx(i)).successors     = [STATES.IGE];
  state_model(elw_idx(i)).trans_scores   = [1];
  
  state_model(elc_idx(i)).successors     = [STATES.IGE];
  state_model(elc_idx(i)).trans_scores   = [1];


  
  state_model(ifw_idx(i)).successors     = [iiw_idx(i)];
  state_model(ifw_idx(i)).trans_scores   = [1];
  
  state_model(ifc_idx(i)).successors     = [iic_idx(i)];
  state_model(ifc_idx(i)).trans_scores   = [1];
  
  state_model(iiw_idx(i)).successors     = [iiw_idx(i) ilw_idx(i)];
  state_model(iiw_idx(i)).trans_scores   = [1 1];
  
  state_model(iic_idx(i)).successors     = [iic_idx(i) ilc_idx(i)];
  state_model(iic_idx(i)).trans_scores   = [1 1];

  state_model(ilw_idx(i)).successors     = [eiw_idx(i)];
  state_model(ilw_idx(i)).trans_scores   = [1];
  
  state_model(ilc_idx(i)).successors     = [eic_idx(i)];
  state_model(ilc_idx(i)).trans_scores   = [1];
end
% transitions from level i to i+1 and vice versa
for i=1:PAR.num_levels-1,
  state_model(ilw_idx(i)).successors(end+1)     = eiw_idx(i+1);
  state_model(ilw_idx(i)).trans_scores(end+1)   = 1;

  state_model(ilc_idx(i)).successors(end+1)     = eic_idx(i+1);
  state_model(ilc_idx(i)).trans_scores(end+1)   = 1;

  state_model(ilw_idx(i+1)).successors(end+1)   = eiw_idx(i);
  state_model(ilw_idx(i+1)).trans_scores(end+1) = 1;

  state_model(ilc_idx(i+1)).successors(end+1)   = eic_idx(i);
  state_model(ilc_idx(i+1)).trans_scores(end+1) = 1;
end
% make sure that no transition directly connects the models for different strands
for i=1:length(w_idx),
  assert(isempty(intersect(state_model(w_idx(i)).successors, c_idx)));
end
for i=1:length(c_idx),
  assert(isempty(intersect(state_model(c_idx(i)).successors, w_idx)));
end


%%% specify whether feature scoring functions are learned
%%% expected is a 0/1 vector with nonzero entries for the features to be
%%% scored by functions included in the learning process
FEATS = get_feature_set_mTIM();

intron_span = 'intron_span';
if isfield(PAR,'ige_intron_span') && PAR.ige_intron_span==0,
    intron_span = '';
end

state_model(STATES.IGE).learn_scores = ismember(FEATS, {'exon_cover','total_cover', ...
    intron_span,'low_cover_blocks','repeats','cufflinks_feat_W','cufflinks_feat_C','pair_span'});
    

for i=1:PAR.num_levels,
  state_model(efw_idx(i)).learn_scores ...
    = ismember(FEATS, {'cover_grad'});
  state_model(eiw_idx(i)).learn_scores ...
    = ismember(FEATS, {'exon_cover', 'exon_diff','repeats','cufflinks_feat_W','cufflinks_feat_C'});
  state_model(elw_idx(i)).learn_scores ...
    = ismember(FEATS, {'cover_grad'});

  state_model(efc_idx(i)).learn_scores ...
    = ismember(FEATS, {'cover_grad'});
  state_model(eic_idx(i)).learn_scores ...
    = ismember(FEATS, {'exon_cover', 'exon_diff','repeats','cufflinks_feat_W','cufflinks_feat_C'});
  state_model(elc_idx(i)).learn_scores ...
      = ismember(FEATS, {'cover_grad'});

  state_model(ifw_idx(i)).learn_scores ...
    = ismember(FEATS, {'intron_start_W', 'don_pred_W'});
  state_model(iiw_idx(i)).learn_scores ...
    = ismember(FEATS, {'pair_span','intron_span', 'intron_diff','low_cover_blocks','repeats','intron_span_low','intron_span_med','intron_span_high','cufflinks_feat_W','cufflinks_feat_C'});
  state_model(ilw_idx(i)).learn_scores ...
    = ismember(FEATS, {'intron_end_W', 'acc_pred_W'});
  
  state_model(ifc_idx(i)).learn_scores ...
    = ismember(FEATS, {'intron_start_C', 'acc_pred_C'});
  state_model(iic_idx(i)).learn_scores ...
    = ismember(FEATS, {'pair_span','intron_span', 'intron_diff','low_cover_blocks','repeats','intron_span_low','intron_span_med','intron_span_high','cufflinks_feat_W','cufflinks_feat_C'});
  state_model(ilc_idx(i)).learn_scores ...
    = ismember(FEATS, {'intron_end_C', 'don_pred_C'});
end



%%% specify whether scoring functions should be shared between several
%%% states as a matrix k x 2, where k is equal to the number of nonzeros
%%% in learn_scores of the same state
%%% first column is a feature index and second column indicates the state
%%% id  to which the scoring parameters correspond
for i=1:length(state_model),
  lidx = find(state_model(i).learn_scores)';
  state_model(i).feature_scores ...
      = [lidx, i*ones(size(lidx))];
end

% where possible re-use W-strand feature scores for C strand states
for i=1:PAR.num_levels,
  state_model(efc_idx(i)).feature_scores(:,2) = ...
      state_model(efw_idx(i)).feature_scores(:,2);
  % plif for exon diff feature will not be shared
  fidx = 1:length(state_model(eiw_idx(i)).feature_scores);
  fidx(2) = [];
  state_model(eic_idx(i)).feature_scores(fidx,2) = ...
      state_model(eiw_idx(i)).feature_scores(fidx,2);
  state_model(elc_idx(i)).feature_scores(:,2) = ...
      state_model(elw_idx(i)).feature_scores(:,2);

  % plif for intron diff feature will not be shared
  fidx = 1:length(state_model(iiw_idx(i)).feature_scores);
  fidx(2) = [];
  state_model(iic_idx(i)).feature_scores(fidx,2) = ...
      state_model(iiw_idx(i)).feature_scores(fidx,2);
end


%%% specify monotonicity constraints for feature scoring functions
%%% as a vector of length k containing +1 (monotonically increasing
%%% scoring function), -1 (monotonically decreasing) and 0 (no
%%% monotonicity desired) entries where k is equal to the
%%% number of nonzeros in learn_scores of the same state.
%%% will not be considered when scoring functions are shared with
%%% another states
for i=1:length(state_model),
  n = sum(state_model(i).learn_scores);
  state_model(i).monot_scores = zeros(1, n);
end

if (isfield(PAR,'enf_monot_score_funcs_IGE') && PAR.enf_monot_score_funcs_IGE),
  state_model(STATES.IGE).monot_scores   = [-1 -1 0 0 0 0 0];
end

if PAR.enf_monot_score_funcs,
  state_model(STATES.IGE).monot_scores   = [-1 -1 0 0 0 0 0];
  for i=1:PAR.num_levels,
    state_model(efw_idx(i)).monot_scores = [+1];
    state_model(eiw_idx(i)).monot_scores = [ 0 +1 0 0 0];
    state_model(elw_idx(i)).monot_scores = [+1];
    state_model(efc_idx(i)).monot_scores = [+1];
    state_model(eic_idx(i)).monot_scores = [ 0 -1 0 0 0];
    state_model(elc_idx(i)).monot_scores = [+1];
    
    state_model(ifw_idx(i)).monot_scores = [+1 +1];
    state_model(iiw_idx(i)).monot_scores = [ 0 +1 0 0 0 0 0 0];
    state_model(ilw_idx(i)).monot_scores = [+1 +1];
    
    state_model(ifc_idx(i)).monot_scores = [+1 +1];
    state_model(iic_idx(i)).monot_scores = [ 0 -1 0 0 0 0 0 0];
    state_model(ilc_idx(i)).monot_scores = [+1 +1];
  end
end


%%% specify whether feature scoring functions will be coupled via
%%% regularization terms to those of other states as a k x 2 matrix where
%%% k is equal to the number of nonzeros in learn_scores of the same state. 
%%% first column is a feature index and second column indicates the state
%%% id to which the scoring parameters correspond (both should be zero
%%% if no coupling is desired)
%%% AVOID TO COUPLE the same pair of states twice as (i,j) and (j,i).
%%% only feature scoring functions which are not shared between states
%%% can be coupled
for i=1:length(state_model),
  state_model(i).score_coupling ...
      = zeros(sum(state_model(i).learn_scores),2);
end


% coupling between expression levels
for i=1:PAR.num_levels-1,
  % couple all plifs of W strand states to next level
  cidx = find(state_model(efw_idx(i)).learn_scores)';
  state_model(efw_idx(i)).score_coupling = [cidx, efw_idx(i+1)*ones(size(cidx))];
  cidx = find(state_model(eiw_idx(i)).learn_scores)';
  state_model(eiw_idx(i)).score_coupling = [cidx, eiw_idx(i+1)*ones(size(cidx))];
  cidx = find(state_model(elw_idx(i)).learn_scores)';
  state_model(elw_idx(i)).score_coupling = [cidx, elw_idx(i+1)*ones(size(cidx))];

  % for C states, only couple those plifs that are not shared with the W strand states 
  cidx = find(state_model(eic_idx(i)).feature_scores(:,2) == eic_idx(i))';
  state_model(eic_idx(i)).score_coupling(cidx,:) ...
      = [state_model(eic_idx(i)).feature_scores(cidx,1), eic_idx(i+1)*ones(length(cidx),1)];

  % couple all plifs of W strand states to next level
  cidx = find(state_model(ifw_idx(i)).learn_scores)';
  state_model(ifw_idx(i)).score_coupling = [cidx, ifw_idx(i+1)*ones(size(cidx))];
  cidx = find(state_model(iiw_idx(i)).learn_scores)';
  state_model(iiw_idx(i)).score_coupling = [cidx, iiw_idx(i+1)*ones(size(cidx))];
  cidx = find(state_model(ilw_idx(i)).learn_scores)';
  state_model(ilw_idx(i)).score_coupling = [cidx, ilw_idx(i+1)*ones(size(cidx))];

  % for C states, only couple those plifs that are not shared with the W strand states 
  cidx = find(state_model(ifc_idx(i)).feature_scores(:,2) == ifc_idx(i))';
  state_model(ifc_idx(i)).score_coupling(cidx,:) ...
      = [state_model(ifc_idx(i)).feature_scores(cidx,1), ifc_idx(i+1)*ones(length(cidx),1)];
  cidx = find(state_model(iic_idx(i)).feature_scores(:,2) == iic_idx(i))';
  state_model(iic_idx(i)).score_coupling(cidx,:) ...
      = [state_model(iic_idx(i)).feature_scores(cidx,1), iic_idx(i+1)*ones(length(cidx),1)];
  cidx = find(state_model(ilc_idx(i)).feature_scores(:,2) == ilc_idx(i))';
  state_model(ilc_idx(i)).score_coupling(cidx,:) ...
      = [state_model(ilc_idx(i)).feature_scores(cidx,1), ilc_idx(i+1)*ones(length(cidx),1)];
end

% eof
