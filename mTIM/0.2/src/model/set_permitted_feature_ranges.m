function feat_ranges = set_permitted_feature_ranges(CFG)

STATES = get_state_set_mTIM(CFG.PAR);
fn = fieldnames(STATES);
ifw_states = strmatch('IFW', fn)';
ilw_states = strmatch('ILW', fn)';

ifc_states = strmatch('IFC', fn)';
ilc_states = strmatch('ILC', fn)';

FEATS = get_feature_set_mTIM();
LOW_COVER_FEAT = strmatch('low_cover_blocks', FEATS, 'exact');
ACC_W_FEAT = strmatch('acc_pred_W', FEATS, 'exact');
DON_W_FEAT = strmatch('don_pred_W', FEATS, 'exact');
ACC_C_FEAT = strmatch('acc_pred_C', FEATS, 'exact');
DON_C_FEAT = strmatch('don_pred_C', FEATS, 'exact');

min_feat_ranges = -inf*ones(length(fn), length(FEATS));
max_feat_ranges = +inf*ones(length(fn), length(FEATS));

if  CFG.PAR.gene_states_low_cover_cutoff,
  gene_states = setdiff([1:length(fn)], STATES.IGE);
  max_feat_ranges(gene_states, LOW_COVER_FEAT) ...
      = CFG.PAR.gene_states_low_cover_cutoff;
end

if CFG.PAR.enforce_splice_site_consensus,
  min_feat_ranges(ifw_states, DON_W_FEAT) = 0;
  min_feat_ranges(ilw_states, ACC_W_FEAT) = 0;

  min_feat_ranges(ifc_states, ACC_C_FEAT) = 0;
  min_feat_ranges(ilc_states, DON_C_FEAT) = 0;
end


feat_ranges = {min_feat_ranges, max_feat_ranges};
