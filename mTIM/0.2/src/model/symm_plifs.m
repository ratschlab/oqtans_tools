function [score_plifs transition_scores] = symm_plifs(score_plifs, transition_scores, PAR)

%

warning('symm_plifs is NOT prepared for turning feature-groups on and off.');

STATES = get_state_set_mTIM(PAR);
FEATS = get_feature_set_mTIM();
symm_feats = find(ismember(FEATS, {'exon_cover', ... 
                    'cover_grad', ...
                    'intron_span', ...
                    'pair_span', ...
                    'total_cover', ...
                    'low_cover_blocks', ...
                    'repeats'}));

symm_feat_pairs = [[symm_feats', symm_feats']; ...
                   find(ismember(FEATS, {'intron_start_W', 'intron_start_C'})); ...
                   find(ismember(FEATS, {'intron_end_W', 'intron_end_C'})); ...
                   find(ismember(FEATS, {'acc_pred_W', 'don_pred_C'})); ...
                   find(ismember(FEATS, {'don_pred_W', 'acc_pred_C'})) ...
                  ];

antisymm_feats = find(ismember(FEATS, {'exon_diff', 'intron_diff'}));

state_model = make_model_mTIM(PAR);

state_pairs = zeros(0, 2);
for l=1:PAR.num_levels,
  % TODO

  eiw_state = getfield(STATES, sprintf('EIW_%02i', l));
  eic_state = getfield(STATES, sprintf('EIC_%02i', l));
  state_pairs = [state_pairs; [eiw_state, eic_state]];

  ifw_state = getfield(STATES, sprintf('IFW_%02i', l));
  ifc_state = getfield(STATES, sprintf('IFC_%02i', l));
  state_pairs = [state_pairs; [ifw_state, ifc_state]];

  iiw_state = getfield(STATES, sprintf('IIW_%02i', l));
  iic_state = getfield(STATES, sprintf('IIC_%02i', l));
  state_pairs = [state_pairs; [iiw_state, iic_state]];

  ilw_state = getfield(STATES, sprintf('ILW_%02i', l));
  ilc_state = getfield(STATES, sprintf('ILC_%02i', l));
  state_pairs = [state_pairs; [ilw_state, ilc_state]];
end

for p=1:size(state_pairs,1),
  i = state_pairs(p,1);
  j = state_pairs(p,2);
  state_i = state_model(i);
  state_j = state_model(j);
  
%  fprintf('Symmetrizing %s and %s\n', state_i.name, state_j.name);

  % make transition scores symmetrical
  t = (transition_scores(state_i.trans_scores) + ...
       transition_scores(state_j.trans_scores)) / 2;
  transition_scores(state_i.trans_scores) = t;
  transition_scores(state_j.trans_scores) = t;

  % make PLiFs symmetrical
  for k=1:size(symm_feat_pairs,1),
    f1 = symm_feat_pairs(k,1);
    f2 = symm_feat_pairs(k,2);
    if any(state_i.feature_scores(:,1) == f1);
      assert(any(state_j.feature_scores(:,1) == f2));

      assert(isequal(score_plifs(f1,i).limits, score_plifs(f2,j).limits));
      sc = mean([score_plifs(f1,i).scores; score_plifs(f2,j).scores]);
      score_plifs(f1,i).scores = sc;
      score_plifs(f2,j).scores = sc;
%      fprintf('  + features %s and %s\n', FEATS{f1}, FEATS{f2});
    end
  end

  % make PLiFs anti-symmetrical
  for k=1:length(antisymm_feats),
    f = antisymm_feats(k);
    if any(state_i.feature_scores(:,1) == f);
      assert(any(state_j.feature_scores(:,1) == f));

      assert(isequal(score_plifs(f,i).limits, score_plifs(f,j).limits));
      assert(isequal(score_plifs(f,i).limits, -score_plifs(f,j).limits(end:-1:1)));
      sc = mean([score_plifs(f,i).scores; score_plifs(f,j).scores(end:-1:1)]);
      score_plifs(f,i).scores = sc;
      score_plifs(f,j).scores(end:-1:1) = sc;
%      fprintf('  - feature %s \n', FEATS{f});
    end
  end
end


%view_model(state_model, score_plifs, PAR, transition_scores);
%keyboard
