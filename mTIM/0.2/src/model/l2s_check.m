function [state_seq problems] = l2s_check(label_seq, state_model, signal, PAR)

% debug
% convert ends of the label sequence to ige labels
%
label_seq(1:5) = 0;
label_seq(end:-1:end-5) = 0;


% [state_seq problems] = labels_to_states(label_seq, state_model, signal, PAR)
%
% Converts a label sequence into a state sequence.
%
% label_seq -- a sequence of labels (see get_label_set_mTIM.m)
% state_model -- graphical model (see make_model_mTIM.m)
% signal -- the feature matrix (sequence of observations) of size m x n
%   where m is equal to the number of features and n the combined length
%   of the training sequences
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
% returns the sequence of states corresponding to the given label
%   sequence, second return value indicates whether conversion was problematic 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2011

LABELS = get_label_set_mTIM();
STATES = get_state_set_mTIM(PAR);
fn = fieldnames(STATES);

ef_states = sort([strmatch('EFW', fn)', strmatch('EFC', fn)']);
ei_states = sort([strmatch('EIW', fn)', strmatch('EIC', fn)']);
el_states = sort([strmatch('ELW', fn)', strmatch('ELC', fn)']);

efw_states = strmatch('EFW', fn)';
eiw_states = strmatch('EIW', fn)';
elw_states = strmatch('ELW', fn)';
efc_states = strmatch('EFC', fn)';
eic_states = strmatch('EIC', fn)';
elc_states = strmatch('ELC', fn)';
ifw_states = strmatch('IFW', fn)';
iiw_states = strmatch('IIW', fn)';
ilw_states = strmatch('ILW', fn)';
ifc_states = strmatch('IFC', fn)';
iic_states = strmatch('IIC', fn)';
ilc_states = strmatch('ILC', fn)';

% feature indices
FEATS = get_feature_set_mTIM();
IFW_FEAT = strmatch('intron_start_W', FEATS, 'exact');
ILW_FEAT = strmatch('intron_end_W', FEATS, 'exact');
IFC_FEAT = strmatch('intron_start_C', FEATS, 'exact');
ILC_FEAT = strmatch('intron_end_C', FEATS, 'exact');

MIN_SS_SCORE = 10;



state_seq = nan(size(label_seq));
state_seq(label_seq==LABELS.intergenic) = STATES.IGE;

EC = strmatch('exon_cover', FEATS, 'exact');
exon_cover = signal(EC,:);

% determine each gene's expression level in order to assign states of the
% corresponding discrete expression level
gene_blocks = find_blocks(label_seq ~= LABELS.intergenic);
N = size(gene_blocks,2);
for g=1:N,
  g_idx = gene_blocks(1,g):gene_blocks(2,g);
  exw_idx = g_idx(1) - 1 + find(label_seq(g_idx)==LABELS.exon_W);
  exc_idx = g_idx(1) - 1 + find(label_seq(g_idx)==LABELS.exon_C);

  if ~isempty(exw_idx),
    median_exo_cover = median(exon_cover(exw_idx));
    expr_level_W = find(PAR.expression_bins(:,1) <= median_exo_cover ...
                        & median_exo_cover < PAR.expression_bins(:,2));

    eiw_state = eiw_states(expr_level_W);
    state_seq(exw_idx) = eiw_state;
    
    iiw_state = iiw_states(expr_level_W);
    inw_idx = g_idx(1) - 1 + find(label_seq(g_idx)==LABELS.intron_W);
    state_seq(inw_idx) = iiw_state;
  end
  if ~isempty(exc_idx),
    median_exo_cover = median(exon_cover(exc_idx));
    expr_level_C = find(PAR.expression_bins(:,1) <= median_exo_cover ...
                        & median_exo_cover < PAR.expression_bins(:,2));
    eic_state = eic_states(expr_level_C);
    state_seq(exc_idx) = eic_state;
    
    iic_state = iic_states(expr_level_C);
    inc_idx = g_idx(1) - 1 + find(label_seq(g_idx)==LABELS.intron_C);
    state_seq(inc_idx) = iic_state;
  end
end

% intensity thresholds used to resolve ambiguous states
% 130114: 0.5 -> 0.4
ige_exo_threshold = 0.35 * (mean(exon_cover(label_seq==LABELS.exon_W | label_seq==LABELS.exon_C)) ...
                           + mean(exon_cover(label_seq==LABELS.intergenic)));
ino_exo_threshold = 0.35 * (mean(exon_cover(label_seq==LABELS.exon_W | label_seq==LABELS.exon_C)) ...
                           + mean(exon_cover(label_seq==LABELS.intron_W | label_seq==LABELS.intron_C)));
if isnan(ino_exo_threshold),
  ino_exo_threshold = ige_exo_threshold;
end

% resolve ambiguities
amb_blocks = find_blocks(isnan(state_seq));
problems = zeros(0,2);
for i=1:size(amb_blocks,2),
  l = amb_blocks(1,i)-1;
  r = amb_blocks(2,i)+1;

  if label_seq(l) == LABELS.intergenic && label_seq(r) == LABELS.exon_W,
    % decide where to split the ambiguous block between intergenic region
    % and exon
    for j=r-1:-1:l+1,
      if exon_cover(j) > ige_exo_threshold,
        state_seq(j) = state_seq(r);
      else
        state_seq(l+1:j) = STATES.IGE;
        break
      end
    end
  elseif label_seq(l) == LABELS.intergenic && label_seq(r) == LABELS.exon_C,
    % decide where to split the ambiguous block between intergenic region
    % and exon
    for j=r-1:-1:l+1,
      if exon_cover(j) > ige_exo_threshold,
        state_seq(j) = state_seq(r);
      else
        state_seq(l+1:j) = STATES.IGE;
        break
      end
    end
  elseif label_seq(l) == LABELS.exon_W && label_seq(r) == LABELS.intergenic,
    % decide where to split the ambiguous block between intergenic region
    % and exon
    for j=l+1:1:r-1,
      if exon_cover(j) > ige_exo_threshold,
        state_seq(j) = state_seq(l);
      else
        state_seq(j:r-1) = STATES.IGE;
        break
      end
    end
  elseif label_seq(l) == LABELS.exon_C && label_seq(r) == LABELS.intergenic,
    % decide where to split the ambiguous block between intergenic region
    % and exon
    for j=l+1:1:r-1,
      if exon_cover(j) > ige_exo_threshold,
        state_seq(j) = state_seq(l);
      else
        state_seq(j:r-1) = STATES.IGE;
        break
      end
    end
  elseif label_seq(l) == LABELS.exon_W && label_seq(r) == LABELS.intron_W,
    % decide where to put the splice site in the ambiguous block between
    % exon and intron
    [ss_score ss_idx] = max(signal(IFW_FEAT,l+1:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l;
      assert(l<ss_idx & ss_idx<=r);
      state_seq(l+1:1:ss_idx-1) = state_seq(l);
      state_seq(ss_idx:1:r-1) = state_seq(r);
    else
      % AIN'T A PROBLEM
      % problems = [problems; l+1, r-1];
      
      for j=l+1:1:r-1,
        if exon_cover(j) > ino_exo_threshold,
          state_seq(j) = state_seq(l);
        else
          state_seq(j:1:r-1) = state_seq(r);
          break
        end
      end
    end
  elseif label_seq(l) == LABELS.exon_C && label_seq(r) == LABELS.intron_C,
    % decide where to put the splice site in the ambiguous block between
    % exon and intron
    [ss_score ss_idx] = max(signal(IFC_FEAT,l+1:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l;
      assert(l<ss_idx & ss_idx<=r);
      state_seq(l+1:1:ss_idx-1) = state_seq(l);
      state_seq(ss_idx:1:r-1) = state_seq(r);
    else
      % AIN'T A PROBLEM
      % problems = [problems; l+1, r-1];

      for j=l+1:1:r-1,
        if exon_cover(j) > ino_exo_threshold,
          state_seq(j) = state_seq(l);
        else
          state_seq(j:1:r-1) = state_seq(r);
          break
        end
      end
    end
  elseif label_seq(l) == LABELS.intron_W && label_seq(r) == LABELS.exon_W,
    % decide where to put the splice site in the ambiguous block between
    % intron and exon
    [ss_score ss_idx] = max(signal(ILW_FEAT,l:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l - 1;
      assert(l<=ss_idx & ss_idx<r);
      state_seq(l+1:1:ss_idx) = state_seq(l);
      state_seq(ss_idx+1:1:r-1) = state_seq(r);
    else
      % AIN'T A PROBLEM
      % problems = [problems; l+1, r-1];
      
      for j=r-1:-1:l+1,
        if exon_cover(j) > ino_exo_threshold,
          state_seq(j) = state_seq(r);
        else
          state_seq(l+1:1:j) = state_seq(l);
          break
        end
      end
    end
  elseif label_seq(l) == LABELS.intron_C && label_seq(r) == LABELS.exon_C,
    % decide where to put the splice site in the ambiguous block between
    % intron and exon
    [ss_score ss_idx] = max(signal(ILC_FEAT,l:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l - 1;
      assert(l<=ss_idx & ss_idx<r);
      state_seq(l+1:1:ss_idx) = state_seq(l);
      state_seq(ss_idx+1:1:r-1) = state_seq(r);
    else
      % AIN'T A PROBLEM
      % problems = [problems; l+1, r-1];
 
      for j=r-1:-1:l+1,
        if exon_cover(j) > ino_exo_threshold,
          state_seq(j) = state_seq(r);
        else
          state_seq(l+1:1:j) = state_seq(l);
          break
        end
      end

    end
  elseif label_seq(l) == LABELS.exon_W && label_seq(r) == LABELS.exon_W,
    % fill in the whole ambiguous block with exon label
    if all(exon_cover(l+1:1:r-1) > ino_exo_threshold),
      state_seq(l:1:r) = state_seq(l);
    else
      [is_score is_idx] = max(signal(IFW_FEAT,l:r));
      [ie_score ie_idx] = max(signal(ILW_FEAT,l:r));
      is_idx = is_idx + l - 1;
      ie_idx = ie_idx + l - 1;
      level = str2num(fn{state_seq(l)}(end-1:end));
      if level == str2num(fn{state_seq(r)}(end-1:end)) && is_idx < ie_idx ...
                && is_score >= MIN_SS_SCORE && ie_score >= MIN_SS_SCORE ...
                && all(exon_cover(is_idx:1:ie_idx) < ino_exo_threshold),
        state_seq(l:1:is_idx-1) = state_seq(l);
        state_seq(ie_idx+1:1:r) = state_seq(r);        
        state_seq(is_idx:ie_idx) = iiw_states(level);
      else
        state_seq(l:1:r) = state_seq(l);
        problems = [problems; l+1, r-1];

4

      end
    end
  elseif label_seq(l) == LABELS.exon_C && label_seq(r) == LABELS.exon_C,
    % fill in the whole ambiguous block with exon label
    if all(exon_cover(l+1:1:r-1) > ino_exo_threshold),
      state_seq(l:1:r) = state_seq(l);
    else
      [is_score is_idx] = max(signal(IFC_FEAT,l:r));
      [ie_score ie_idx] = max(signal(ILC_FEAT,l:r));
      is_idx = is_idx + l - 1;
      ie_idx = ie_idx + l - 1;
      level = str2num(fn{state_seq(l)}(end-1:end));
      if level == str2num(fn{state_seq(r)}(end-1:end)) && is_idx < ie_idx ...
                && is_score >= MIN_SS_SCORE && ie_score >= MIN_SS_SCORE ...
                && all(exon_cover(is_idx:1:ie_idx) < ino_exo_threshold),
        state_seq(l:1:is_idx-1) = state_seq(l);
        state_seq(ie_idx+1:1:r) = state_seq(r);        
        state_seq(is_idx:ie_idx) = iic_states(level);
      else
        state_seq(l:1:r) = state_seq(l);
        problems = [problems; l+1, r-1];
 
 5
 
      end
    end
  elseif label_seq(l) == LABELS.intron_W && label_seq(r) == LABELS.intron_W,
    if all(exon_cover(l+1:1:r-1) <= ino_exo_threshold),
      % fill in the whole ambiguous block with intron label
      state_seq(l:1:r) = state_seq(l);
    else
      % try to heuristically recognize and label one exon within the ambiguous block
      [es_score es_idx] = max(signal(ILW_FEAT,l:r));
      [ee_score ee_idx] = max(signal(IFW_FEAT,l:r));
      es_idx = es_idx + l;
      ee_idx = ee_idx + l - 2;
      level = str2num(fn{state_seq(l)}(end-1:end));
      if level == str2num(fn{state_seq(r)}(end-1:end)) && es_idx < ee_idx ...
                && es_score >= MIN_SS_SCORE && ee_score >= MIN_SS_SCORE ...
                && all(exon_cover(es_idx:1:ee_idx) > ino_exo_threshold),
        state_seq(l:1:es_idx-1) = state_seq(l);
        state_seq(ee_idx+1:1:r) = state_seq(r);        
        state_seq(es_idx:ee_idx) = eiw_states(level);
      else
        state_seq(l:1:r) = state_seq(l);
        problems = [problems; l+1, r-1];

6

      end
    end
  elseif label_seq(l) == LABELS.intron_C && label_seq(r) == LABELS.intron_C,
    if all(exon_cover(l+1:1:r-1) <= ino_exo_threshold),
      % fill in the whole ambiguous block with intron label
      state_seq(l:1:r) = state_seq(l);
    else
      % try to heuristically recognize and label one exon within the ambiguous block
      [es_score es_idx] = max(signal(ILC_FEAT,l:r));
      [ee_score ee_idx] = max(signal(IFC_FEAT,l:r));
      es_idx = es_idx + l;
      ee_idx = ee_idx + l - 2;
      level = str2num(fn{state_seq(l)}(end-1:end));
      if level == str2num(fn{state_seq(r)}(end-1:end)) && es_idx < ee_idx ...
                && es_score >= MIN_SS_SCORE && ee_score >= MIN_SS_SCORE ...
                && all(exon_cover(es_idx:1:ee_idx) > ino_exo_threshold),
        state_seq(l:1:es_idx-1) = state_seq(l);
        state_seq(ee_idx+1:1:r) = state_seq(r);        
        state_seq(es_idx:ee_idx) = eic_states(level);
      else
        state_seq(l:1:r) = state_seq(l);
        problems = [problems; l+1, r-1];

7

      end
    end
  elseif label_seq(l) == LABELS.intergenic && label_seq(r) == LABELS.intron_W,
    % introduce an initial exon in ambiguous block between intergenic
    % region and intron
    % assert that there is at least one genic position that can be converted to exonic
    assert(l+1 <= r-1);
    exw_state = state_seq(r + find(label_seq(r+1:end)==LABELS.exon_W, 1, 'first'));
    [ss_score ss_idx] = max(signal(IFW_FEAT,l+1:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l;
      % assert that there is at least one genic position that can be converted to exonic
      assert(l+1<ss_idx & ss_idx<=r);
      state_seq(l+1:1:ss_idx-1) = exw_state;
      state_seq(ss_idx:1:r-1) = state_seq(r);
    else
      state_seq(l+1:1:r-1) = exw_state;      
      problems = [problems; l+1, r-1];

8

    end
  elseif label_seq(l) == LABELS.intergenic && label_seq(r) == LABELS.intron_C,
    % introduce an initial exon in ambiguous block between intergenic
    % region and intron
    % assert that there is at least one genic position that can be converted to exonic
    assert(l+1 <= r-1);
    exc_state = state_seq(r + find(label_seq(r+1:end)==LABELS.exon_C, 1, 'first'));
    [ss_score ss_idx] = max(signal(IFC_FEAT,l+1:r));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l;
      % assert that there is at least one genic position that can be converted to exonic
      assert(l+1<ss_idx & ss_idx<=r);
      state_seq(l+1:1:ss_idx-1) = exc_state;
      state_seq(ss_idx:1:r-1) = state_seq(r);
    else
      state_seq(l+1:1:r-1) = exc_state;      
      problems = [problems; l+1, r-1];
 
 9
 
    end
  elseif label_seq(l) == LABELS.intron_W && label_seq(r) == LABELS.intergenic,
    % introduce a terminal exon in ambiguous block between intron and intergenic region
    % assert that there is at least one genic position that can be converted to exonic 
    assert(l+1 <= r-1); 
    exw_state = state_seq(find(label_seq(1:l-1)==LABELS.exon_W, 1, 'last'));
    [ss_score ss_idx] = max(signal(ILW_FEAT,l:r-1));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l - 1;
      % assert that there is at least one genic position that can be converted to exonic
      assert(l<=ss_idx & ss_idx<r-1);
      state_seq(l+1:1:ss_idx) = state_seq(l);
      state_seq(ss_idx+1:1:r-1) = exw_state;            
    else
      state_seq(l+1:1:r-2) = exw_state;      

10      
      
      problems = [problems; l+1, r-1];
    end
  elseif label_seq(l) == LABELS.intron_C && label_seq(r) == LABELS.intergenic,
    % introduce a terminal exon in ambiguous block between intron and intergenic region
    % assert that there is at least one genic position that can be converted to exonic
    assert(l+1 <= r-1); 
    exc_state = state_seq(find(label_seq(1:l-1)==LABELS.exon_C, 1, 'last'));
    [ss_score ss_idx] = max(signal(ILC_FEAT,l:r-1));
    if ss_score >= MIN_SS_SCORE,
      ss_idx = ss_idx + l - 1;
      % assert that there is at least one genic position that can be converted to exonic
      assert(l<=ss_idx & ss_idx<r-1);
      state_seq(l+1:1:ss_idx) = state_seq(l);
      state_seq(ss_idx+1:1:r-1) = exc_state;
    else
      state_seq(l+1:1:r-2) = exc_state;

11      
      
      problems = [problems; l+1, r-1];
    end
  elseif label_seq(l) == LABELS.exon_W && label_seq(r) == LABELS.exon_C,
    % overlapping genes in antisense orientation, introduce a small
    % intergenic block
    % assert that there is are at least 3 positions that can be converted to intergenic
    if r-l <= 3,
      state_seq(l:r) = STATES.IGE;
    end 
    [mn min_idx] = min(exon_cover(l+2:r-2));
    if mn > ige_exo_threshold,

12        
        
        problems = [problems; l+1, r-1];
    end
    min_idx = min_idx + l+1;
    state_seq(min_idx-1:min_idx+1) = STATES.IGE;
    for j=min_idx-2:-1:l+1,
      if exon_cover(j) <= ige_exo_threshold,
        state_seq(j) = STATES.IGE;
      else
        state_seq(l+1:1:j) = state_seq(l);
        break
      end
    end
    for j=min_idx+2:1:r-1,
      if exon_cover(j) <= ige_exo_threshold,
        state_seq(j) = STATES.IGE;
      else
        state_seq(j:1:r-1) = state_seq(r);
        break
      end
    end
  elseif label_seq(l) == LABELS.exon_C && label_seq(r) == LABELS.exon_W,
    % overlapping genes in antisense orientation, introduce a small
    % intergenic block
    % assert that there are at least 3 positions that can be converted to intergenic
    if r-l <= 3,
      state_seq(l:r) = STATES.IGE;
    end 
    [mn min_idx] = min(exon_cover(l+2:r-2));
    if mn > ige_exo_threshold,

13        
        
        problems = [problems; l+1, r-1];
    end
    min_idx = min_idx + l+1;
    state_seq(min_idx-1:min_idx+1) = STATES.IGE;
    for j=min_idx-2:-1:l+1,
      if exon_cover(j) <= ige_exo_threshold,
        state_seq(j) = STATES.IGE;
      else
        state_seq(l+1:1:j) = state_seq(l);
        break
      end
    end
    for j=min_idx+2:1:r-1,
      if exon_cover(j) <= ige_exo_threshold,
        state_seq(j) = STATES.IGE;
      else
        state_seq(j:1:r-1) = state_seq(r);
        break
      end
    end
  else
%    warning('could not resolve ambiguities in label sequence');

14

    problems = [problems; l+1, r-1];
  end
end

% insert first and last intron states into the outermost positions of
% intron blocks
for l=1:PAR.num_levels,
  ifw_state = ifw_states(l);
  iiw_state = iiw_states(l);
  ilw_state = ilw_states(l);
  ino_blocks = find_blocks(state_seq==iiw_state);
  state_seq(ino_blocks(1,:)) = ifw_state;
  state_seq(ino_blocks(2,:)) = ilw_state;
  % check whether intron blocks are properly surrounded by corresponding
  % exon states
  if isempty(problems),
    eiw_state = eiw_states(l);
    assert(all(state_seq(ino_blocks(1,:)-1) == eiw_state));
    assert(all(state_seq(ino_blocks(2,:)+1) == eiw_state));
  end

  ifc_state = ifc_states(l);
  iic_state = iic_states(l);
  ilc_state = ilc_states(l);
  ino_blocks = find_blocks(state_seq==iic_state);
  state_seq(ino_blocks(1,:)) = ifc_state;
  state_seq(ino_blocks(2,:)) = ilc_state;
  % check whether intron blocks are properly surrounded by corresponding
  % exon states
  if isempty(problems),
    eic_state = eic_states(l);
%    assert(all(state_seq(ino_blocks(1,:)-1) == eic_state));
%    assert(all(state_seq(ino_blocks(2,:)+1) == eic_state));
  end
end
if isempty(problems),
  assert(~any(isnan(state_seq)));
end


% insert first and last exon states into the outermost positions of exon
% blocks 
g_start = [];
g_stop = [];
for l=1:PAR.num_levels,
  efw_state = efw_states(l);
  eiw_state = eiw_states(l);
  elw_state = elw_states(l);
  exo_blocks = find_blocks(state_seq==eiw_state);
  for e=1:size(exo_blocks,2),
    if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,
 
15        
        problems = [problems; exo_blocks(:,e)'];
    end

    if state_seq(exo_blocks(1,e)-1) == STATES.IGE,
      if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,
        state_seq(exo_blocks(1,e)-1) = efw_state;
        g_start = [g_start, exo_blocks(1,e)-1];
      else
        state_seq(exo_blocks(1,e)) = efw_state;
        g_start = [g_start, exo_blocks(1,e)];
      end
    end
    if state_seq(exo_blocks(2,e)+1) == STATES.IGE,
      if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,
        state_seq(exo_blocks(2,e)+1) = elw_state;
        g_stop = [g_stop, exo_blocks(2,e)+1];
      else
        state_seq(exo_blocks(2,e)) = elw_state;
        g_stop = [g_stop, exo_blocks(2,e)];
      end
    end

  end

  efc_state = efc_states(l);
  eic_state = eic_states(l);
  elc_state = elc_states(l);
  exo_blocks = find_blocks(state_seq==eic_state);
  for e=1:size(exo_blocks,2),
    if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,

16        
        
        problems = [problems; exo_blocks(:,e)'];
    end

    if state_seq(exo_blocks(1,e)-1) == STATES.IGE,
      if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,
        state_seq(exo_blocks(1,e)-1) = efc_state;
        g_start = [g_start, exo_blocks(1,e)-1];
      else
        state_seq(exo_blocks(1,e)) = efc_state;
        g_start = [g_start, exo_blocks(1,e)];
      end
    end
    if state_seq(exo_blocks(2,e)+1) == STATES.IGE,
      if exo_blocks(2,e)-exo_blocks(1,e)+1 < 3,
        state_seq(exo_blocks(2,e)+1) = elc_state;
        g_stop = [g_stop, exo_blocks(2,e)+1];
      else
        state_seq(exo_blocks(2,e)) = elc_state;
        g_stop = [g_stop, exo_blocks(2,e)];
      end
    end

  end
end
g_start = sort(g_start);
g_stop  = sort(g_stop);

%% CHECKS
if isempty(problems),
  ef_idx = find(ismember(state_seq, ef_states));
  assert(all(state_seq(ef_idx-1) == STATES.IGE));
  assert(all(ismember(state_seq(ef_idx+1),ei_states)));
  el_idx = find(ismember(state_seq, el_states));
  assert(all(ismember(state_seq(el_idx-1),ei_states)));
  assert(all(state_seq(el_idx+1) == STATES.IGE));

  ige_blocks = find_blocks(state_seq==STATES.IGE);
  assert(size(ige_blocks,2) >= 1);
  for b=2:size(ige_blocks,2),
%    assert(any(el_states == state_seq(ige_blocks(1,b)-1)));
  end
  for b=1:size(ige_blocks,2)-1,
 %   assert(any(ef_states == state_seq(ige_blocks(2,b)+1)));
  end
end


% adjust outermost exon boundaries to maximise agreement with coverage
% gradient feature
F_cover_grad = strmatch('cover_grad', get_feature_set_mTIM());
w = 100;
w = 20;
min_cover_grad = 0.8;
exon_states = sort([strmatch('EFW', fn)', strmatch('EIW', fn)', strmatch('ELW', fn)', ...
                    strmatch('EFC', fn)', strmatch('EIC', fn)', strmatch('ELC', fn)']);
for t=1:length(g_start),
  b = max(1, g_start(t) - w);
  % check whether b should be restricted because of close neighboring genes
  br = find(state_seq(1:g_start(t)-1) ~= STATES.IGE, 1, 'last') + 1;
  if ~isempty(br),
    assert(state_seq(br) == STATES.IGE)
    if  state_seq(br+2) == STATES.IGE,
      br = br + 2;
    end
    if br > b,
      b = br;
    end
  end

  e = min(length(state_seq), g_start(t) + w);
  % check whether e should be restricted to not extend beyond the first exon
  er = find(~ismember(state_seq(g_start(t)+1:end), ei_states), 1, 'first') ...
       - 1 + g_start(t);
  if ~isempty(er),
    assert(ismember(state_seq(er), ei_states))
    if ismember(state_seq(er-2), ei_states),
      er = er - 2;
    end
    if er < e,
      e = er;
    end
  end
  assert(state_seq(b) == STATES.IGE);
  assert(ismember(state_seq(e), ei_states));
  assert(b+2 <= e);

  % only relocate transcript ends if there is enough read evidence
  [mx mx_idx] = max(signal(F_cover_grad,b+1:e-1));  
  if mx > min_cover_grad,
    mx_idx = mx_idx + b;
    assert(abs(g_start(t)-mx_idx) <= w);
    state_seq(mx_idx) = state_seq(g_start(t));
    state_seq(b:mx_idx-1) = state_seq(b);
    state_seq(mx_idx+1:e) = state_seq(e);
  end
end

for t=1:length(g_stop),
  b = max(1, g_stop(t) - w);
  % check whether b should be restricted to not extend beyond the last exon
  br = find(~ismember(state_seq(1:g_stop(t)-1), ei_states), 1, 'last') + 1;
  if ~isempty(br),
    assert(ismember(state_seq(br), ei_states))
    if ismember(state_seq(br+2), ei_states),
      br = br + 2;
    end
    if br > b,
      b = br;
    end
  end

  e = min(length(state_seq), g_stop(t) + w);
  % check whether e should be restricted because of close neighboring genes
  er = find(state_seq(g_stop(t)+1:end) ~= STATES.IGE, 1, 'first') ...
       - 1 + g_stop(t);
  if ~isempty(er),
    assert(state_seq(er) == STATES.IGE)
    if state_seq(er-2) == STATES.IGE,
      er = er - 2;
    end
    if er < e,
      e = er;
    end
  end
  assert(ismember(state_seq(b), ei_states));
  assert(state_seq(e) == STATES.IGE);
  assert(b+2 <= e);

  % only relocate transcript ends if there is enough read evidence
  [mn mn_idx] = min(signal(F_cover_grad,b+1:e-1));
  if mn < -min_cover_grad,
    mn_idx = mn_idx + b;
    assert(abs(g_stop(t)-mn_idx) <= w);
    state_seq(mn_idx) = state_seq(g_stop(t));
    state_seq(b:mn_idx-1) = state_seq(b);
    state_seq(mn_idx+1:e) = state_seq(e);
  end
end


% assert that intergenic blocks are properly surrounded by EL and EF
% states
if isempty(problems),
  ef_idx = find(ismember(state_seq, ef_states));
  assert(all(state_seq(ef_idx-1) == STATES.IGE));
  assert(all(ismember(state_seq(ef_idx+1),ei_states)));
  el_idx = find(ismember(state_seq, el_states));
  assert(all(ismember(state_seq(el_idx-1),ei_states)));
  assert(all(state_seq(el_idx+1) == STATES.IGE));

  ige_blocks = find_blocks(state_seq==STATES.IGE);
  assert(size(ige_blocks,2) >= 1);
  for b=2:size(ige_blocks,2),
    if ~any(el_states == state_seq(ige_blocks(1,b)-1)),
      disp(state_seq(ige_blocks(1,b)-5:ige_blocks(1,b)+5));
    end
  end
  for b=1:size(ige_blocks,2)-1,
    if ~any(ef_states == state_seq(ige_blocks(2,b)+1)),
%       disp(state_seq(ige_blocks(2,b)-5:ige_blocks(2,b)+5));
    end
  end
end

if 0%~isempty(problems),
  clf; hold on; 
  plot(signal(1,:), 'b'); 
  plot(10*signal(3,:), 'r');
  plot(state_seq, 'k');
  plot(label_seq, ':k');
  for p=1:size(problems,1),
    x = [problems(p,1); problems(p,1); problems(p,2); problems(p,2)];
    y = [-90; -70; -70; -90];
    fill(x, y, 'r');
  end

  keyboard

end

% eof
