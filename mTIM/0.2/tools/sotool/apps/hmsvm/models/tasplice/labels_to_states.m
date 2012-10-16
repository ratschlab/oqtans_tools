function [state_seq, problems] = labels_to_states(label_seq, state_model, signal, PAR)
% state_seq = labels_to_states(label_seq, state_model, signal, PAR)
% converts a label sequence into a state sequence

% written by Georg Zeller, MPI Tuebingen, Germany
problems = [];

LABELS = get_label_set();
STATES = get_state_set(PAR);

state_seq = nan(size(label_seq));

state_seq(label_seq==LABELS.intergenic) = STATES.ige;
state_seq(label_seq==LABELS.exonic)     = STATES.exo;
state_seq(label_seq==LABELS.intronic)   = STATES.ino;

state_seq(label_seq==LABELS.don_ss) = STATES.don;
state_seq(label_seq==LABELS.acc_ss) = STATES.acc;

% intensity thresholds used to resolve ambiguous states
ige_exo_threshold = 0.5 * (mean(signal(1,label_seq==LABELS.exonic)) ...
                           + mean(signal(1,label_seq== ...
                                         LABELS.intergenic)));
ino_exo_threshold = 0.5 * (mean(signal(1,label_seq==LABELS.exonic)) ...
                           + mean(signal(1,label_seq==LABELS.intronic)));

% resolve ambiguities for hybridization positions
amb_blocks = find_blocks(isnan(state_seq(1:2:end)));
amb_blocks = 2*(amb_blocks-1)+1;
for i=1:size(amb_blocks,2),
  l = amb_blocks(1,i)-2;
  r = amb_blocks(2,i)+2;

  if label_seq(l) == LABELS.intergenic & label_seq(r) == LABELS.exonic,
    for j=r-2:-2:l+2,
      if signal(1,j) > ige_exo_threshold,
        state_seq(j) = STATES.exo;
      else
        state_seq(l+2:2:j) = STATES.ige;
        break
      end
    end
  elseif label_seq(l) == LABELS.exonic & label_seq(r) == LABELS.intergenic,
    for j=l+2:2:r-2,
      if signal(1,j) > ige_exo_threshold,
        state_seq(j) = STATES.exo;
      else
        state_seq(j:2:r-2) = STATES.ige;
        break
      end
    end
  elseif label_seq(l) == LABELS.exonic & label_seq(r) == LABELS.intronic,
    don_idx = find(label_seq(l:r) == LABELS.don_ss) + l - 1;
    if isempty(don_idx),
      for j=l+2:2:r-2,
        if signal(1,j) > ino_exo_threshold,
          state_seq(j) = STATES.exo;
        else
          state_seq(j:2:r-2) = STATES.ino;
          break
        end
      end
    else
      assert(length(don_idx)==1);
      assert(l<don_idx & don_idx<r);
      state_seq(l+2:2:don_idx-1) = STATES.exo;
      state_seq(don_idx+1:2:r-2) = STATES.ino;
    end
  elseif label_seq(l) == LABELS.intronic & label_seq(r) == LABELS.exonic,
    acc_idx = find(label_seq(l:r) == LABELS.acc_ss) + l - 1;
    if isempty(acc_idx),
      for j=r-2:-2:l+2,
        if signal(1,j) > ino_exo_threshold,
          state_seq(j) = STATES.exo;
        else
          state_seq(l+2:2:j) = STATES.ino;
          break
        end
      end
    else
      assert(length(acc_idx)==1);
      assert(l<acc_idx & acc_idx<r);
      state_seq(l+2:2:acc_idx-1) = STATES.ino;
      state_seq(acc_idx+1:2:r-2) = STATES.exo;
    end
  elseif label_seq(l) == LABELS.exonic & label_seq(r) == LABELS.exonic,
    don_idx = find(label_seq(l:r) == LABELS.don_ss) + l - 1;
    acc_idx = find(label_seq(l:r) == LABELS.acc_ss) + l - 1;
    if ~isempty(don_idx) && ~isempty(acc_idx) && don_idx<acc_idx,
      assert(l<don_idx & don_idx<r);
      assert(l<acc_idx & acc_idx<r);
      state_seq(l+2:2:don_idx-1) = STATES.exo;
      state_seq(don_idx+1:2:acc_idx-1) = STATES.ino;
      state_seq(acc_idx+1:2:r-2) = STATES.exo;
    else
      state_seq(l+2:2:r) = STATES.exo;
    end
  elseif label_seq(l) == LABELS.intronic & label_seq(r) == LABELS.intronic,
    acc_idx = find(label_seq(l:r) == LABELS.acc_ss) + l - 1;
    don_idx = find(label_seq(l:r) == LABELS.don_ss) + l - 1;
    if ~isempty(acc_idx) && ~isempty(don_idx) && acc_idx<don_idx,
      assert(l<acc_idx & acc_idx<r);
      assert(l<don_idx & don_idx<r);
      state_seq(l+2:2:acc_idx-1) = STATES.ino;
      state_seq(acc_idx+1:2:don_idx-1) = STATES.exo;
      state_seq(don_idx+1:2:r-2) = STATES.ino;
    else
      state_seq(l+2:2:r) = STATES.ino;
    end

  else
    label_seq(l)
    label_seq(r)
    l
    r
    keyboard
  end
end

% resolve ambiguities for splice site positions
idx = find(label_seq==LABELS.no_ss | label_seq==LABELS.double_ss);
assert(all(mod(idx,2)==0));
for i=1:length(idx),
  l = idx(i) - 1;
  assert(~isnan(state_seq(l)));
  r = idx(i) + 1;
  assert(~isnan(state_seq(r)));

  if state_seq(l) == STATES.ige & state_seq(r) == STATES.ige,
    state_seq(idx(i)) = STATES.ige_ss;
  elseif state_seq(l) == STATES.exo & state_seq(r) == STATES.exo,
    state_seq(idx(i)) = STATES.exo_ss;
  elseif state_seq(l) == STATES.ino & state_seq(r) == STATES.ino,
    state_seq(idx(i)) = STATES.ino_ss;
  elseif state_seq(l) == STATES.ige & state_seq(r) == STATES.exo,
    state_seq(idx(i)) = STATES.trans_start;
  elseif state_seq(l) == STATES.exo & state_seq(r) == STATES.ige,
    state_seq(idx(i)) = STATES.trans_end;
  elseif state_seq(l) == STATES.exo & state_seq(r) == STATES.ino,
    state_seq(idx(i)) = STATES.don;
%    don_idx = find(label_seq(l:r) == LABELS.don_ss);
%    if idx(i)<don_idx+l-1,
%      state_seq(idx(i)) = STATES.exo_ss;
%    else
%      state_seq(idx(i)) = STATES.ino_ss;
%    end
  elseif state_seq(l) == STATES.ino & state_seq(r) == STATES.exo,
    state_seq(idx(i)) = STATES.acc;
%    acc_idx = find(label_seq(l:r) == LABELS.acc_ss);
%    if idx(i)<acc_idx+l-1,
%      state_seq(idx(i)) = STATES.ino_ss;
%    else
%      state_seq(idx(i)) = STATES.exo_ss;
%    end
  else
    label_seq(l)
    label_seq(r)
    l
    r
    keyboard
  end
end


assert(~any(isnan(state_seq)));
% eof