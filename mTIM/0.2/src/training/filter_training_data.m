function train_exm_id = filter_training_data(label, signal, exm_id_intervals, CFG)
%
% train_exm_id = filter_training_data(label, signal, exm_id_intervals, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

MIN_EXON_PROP = 0;%0.01;
MIN_IGE_BLOCKS = 1;%2
MAX_IGE_BLOCKS = 10; %10;
MAX_AMB_BLOCK = 1000;
MIN_EXPR_LEVEL = 0;% -inf; %1;
MAX_IGE_COVER = [100, 25];

fprintf('Filtering training examples...\n');
LABELS = get_label_set_mTIM();
FEATS = get_feature_set_mTIM();
len = exm_id_intervals(:,3) - exm_id_intervals(:,2) + 1;
is_good = (len <= CFG.max_train_chunk_len);
filtered_max_len = sum(~is_good);


filtered_exon_prop = 0;
filtered_ige_blocks = 0;
filtered_ige_cover = 0;
filtered_amb_label = 0;
filtered_adj_antisense = 0;
filtered_label_noise = 0;
filtered_label_problems = 0;
tic
%fh = figure
N = size(exm_id_intervals,1);
for i=1:N,
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);
  exo_idx = exm_id_intervals(i,2) - 1 + find(label(idx)==LABELS.exon_W | label(idx)==LABELS.exon_C);
  ino_idx = exm_id_intervals(i,2) - 1 + find(label(idx)==LABELS.intron_W | label(idx)==LABELS.intron_C);

  true_label_seq = label(idx);
  obs_seq = signal(:,idx);
  state_model = make_model_mTIM(CFG.PAR);
  [true_state_seq, problems] = labels_to_states(true_label_seq, state_model, obs_seq, CFG.PAR);
  if ~isempty(problems),
    is_good(i) = 0;
    filtered_label_problems = filtered_label_problems + 1;
  end

  if length(exo_idx) / length(idx) < MIN_EXON_PROP,
    is_good(i) = 0;
%    warning('proportion of annotated exons is %.2f%%', 100*length(exo_idx) / length(idx));
%    keyboard
    filtered_exon_prop = filtered_exon_prop + 1;
  end
  
  ige_blocks = find_blocks(label(idx) == LABELS.intergenic);
  if size(ige_blocks,2) < MIN_IGE_BLOCKS ...
        || size(ige_blocks,2) > MAX_IGE_BLOCKS,
    is_good(i) = 0;
%    warning('number of intergenic blocks is %i', size(ige_blocks,2));
%    keyboard
    filtered_ige_blocks = filtered_ige_blocks + 1;
  end
  
  if (~isfield(CFG.PAR,'train_filter_ige_cover') || ...
        (isfield(CFG.PAR,'train_filter_ige_cover') && ...
        ~CFG.PAR.train_filter_ige_cover)),

      FTC = strmatch('total_cover', FEATS, 'exact');
      total_cover = signal(FTC,idx);
      for b=1:size(ige_blocks,2),
        ige_cover_blocks = ige_blocks(1,b) - 1 + ...
            find_blocks(total_cover(ige_blocks(1,b):ige_blocks(2,b)) > 0);
        for c=1:size(ige_cover_blocks,2),
          if ige_cover_blocks(2,c)-ige_cover_blocks(1,c)+1 >= MAX_IGE_COVER(1) ...
              && mean(total_cover(ige_cover_blocks(1,c):ige_cover_blocks(2,c))) >= MAX_IGE_COVER(2),
            is_good(i) = 0;
            filtered_ige_cover = filtered_ige_cover + 1;
            break
          end
          if is_good(i) == 0,
            break
          end
        end
      end
  end
  
  % filter out ambiguous blocks with exon on one strand and intron on the
  % other strand flanking the ambiguous region
  if (~isfield(CFG.PAR,'train_filter_amb_label') || ...
        (isfield(CFG.PAR,'train_filter_amb_label') && ...
        ~CFG.PAR.train_filter_amb_label)),

      amb_blocks = find_blocks(label(idx) == LABELS.ambiguous);
      for j=1:size(amb_blocks,2),
        l = amb_blocks(1,j);
        r = amb_blocks(2,j);
        if (r - l + 1 > MAX_AMB_BLOCK) ...
              || (label(l)==LABELS.exon_W && label(r)==LABELS.intron_C) ...
              || (label(l)==LABELS.exon_C && label(r)==LABELS.intron_W) ...
              || (label(l)==LABELS.intron_W && label(r)==LABELS.exon_C) ...
              || (label(l)==LABELS.intron_C && label(r)==LABELS.exon_W) ...
              || (label(l)==LABELS.intron_W && label(r)==LABELS.intron_C) ...
              || (label(l)==LABELS.intron_C && label(r)==LABELS.intron_W),

          is_good(i) = 0;
          filtered_amb_label = filtered_amb_label + 1;
          break
        end
      end
  end

  % filter out directly adjacent labelings of gene models on different
  % strands
  W_pos = find(label(idx)==LABELS.exon_W | label(idx)==LABELS.intron_W);
  C_pos = find(label(idx)==LABELS.exon_C | label(idx)==LABELS.intron_C);
  if ~isempty(intersect(W_pos, C_pos+1)) ...
        || ~isempty(intersect(W_pos+1, C_pos)),
    is_good(i) = 0;
    filtered_adj_antisense = filtered_adj_antisense + 1;
  end

  % filter out examples with large discrepancies between label and read coverage 
  % also filter out genes which are too lowly expressed
  if (~isfield(CFG.PAR,'train_filter_label_noise') || ...
        (isfield(CFG.PAR,'train_filter_label_noise') && ...
        ~CFG.PAR.train_filter_label_noise)),
      
      FEC = strmatch('exon_cover', FEATS, 'exact');
      FIC = strmatch('intron_span', FEATS, 'exact');
      if mean(signal(FEC,exo_idx))    <= 2*mean(signal(FEC,setdiff(idx,exo_idx))) || ...
            mean(signal(FIC,ino_idx)) <= 2*mean(signal(FIC,setdiff(idx,ino_idx))) || ...
            mean(signal(FEC,exo_idx)) < MIN_EXPR_LEVEL,
        is_good(i) = 0;
        filtered_label_noise = filtered_label_noise + 1;
      end
  end 

%  is_good(i)
%  view_label_seqs(fh, signal(:,idx), label(idx));
%  keyboard

  if mod(i,100) == 0,
    fprintf('  %2.1f%% (%.1f sec)\r', 100*i/N, toc);
  end
end
fprintf('  %2.1f%% (%.1f sec)\n', 100*i/N, toc);

% TODO
filtered_max_len
filtered_exon_prop
filtered_ige_blocks
filtered_ige_cover
filtered_amb_label
filtered_adj_antisense
filtered_label_noise
filtered_label_problems

fprintf('  retained %i (%2.1f%%) potential training examples\n', ...
        sum(is_good), 100*sum(is_good)/length(is_good));

train_exm_id = exm_id_intervals(is_good,1);

fn = fieldnames(LABELS);
for f=1:length(fn),
  cnt = 0;
  total = 0;
  for i=1:length(train_exm_id),
    id = find(exm_id_intervals(:,1) == train_exm_id(i));
    assert(length(id) == 1);
    idx = exm_id_intervals(id,2):exm_id_intervals(id,3);
    cnt = cnt + sum(label(idx) == getfield(LABELS, fn{f}));
    total = total + length(idx);
  end
  fprintf('Labeled %12i positions in training sequences as %10s (%2.1f%%)\n', ...
          cnt, fn{f}, 100*cnt/total);
end




% eof
