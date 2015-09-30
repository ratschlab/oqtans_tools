function train_exm_id = min_filter_training_data(label, signal, exm_id_intervals, CFG)
%
% train_exm_id = filter_training_data(label, signal, exm_id_intervals, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009
% adapted by Nico Görnitz, TU Berlin, 2013


fprintf('Filtering training examples...\n');
LABELS = get_label_set_mTIM();
FEATS = get_feature_set_mTIM();
len = exm_id_intervals(:,3) - exm_id_intervals(:,2) + 1;
is_good = (len <= CFG.max_train_chunk_len);


filtered_label_problems = 0;
tic

N = size(exm_id_intervals,1);
for i=1:N,
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);
  exo_idx = exm_id_intervals(i,2) - 1 + find(label(idx)==LABELS.exon_W | label(idx)==LABELS.exon_C);
  ino_idx = exm_id_intervals(i,2) - 1 + find(label(idx)==LABELS.intron_W | label(idx)==LABELS.intron_C);

  true_label_seq = label(idx);
  obs_seq = signal(:,idx);
  state_model = make_model_mTIM(CFG.PAR);
  [true_state_seq, problems] = labels_to_states(true_label_seq, state_model, obs_seq, CFG.PAR);
  % any non-decodable sequence?
  if any(isnan(true_state_seq)),
    is_good(i) = 0;
    filtered_label_problems = filtered_label_problems + 1;
  end
  if mod(i,100) == 0,
    fprintf('  %2.1f%% (%.1f sec)\r', 100*i/N, toc);
  end
end
fprintf('  %2.1f%% (%.1f sec)\n', 100*i/N, toc);

% TODO
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
