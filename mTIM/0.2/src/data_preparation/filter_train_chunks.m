function filter_train_chunks(CFG)

% filter_train_chunks(CFG)
%
% Throw away some training chunks to avoid filling efforts
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011

sample_size = inf;%100%inf;

for fold=1:CFG.num_xval_folds,
  fprintf('Filtering data blocks for cross-validation fold %i\n', fold);

  d = CFG.xval_dirs{fold};
  
  fn = [d 'train_data_complete.mat'];
  if ~exist(fn, 'file'),
    fn = [d 'train_data.mat'];
  end
  fprintf('Loading training data from %s', fn);

  load(fn, 'train_chunks', 'train_chunk_ids');
  size(train_chunks)
  iss = issorted(train_chunks, 'rows');
  assert(iss);
  fn = [d 'train_data_complete.mat'];
  save(fn, 'train_chunks', 'train_chunk_ids', 'CFG', 'fold');
  
  L_orig = sum(train_chunks(:,3)-train_chunks(:,2)+1);
  L_disc = 0;
  L_keep = 0;
  keep_chunks = ones(size(train_chunks,1), 1);
  for c=1:size(train_chunks,1),
    len = train_chunks(c,3) - train_chunks(c,2) + 1;
    if len > CFG.max_train_chunk_len,
      keep_chunks(c) = 0;
      L_disc = L_disc + len;
    else
      L_keep = L_keep + len;
    end
  end
  fprintf('  discarded  %5i chunks with a total length of %9i (%2.1f%%)\n', ...
          sum(keep_chunks==0), L_disc, 100*L_disc/L_orig);
  fprintf('  retained   %5i chunks with a total length of %9i (%2.1f%%)\n', ...
          sum(keep_chunks==1), L_keep, 100*L_keep/L_orig);
  train_chunks = train_chunks(find(keep_chunks),:);
  train_chunk_ids = train_chunk_ids(find(keep_chunks));
  
  if length(train_chunk_ids) > sample_size,
    r = randperm(length(train_chunk_ids));
    r = r(1:sample_size);
    train_chunks = train_chunks(r,:);
    train_chunk_ids = train_chunk_ids(r);
  end
  L_sample = sum(train_chunks(:,3)-train_chunks(:,2)+1);
  fprintf('  subsampled %5i chunks with a total length of %9i (%2.1f%%)\n', ...
          sample_size, L_sample, 100*L_sample/L_orig);
  
  [train_chunks perm] = sortrows(train_chunks);
  train_chunk_ids = train_chunk_ids(perm);
  
  fn = [d 'train_data.mat'];
  save(fn, 'train_chunks', 'train_chunk_ids', 'CFG', 'fold');
  
  fprintf('\n');
end

% eof
