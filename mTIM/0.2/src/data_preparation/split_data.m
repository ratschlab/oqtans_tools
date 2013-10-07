function chunks = split_data(CFG, chunks)

% chunks = split_data(CFG, chunks)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

% seed for random number generation
rand('seed', 11081979);

%%% assign chunks to subsets

% chunk subset index will be the fifth column vector
chunks(:,5) = zeros(size(chunks,1),1);
% chunk ids are in the fourth column vector
cid = chunks(:,4);
assert(isequal(unique(cid), cid));
assert(~any(cid == 0));

% assign several adjacent chunks to the same subset in oder to avoid to
% many breakpoints between training and test chunks
cid = zeros(size(cid));
cid(1) = 1;
unit_length = ceil(sum(CFG.chr_lens) / (CFG.num_xval_folds * 100) / 1000) * 1000;
len = chunks(1,3) - chunks(1,2) + 1;
for c=2:size(chunks,1),
  if chunks(c,1) == chunks(c-1,1), % same chr
    len = len + chunks(c,3) - chunks(c,2) + 1;
    if len < unit_length,
      cid(c) = cid(c-1);
    else
      len = chunks(c,3) - chunks(c,2) + 1;
      cid(c) = cid(c-1) + 1;
    end
  else
    len = chunks(c,3) - chunks(c,2) + 1;
    cid(c) = cid(c-1) + 1;
  end
end
ucid = unique(cid);
ucid = ucid(randperm(length(ucid)));
subset_ends = round(linspace(0, length(ucid), CFG.num_data_subsets+1));
for i=1:CFG.num_data_subsets,
  subset_ids{i} = sort(ucid(subset_ends(i)+1:subset_ends(i+1)))';
  idx = find(ismember(cid, subset_ids{i}));
  chunks(idx,5) = i;
end
%cid = cid(randperm(length(cid)));
%subset_ends = round(linspace(0, length(cid), CFG.num_data_subsets+1));
%for i=1:CFG.num_data_subsets,
%  subset_ids{i} = sort(cid(subset_ends(i)+1:subset_ends(i+1)))';
%  idx = find(ismember(chunks(:,4), subset_ids{i}));
%  chunks(idx,5) = i;
%end
assert(isequal(sort([subset_ids{:}])', sort(ucid)));
assert(~any(chunks(:,5) == 0));
assert(isequal(unique(chunks(:,5)'), 1:CFG.num_data_subsets));

for i=1:CFG.num_data_subsets,
  idx = find(chunks(:,5) == i);
  fprintf('  subset %i contains %i chunks w/ a total length of %i\n', ...
          i, length(idx), sum(chunks(idx,3)-chunks(idx,2)+1));
end
clear subset_ids idx cid

%%% assign examples to training, validation and test sets for cross-validation

assert(CFG.num_xval_folds >= 3);  % otherwise the following data splitting
                                  % will not work!
perm  = 1:CFG.num_data_subsets;   % a permutation of subsets
train = 1:CFG.num_data_subsets-2; % use N-2 subsets for training
vald  = CFG.num_data_subsets-1;   % use 1 subset for validation
test  = CFG.num_data_subsets;     % use 1 subset for testing

for fold=1:CFG.num_xval_folds,
  train_subsets = perm(train);
  train_chunk_ids = find(ismember(chunks(:,5), train_subsets));
  train_chunks = chunks(train_chunk_ids,:);
  train_chunk_ids = chunks(train_chunk_ids,4);

  vald_subsets = perm(vald);
  vald_chunk_ids = find(ismember(chunks(:,5), vald_subsets));
  vald_chunks = chunks(vald_chunk_ids,:);
  vald_chunk_ids = chunks(vald_chunk_ids,4);
  
  test_subsets = perm(test);
  test_chunk_ids = find(ismember(chunks(:,5), test_subsets));
  test_chunks = chunks(test_chunk_ids,:);
  test_chunk_ids = chunks(test_chunk_ids,4);
  
  assert(isempty(intersect(train_chunk_ids, vald_chunk_ids)));
  assert(isempty(intersect(test_chunk_ids, [train_chunk_ids; vald_chunk_ids])));

  d = CFG.xval_dirs{fold};

  fn = [d 'train_data.mat'];
  save(fn, 'train_chunks', 'train_chunk_ids', 'CFG', 'fold');
  fn = [d 'vald_data.mat'];
  save(fn, 'vald_chunks', 'vald_chunk_ids', 'CFG', 'fold');
  fn = [d 'test_data.mat'];
  save(fn, 'test_chunks', 'test_chunk_ids', 'CFG', 'fold');
  
  tmp = perm(end);
  perm(2:end) = perm(1:end-1);
  perm(1) = tmp;
end

