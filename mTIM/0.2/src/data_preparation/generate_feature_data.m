function generate_feature_data(CFG)

% generate_feature_data(CFG)
%
% Fills chunks with label and feature sequences.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011


for fold=1:CFG.num_xval_folds,
  fprintf('Filling data blocks for cross-validation fold %i\n', fold);

  d = CFG.xval_dirs{fold}

  % fill training chunks
  fn = [d 'train_data'];
  if ~exist(fn, 'file'),
    load([fn '.mat'],  'train_chunks', 'train_chunk_ids');
  else,
    load(fn,  'train_chunks', 'train_chunk_ids');
  end
  [signal exm_id_intervals label] = fill_chunks(train_chunks, CFG);
  save(fn, 'train_chunks', 'train_chunk_ids', ...
         'label', 'exm_id_intervals', 'CFG', 'fold');
  
  clear label exm_id_intervals train_chunks train_chunk_ids
  fn = [d 'train_data_signal'];
  save_struct(fn, signal, 'signal', 10e8);
  clear signal
  
  if CFG.fill_vald_chunks,
    % fill validation chunks
    fn = [d 'vald_data.mat'];
    load(fn, 'vald_chunks', 'vald_chunk_ids');
    [signal exm_id_intervals label] = fill_chunks(vald_chunks, CFG);
    save(fn, 'vald_chunks', 'vald_chunk_ids', ...
         'label', 'exm_id_intervals', 'CFG', 'fold');

    clear label exm_id_intervals vald_chunks vald_chunk_ids
    fn = [d 'vald_data_signal'];
    save_struct(fn, signal, 'signal', 10e8);
    clear signal
  end
  
  if CFG.fill_test_chunks,
    % fill test chunks
    fn = [d 'test_data.mat'];
    load(fn, 'test_chunks', 'test_chunk_ids');
    [signal exm_id_intervals label] = fill_chunks(test_chunks, CFG);
    save(fn, 'test_chunks', 'test_chunk_ids', ...
         'label', 'exm_id_intervals', 'CFG', 'fold');

    clear label exm_id_intervals test_chunks test_chunk_ids
    fn = [d 'test_data_signal'];
    save_struct(fn, signal, 'signal', 10e8);
    clear signal
  end
  fprintf('\n');
end

% eof
