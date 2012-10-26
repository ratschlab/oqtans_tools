function fn_data = save_reduced_data(fn_data, exm_ids, CFG);

% fn_data = save_reduced_data(fn_data, exm_ids, CFG);
%
% Subselects training examples (specified by exm_ids) and re-saves the
% reduced data set dt fn_data (which is also returned).
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

tmp = load(fn_data);
fold = tmp.fold;
assert(all(ismember(exm_ids, tmp.train_chunk_ids)));
assert(all(ismember(exm_ids, tmp.train_chunks(:,4))));

if ~isfield(tmp, 'signal'),
  if isequal(fn_data(end-3:end), '.mat'),
    fn_data_signal = fn_data(1:end-4);
  else
    fn_data_signal = fn_data;
  end
  fn_data_signal = [fn_data_signal '_signal'];
  if(~exist(fn_data_signal)),
     tmp.signal = load_struct([fn_data_signal '.mat'], 'signal');
  else
     tmp.signal = load_struct(fn_data_signal, 'signal');
  end
  
end

fprintf('Subselecting %i training examples (out of %i, %2.1f%%)\n', ...
        length(exm_ids), length(tmp.train_chunk_ids), ...
        100*length(exm_ids)/length(tmp.train_chunk_ids));

% subselect the sequences actually needed for training
idx = find(ismember(tmp.train_chunk_ids, exm_ids));
exm_id_intervals = tmp.exm_id_intervals(idx,:);
train_chunks = tmp.train_chunks(idx,:);
train_chunk_ids = tmp.train_chunk_ids(idx);

N = sum(exm_id_intervals(:,3)-exm_id_intervals(:,2)+1);
label = zeros(1, N, 'int8');
signal = zeros(size(tmp.signal,1), N);
cnt = 0;
for i=1:size(exm_id_intervals,1),
  % idx points into the old sequence structures
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);
  l = tmp.label(idx);
  s = tmp.signal(:,idx);

  % update exm_id_intervals
  exm_id_intervals(i,2) = cnt + 1;
  exm_id_intervals(i,3) = cnt + length(idx);
  cnt = cnt + length(idx);
  
  % now idx points into the new structures
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);
  label(idx) = l;
  signal(:,idx) = s;
end

% save select data
if isequal(fn_data(end-3:end), '.mat'),
  fn_data = fn_data(1:end-4);
end
whosinfo = whos('signal');
if  whosinfo.bytes > 1e9, % about 1Gb
  save([fn_data '_select_' CFG.start_time '.mat'], 'CFG', 'fold', ...
       'exm_id_intervals', 'train_chunks', 'train_chunk_ids', 'label');
  save_struct([fn_data '_select_' CFG.start_time '_signal'], signal, ...
              'signal');
  fn_data = [fn_data '_select_' CFG.start_time '.mat'];
else
  fn_data = [fn_data '_select_' CFG.start_time '.mat'];
  save(fn_data, 'CFG', 'fold', 'exm_id_intervals', 'train_chunks', ...
       'train_chunk_ids', 'label', 'signal');
end
% eof
