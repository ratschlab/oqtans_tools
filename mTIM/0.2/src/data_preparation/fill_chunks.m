function [signal exm_id_intervals label] = fill_chunks(chunks, CFG)

% [signal exm_id_intervals label] = fill_chunks(chunks, CFG)
%
% Fills chunks with label and feature sequences.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

% set CHECK=1 to enable checking the predictivity of all features, this is
% especially useful after changing feature generation procedures
CHECK = 0;

assert(issorted(chunks, 'rows'));
total_len = sum(chunks(:,3)-chunks(:,2) + 1);

if nargout > 2,
  if ~exist(CFG.label_fn, 'file'),
     load([CFG.label_fn '.mat'], 'LABELS', 'label_map');
  else
     load(CFG.label_fn, 'LABELS', 'label_map');
  end
  [label exm_id_intervals] = label_chunks(chunks, label_map);
  assert(length(label) == total_len);
end

fprintf('Setup empty signal with length %i.\n', total_len);
signal = zeros(0, total_len);
% append features derived from read alignments to feature matrix
signal = append_read_feats(chunks, signal, CFG);
% append predicted splice sites to feature matrix
signal = append_splice_feats(chunks, signal, CFG);
% append pair features
if (CFG.PAR.use_pair_feats),
    signal = append_pair_feats(chunks, signal, CFG);
else
    signal = [signal; zeros(3,size(signal,2))];
    fprintf('fill_chunks: dont generate pair feature.\n');
end
% TODO: reanimante
%  signal = append_repeat_feats(chunks, exp, signal, CFG);
signal = [signal; zeros(1,size(signal,2))];
% add binned read features (here: only 3 binned intron_span features)
signal = append_binned_read_feats(chunks, signal, CFG);
% add the cufflinks feature chunk (TODO:VALIDATE)
signal = append_cufflinks_feats(chunks, signal, CFG);

num_signals = size(signal,1);
approx_mem = (num_signals*total_len*8)/(1024.^3);
fprintf('Finished signal-generation. Signal has now a size of (%i,%i) which is approx %3.2fGB.\n', ...
	num_signals, total_len, approx_mem);

% check predictivity of all features collected
if nargout > 2 && CHECK,
  sample_size = 100000;
  check_features(signal, label, sample_size);
end

% eof
