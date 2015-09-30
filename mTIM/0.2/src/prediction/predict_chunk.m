function [pred_genes pred_label] = predict_chunk(ARGS)

% [pred_genes pred_label] = predict_chunk(ARGS)
%
%
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

assert(size(ARGS.chunk,1) == 1);

% in case of parallel prediction add paths
if (ARGS.CFG.grid_use),
    % include user-specified include paths
    if isfield(ARGS.CFG.PAR, 'include_paths'),
        for i=1:length(ARGS.CFG.PAR.include_paths),
            addpath(ARGS.CFG.PAR.include_paths{i});
        end
    end
end
% init paths for hmsvm
init_hmsvm_paths();

LABELS = get_label_set_mTIM();

% obtain features
if ~ARGS.CFG.grid_use && ARGS.CFG.verbose>=3,
  % additionally load true label for visualization (see below)
  [signal exm_id_interval true_label] = fill_chunks(ARGS.chunk, ARGS.CFG);
else
  signal = fill_chunks(ARGS.chunk, ARGS.CFG);
end


% if this is a zero-shot transfer learning prediction
% and there exist signal cumulants then transform the signal first
fprintf('Check for zero-shot transfer learning...\n');
if isfield(ARGS.CFG.PAR,'transfer'),

    % zero-shot transfer learning?
    if isfield(ARGS.CFG.PAR.transfer,'zero_shot') ...
        && ARGS.CFG.PAR.transfer.zero_shot ...
        && isfield(ARGS.pred_PAR,'transfer'),
        
        fprintf('Zero-shot transfer learning active, hence transforming the signal...\n');
        transfer = ARGS.pred_PAR.transfer;
        signal = set_signal_cumulants(signal, ...
            transfer.sig_means,transfer.sig_stds,transfer.sig_valid);
    end

end


% suppress features
signal = signal(ARGS.pred_PAR.suppressed_feats_map(:,2)>0,:);
%[ARGS.pred_PAR signal ARGS.score_plifs] = suppress_features(ARGS.pred_PAR, signal, ARGS.score_plifs, ARGS.state_model);

% Viterbi decoding
pred_path = decode_viterbi(signal, ARGS.transition_scores, ARGS.score_plifs, ARGS.pred_PAR);

offset = ARGS.chunk(2) - 1;
pred_genes = label_to_genes(pred_path.label_seq, offset);

for g=1:length(pred_genes),
  pred_genes(g).chr_num = ARGS.chunk(1);
  pred_genes(g).chr = ARGS.CFG.chr_names{pred_genes(g).chr_num};
  idx = [pred_genes(g).start:pred_genes(g).stop] - offset;
  exo_idx = find(pred_path.label_seq(idx) == LABELS.exon_W ...
                 | pred_path.label_seq(idx) == LABELS.exon_C);
  assert(~isempty(exo_idx));
  pred_genes(g).expr = mean(states_to_levels(pred_path.state_seq(idx(exo_idx)), ...
                                             ARGS.pred_PAR));
  assert(all(pred_genes(g).exons{1}(:,1) >= pred_genes(g).start));
  assert(all(pred_genes(g).exons{1}(:,2) <= pred_genes(g).stop));
  assert(1 <= pred_genes(g).expr && pred_genes(g).expr <= ARGS.pred_PAR.num_levels);
end

if nargout>1,
  pred_label = pred_path.label_seq;
end

if ~ARGS.CFG.grid_use && ARGS.CFG.verbose>=3,
  fh = figure();
  view_label_seqs(fh, signal, true_label, pred_path.label_seq);
  keyboard
end

% eof
