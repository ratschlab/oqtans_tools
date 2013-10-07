function ev = intron_eval(anno_genes, pred_genes)
% ev = intron_eval(anno_genes, pred_genes)
warning('Works only for chromosome size less than 100 mbp!');

anno_introns = get_unique_introns(anno_genes);
fprintf(' extracted from annotated transcripts.\n', size(anno_introns,1));
pred_introns = get_unique_introns(pred_genes);
fprintf(' extracted from predicted transcripts.\n', size(pred_introns,1));

anno_intron_matches = zeros(size(anno_introns,1),1);
pred_intron_matches = zeros(size(pred_introns,1),1);
strands = [+1, -1];
for s=1:length(strands),
  % TODO hack: this only works for a chromosome size <= 100 Mbp
  s_idx_anno = find(anno_introns(:,2) == strands(s));
  anno_intron_ids = anno_introns(s_idx_anno,3:4) + repmat(anno_introns(s_idx_anno,1).*10^8,1,2);
  s_idx_pred = find(pred_introns(:,2) == strands(s));
  pred_intron_ids = pred_introns(s_idx_pred,3:4) + repmat(pred_introns(s_idx_pred,1).*10^8,1,2);

  assert(issorted(anno_intron_ids, 'rows'));
  assert(issorted(pred_intron_ids, 'rows'));
  tol = 0;
  [idx_anno idx_pred] = compare_intervals_sorted(anno_intron_ids, pred_intron_ids, 'similar', tol);
  anno_intron_matches(s_idx_anno) = ~cellfun(@isempty, idx_anno);
  pred_intron_matches(s_idx_pred) = ~cellfun(@isempty, idx_pred);
end
assert(length(anno_intron_matches) == size(anno_introns,1))
assert(length(pred_intron_matches) == size(pred_introns,1))

fprintf('Compared %i annotated and %i predicted introns\n', size(anno_introns,1), ...
        size(pred_introns,1));
ev.sens = mean(anno_intron_matches);
ev.prec = mean(pred_intron_matches);
ev.f1 = 2 * ev.sens*ev.prec ./ (ev.sens+ev.prec);


