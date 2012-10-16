function [ev, split, fusion] = gene_eval(anno_genes, pred_genes)
% ev = gene_eval(anno_genes, pred_genes)

split.anno_inds = [];
split.pred_inds = [];

fusion.anno_inds = [];
fusion.pred_inds = [];

anno_gene_matches = zeros(length(anno_genes),1);
pred_gene_matches = zeros(length(pred_genes),1);
strands = {'+', '-'};
for s=1:length(strands),
  % First, identify coverlapping genes on the same strand
  % TODO hack: this only works for a chromosome size <= 100 Mbp
  s_idx_anno = find([anno_genes.strand] == strands{s});
  anno_gene_ids = [[anno_genes(s_idx_anno).start]', [anno_genes(s_idx_anno).stop]'] ...
      + repmat([anno_genes(s_idx_anno).chr_num]' * 10^8, 1, 2);
  s_idx_pred = find([pred_genes.strand] == strands{s});
  pred_gene_ids = [[pred_genes(s_idx_pred).start]', [pred_genes(s_idx_pred).stop]'] ...
      + repmat([pred_genes(s_idx_pred).chr_num]' * 10^8, 1, 2);

%  tic
  [anno_gene_ids anno_perm] = sortrows(anno_gene_ids);
  s_idx_anno = s_idx_anno(anno_perm);
  [pred_gene_ids pred_perm] = sortrows(pred_gene_ids);
  s_idx_pred = s_idx_pred(pred_perm);
%  assert(issorted(anno_gene_ids, 'rows'));
%  assert(issorted(pred_gene_ids, 'rows'));
  ovl = 1;
  [ovl_anno ovl_pred] = compare_intervals_sorted(anno_gene_ids, pred_gene_ids, 'overlap', ovl);

  % find predicted genes that were splitted into >1 pieces
  split_lens = cellfun(@length,ovl_anno);
  inds = find(split_lens>1);
  split.anno_inds = [split.anno_inds, s_idx_anno(inds)];
  for i=1:length(inds),
    foo = ovl_anno{inds(i)};
    split.pred_inds = [split.pred_inds, s_idx_pred(foo)];
  end

  % find predicted genes that were merged 
  fusion_lens = cellfun(@length,ovl_pred);
  inds = find(fusion_lens>1);
  fusion.pred_inds = [fusion.pred_inds, s_idx_pred(inds)];
  for i=1:length(inds),
    foo = ovl_pred{inds(i)};
    fusion.anno_inds = [fusion.anno_inds, s_idx_anno(foo)];
  end


%  fprintf('Interval comparison between annotated and predicted genes took %.1f sec\n', toc);
  % Second, compare the intron structure of all of their transcripts
  tic
  for i=1:length(ovl_anno),
    % see whether any pred_gene overlapping with the current anno_gene
    % has a matching transcript (one with identical introns) 
    anno_gene_matches(s_idx_anno(i)) = ...
        match_transcripts(anno_genes(s_idx_anno(i)), pred_genes(s_idx_pred(ovl_anno{i})));
%    fprintf('matched %i of %i annotated genes (%2.1f%%, %.1f sec)\r', i, ...
%            length(ovl_anno), 100*i/length(ovl_anno), toc);
  end
%  fprintf('Sensitivity assessment took %.1f sec                          \n', toc);
  tic
  for i=1:length(ovl_pred),
    % see whether any anno_gene overlapping with the current pred_gene
    % has a matching transcript (one with identical introns) 
    pred_gene_matches(s_idx_pred(i)) = ...
        match_transcripts(pred_genes(s_idx_pred(i)), anno_genes(s_idx_anno(ovl_pred{i})));
%    fprintf('matched %i of %i predicted genes (%2.1f%%, %.1f sec)\r', i, ...
%            length(ovl_pred), 100*i/length(ovl_pred), toc);
  end
%  fprintf('Precision assessment took %.1f sec                            \n', toc);
end
fprintf('Compared %i annotated and %i predicted genes\n', length(anno_genes), ...
        length(pred_genes));
ev.sens = mean(anno_gene_matches);
ev.prec = mean(pred_gene_matches);
ev.f1 = 2 * ev.sens*ev.prec ./ (ev.sens+ev.prec);


split.anno_inds = unique(split.anno_inds);
split.pred_inds = unique(split.pred_inds);

fusion.anno_inds = unique(fusion.anno_inds);
fusion.pred_inds = unique(fusion.pred_inds);


