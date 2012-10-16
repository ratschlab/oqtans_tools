function genes = norm_expr_levels(genes, num_levels)

gene_levels = nan(1,length(genes));
for g=1:length(genes),
  gene_levels(g) = round(mean([genes(g).expr]));
end
levels = unique(round(gene_levels));

if length(levels) > num_levels,
  L = num_levels;
  bins = [];
  val = sort(gene_levels);
  n = L-1;
  while length(bins) < L-1;
    bins = round(linspace(1, length(val), n+2));
    % the following unique statement can result in a reduced number of
    % limits; that's why we increase n in the loop (but in most cases we
    % will just pass through the while loop once)
    bins = bins(2:end-1);
    bins = unique(val(bins));
    n = n + 1;

    if n > 5*L,
      keyboard
    end
  end
  bins = [[-inf; bins'], [bins'; inf]];
  assert(all(bins(1:end-1,1) < bins(2:end,1)));
  
  % assign a discrete expression level to each gene
  discrete_levels = nan(length(gene_levels),1);
  for i=1:size(bins,1),
    idx = find(bins(i,1) <= gene_levels & gene_levels < bins(i,2));
    discrete_levels(idx) = i;
  end
  assert(all(discrete_levels >= 1 & discrete_levels <= L));
  for g=1:length(genes),
    genes(g).epxr = discrete_levels(g);
  end
else
  warning('found only %i levels', length(levels));
end
