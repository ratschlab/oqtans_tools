function bins = discretize_expression(exon_cover, label, exm_id_intervals, L)
%
% bins = discretize_expression(exon_cover, label, exm_id_intervals, L)
%
% Computes bin boundaries for L discrete levels of exonic coverage such
% that each bin contains approximately equally many examples.
%
% exon_cover -- a vector of counts indicating how many reads are mapped
%   to the genome at a given position (of length N)
% label -- a labeling of probes (of length N) used to determine exonic
%   probes
% exm_id_intervals -- a matrix with 3 columns of length N (pos. integers)
%   indicating to which example a label belongs to. First column indicates
%   the example id, the following two columns the range of the corresponding
%   sequences in the coverage and label vector
% L -- the desired number of discrete expression levels
% returns L bin boundaries (as a vector of size L x 2)
%
% see get_label_set_mTIM.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011

fprintf('Computing discrete expression levels across genes...\n');
LABELS = get_label_set_mTIM();

gene_blocks = find_blocks(label ~= LABELS.intergenic);
N = size(gene_blocks,2);
median_exo_cover = nan(N,1);
for g=1:N,
  g_idx = gene_blocks(1,g):gene_blocks(2,g);
  exo_idx = g_idx(1) - 1 + find(label(g_idx)==LABELS.exon_W | label(g_idx)==LABELS.exon_C);
  if length(exo_idx) > 0,
    median_exo_cover(g) = median(exon_cover(exo_idx));
  else
%    warning('no exons labeled')
    median_exo_cover(g) = NaN;
  end
end

% have one bin between -inf and 1 (min_expr_level) and adjust bin boundaries
% for the remainder such that each bin contains approximately the same
% number of genes
min_expr_level = 1;
%min_expr_level = -inf;
idx = find(~isnan(median_exo_cover) & median_exo_cover >= min_expr_level);
val = sort(median_exo_cover(idx),'ascend');
n = L;
bins = [];
while length(bins) < L;
  bins = round(linspace(1, length(val), n));
  % the following unique statement can result in a reduced number of
  % limits; that's why we increase N in the loop (but in most cases we
  % will just pass through the while loop once) 
  bins = unique(val(bins));
  n = n + 1;
end
bins = [[-inf; min_expr_level; bins(2:end-1)], [min_expr_level; bins(2:end-1); inf]];

% HACK if L=1
if (L==1),
    warning('Using only one expression level.');
    bins = [[-inf],[+inf]];
end

%warning('CHANGED EXPRESSION LEVEL CODE');
%bins = [-inf],
%gpl = floor(length(val)/L);
%for i=1:L-1,
%    curr = val(i*gpl);
%    bins = [bins; curr];
%end
%bins = [[bins], [bins(2:end); inf]];
assert(all(bins(1:end-1,1) < bins(2:end,1)));

% assign a discrete expression level to each gene
discrete_level = nan(N,1);
for i=1:L,
  idx = find(~isnan(median_exo_cover) ...
             & bins(i,1) <= median_exo_cover & median_exo_cover < bins(i,2));
  discrete_level(idx) = i;
end
assert(~any(isnan(discrete_level(~isnan(median_exo_cover)))));
assert(all(discrete_level(~isnan(discrete_level)) >= 1 ...
           & discrete_level(~isnan(discrete_level)) <= L));

for i=1:L
  fprintf('  level %i: read coverage between %2.1f and %2.1f (%i genes)\n', ...
          i, bins(i,1), bins(i,2), sum(discrete_level==i));
end
fprintf('\n');

% eof
