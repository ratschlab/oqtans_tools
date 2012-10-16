function levels = get_expression_bins(genes, num_bins, method)
% Sort genes according to their respective (averaged) expression
% and discretize them into 'num_bins' bins. If 'method' equals one,
% the genes are not only discretized by expression levels, instead they
% are also equalized such that a equal amount of genes is in each bin.
%
% HINT:
% - a simple test showed that taking the max or the mean
%   during averaging doesn't make a big difference in the
%   outcome
%
% IN
%   genes     : [1..n] gene structures
%   num_bins  : number of bins
%
%
% OUT
%  levels     : [1..num_bins] cell array. each entry contains a list
%               with the associated gene indices.
%
% written by nico goernitz



% return structure
levels = cell(num_bins,1);

% average expression levels within each gene
fprintf('Averaging expression levels within each gene.\n');
max_diff = 0;
for i=1:length(genes),
    expr = genes(i).expr;
    foo = (max(expr)-min(expr));
    if (foo>max_diff), max_diff=foo; end;
    % for averaging various methods are possible:
    % take min, take max or mean, ..
    % however, for drosophila cufflinks data with 
    % max_diff>4500 the difference of taking max or mean is rather small
    genes(i).expr = max(genes(i).expr);
end
fprintf('Maximum difference of expression between transcripts of a certain gene: %1.2f\n',max_diff);


% sort all genes according to expression
expr_level = unique([genes.expr]);
lvl_inc = floor(length(expr_level)/num_bins);

% ..and equalize 
[foo, inds] = sort([genes.expr],'ascend');
gene_inc = floor(length(genes)/num_bins);

fprintf('%i genes with %i unique expression levels.\n',length(genes),length(expr_level));
for i=1:num_bins,
    
    if (method==1),
      % equalized expression levels
      gene_end = i*gene_inc;
      if (i==num_bins), gene_end=length(inds); end;
      levels{i} = inds((i-1)*gene_inc+1:gene_end);
    else
      % solely expression levels
      expr_start = expr_level((i-1)*lvl_inc + 1);
      expr_stop  = expr_level((i  )*lvl_inc);
      if (i==num_bins), expr_stop = expr_level(end); end;
      levels{i} = find([genes.expr]>=expr_start & [genes.expr]<=expr_stop);
    end

end


