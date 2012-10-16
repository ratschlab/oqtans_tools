function genes = label_to_genes(label_seq, offset)

LABELS = get_label_set_mTIM();

gene_blocks = find_blocks(label_seq ~= LABELS.intergenic);

g = 0;
while g < size(gene_blocks,2),
  g = g + 1;
  idx = gene_blocks(1,g):gene_blocks(2,g);
  exw_blocks = find_blocks(label_seq(idx) == LABELS.exon_W) ...
      + gene_blocks(1,g) - 1;
  exc_blocks = find_blocks(label_seq(idx) == LABELS.exon_C) ...
      + gene_blocks(1,g) - 1;
  
  % split again if exons from both strands are present
  while ~isempty(exw_blocks) && ~isempty(exc_blocks),
    if exw_blocks(1,1) < exc_blocks(1,1),
      e = exc_blocks(1,1) - 1;
      idx = find(exw_blocks(2,:)<=e, 1, 'last');
      b = exw_blocks(2,idx) + 1;
      exw_blocks(:,1:idx) = [];
    else
      e = exw_blocks(1,1) - 1;
      idx = find(exc_blocks(2,:)<=e, 1, 'last');
      b = exc_blocks(2,idx) + 1;
      exc_blocks(:,1:idx) = [];
    end
    split = find(label_seq(b:e) == LABELS.intergenic ...
                 | label_seq(b:e) == LABELS.ambiguous) + b - 1;
    split = split(round(length(split)/2));
    gene_blocks(:,g+2:end+1) = gene_blocks(:,g+1:end);
    gene_blocks(2,g+1) = gene_blocks(2,g);
    gene_blocks(2,g) = split-1;
    gene_blocks(1,g+1) = split;
    g = g + 1;
  end
end

if ~isempty(gene_blocks), % this test is necessary to avoid a seg fault (!!!)
  assert(issorted(gene_blocks', 'rows'));
end

genes = [];
cnt = 1;
for g=1:size(gene_blocks,2),
  genes(cnt).start = gene_blocks(1,g) + offset;
  genes(cnt).stop = gene_blocks(2,g) + offset;

  gene_offset = offset + gene_blocks(1,g) - 1;

  idx = gene_blocks(1,g):gene_blocks(2,g);
  exw_idx = label_seq(idx) == LABELS.exon_W;
  exc_idx = label_seq(idx) == LABELS.exon_C;
  exons_W = find_blocks(exw_idx)' + gene_offset;
  exons_C = find_blocks(exc_idx)' + gene_offset;
  exw_idx = idx(exw_idx);
  exc_idx = idx(exc_idx);
  
  if ~isempty(exons_W),
    assert(isempty(exons_C));
    genes(cnt).strand = '+';
    %genes(cnt).exons{1} = exons_W
    
    genes(cnt).exons = {};
    tmp = genes(cnt).exons;
    tmp{1} = exons_W;
    genes(cnt).exons = tmp;
  else
    genes(cnt).strand = '-';
   %genes(cnt).exons{1} = exons_C
    genes(cnt).exons = {};
    tmp = genes(cnt).exons;
    tmp{1} = exons_C;
    genes(cnt).exons = tmp;
  end
  cnt = cnt + 1;
end
