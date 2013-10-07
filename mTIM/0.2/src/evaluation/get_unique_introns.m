function introns = get_unique_introns(genes)

cnt = 0;
for g=1:length(genes),
  for t=1:length(genes(g).transcripts),
    exons = genes(g).exons{t};
    cnt = cnt + size(exons,1) - 1;
  end
end
introns = zeros(cnt,4);

cnt = 0;
cnt2 = 0;
for g=1:length(genes),

  if (genes(g).strand=='.'),
      continue;
      cnt2=cnt2+1;
  end

  for t=1:length(genes(g).transcripts),
    exons = genes(g).exons{t};
    num_i = size(exons,1) - 1;
    s = +1;
    if genes(g).strand == '-',
      s=-1;
    end
    introns(cnt+1:cnt+num_i,1:2) = repmat([genes(g).chr_num, s], num_i, 1);
    introns(cnt+1:cnt+num_i,3:4) = [exons(1:end-1,2), exons(2:end,1)];
    cnt = cnt + num_i;
  end
end

t = size(introns,1);
introns = unique(introns, 'rows');
u = size(introns,1);
fprintf('%i total and %i unique introns\n', t, u);
fprintf('Discarded %i genes without strand information.\n', cnt2);
