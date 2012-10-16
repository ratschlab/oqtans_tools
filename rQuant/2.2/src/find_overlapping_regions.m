function idx = find_overlapping_regions(regions1, regions2, consider_strands)
% FIND_OVERLAPPING_REGIONS   Finds regions that overlap in two sets of regions.
%
%   idx = find_overlapping_regions(regions1, regions2, consider_strands)
%
%   -- input --
%   regions1:         first set of regions (based on genes struct)
%   regions2:         second set of regions (based on genes struct);
%                     if empty then regions1 is compared against itself
%   consider_strands: if true then merging is done for each strand
%                     separately
%
%   -- output --
%   idx:              n x 2 matrix of indices corresponding to
%                     pairs of overlapping regions
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2010-2011 Jonas Behr, Gunnar Raetsch, Regina Bohnert
%   Copyright (C) 2010-2011 Max Planck Society
%


if isempty(regions1) && isempty(regions2)
  idx = [];
  return;
end

if nargin==2
  consider_strands = 1;
end

if isempty(regions2)
  % find overlapping pairs within one set of regions
  chr_nums = unique([regions1.chr_num]);
  all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums); % fast if there are many contigs
  idx = [];
  if consider_strands
    for chr = chr_nums,
      for s = '+-',
        chr_genes_idx1 = all_gene_chr_idxs1{chr}; 
        str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
        r1_chr_idx = chr_genes_idx1(str_genes_idx1);
        if isempty(r1_chr_idx), continue; end 
        idx = compute_overlap(regions1, r1_chr_idx, regions1, r1_chr_idx, idx);
      end
    end
  else
    for chr = chr_nums,
      r1_chr_idx = all_gene_chr_idxs1{chr}; 
      if isempty(r1_chr_idx), continue; end 
      idx = compute_overlap(regions1, r1_chr_idx, regions1, r1_chr_idx, idx);
    end
  end
  idx(idx(:,1)==idx(:,2),:) = [];
  idx = sort(idx')';
  idx = unique(idx, 'rows');
else
  % find overlapping pairs within two set of regions
  chr_nums = unique([[regions1.chr_num] [regions2.chr_num]]);
  all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums); % fast if there are many contigs
  all_gene_chr_idxs2 = store_gene_chr_idxs(regions2, chr_nums);
  idx = [];
  if consider_strands
    for chr = chr_nums,
      for s = '+-',
        chr_genes_idx1 = all_gene_chr_idxs1{chr}; % precomputed to make it faster for many contigs
        str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
        r1_chr_idx = chr_genes_idx1(str_genes_idx1);
        chr_genes_idx2 = all_gene_chr_idxs2{chr}; 
        str_genes_idx2 = find([regions2(chr_genes_idx2).strand]==s);
        r2_chr_idx = chr_genes_idx2(str_genes_idx2);
        if isempty(r1_chr_idx) || isempty(r2_chr_idx), continue; end 
        idx = compute_overlap(regions1, r1_chr_idx, regions2, r2_chr_idx, idx);
      end
    end
  else
    for chr = chr_nums,
      r1_chr_idx = all_gene_chr_idxs1{chr};
      r2_chr_idx = all_gene_chr_idxs2{chr};
      if isempty(r1_chr_idx) || isempty(r2_chr_idx) continue; end 
      idx = compute_overlap(regions1, r1_chr_idx, regions2, r2_chr_idx, idx);
    end
  end
end
return;


% pre-computes gene indices per chromosome
function all_idxs = store_gene_chr_idxs(genes, chr_nums)
max_chr_num = max(chr_nums);
all_idxs = cell(max_chr_num);
all_num = zeros(1, max_chr_num);
% allocate some mem
for n = 1:length(all_idxs),
  all_idxs{n} = zeros(1, 1000);
end
for g = 1:length(genes),
  chr_num = genes(g).chr_num;
  if all_num(chr_num)+1>length(all_idxs{chr_num})
    % double allocated mem
    all_idxs{chr_num} = [all_idxs{chr_num}, zeros(1, length(all_idxs{chr_num}))];
  end
  all_idxs{chr_num}(all_num(chr_num)+1) = g;
  all_num(chr_num) = all_num(chr_num)+1; 
end
all_genes_ids = [];
for n = 1:length(all_idxs),
  all_idxs{n}(all_num(n)+1:end) = [];
  all_genes_ids = [all_genes_ids, all_idxs{n}];
  assert(all(all_idxs{n}>0));
end
assert(length(all_genes_ids)==length(genes));
assert(length(all_genes_ids)==length(unique(all_genes_ids)));
return;


% computes overlaps between two regions
function idx = compute_overlap(regions1, r1_idx, regions2, r2_idx, idx)
starts1 = [regions1(r1_idx).start];
stops1  = [regions1(r1_idx).stop];
starts2 = [regions2(r2_idx).start];
stops2  = [regions2(r2_idx).stop];
[a b] = interval_overlap(starts1, stops1, starts2, stops2);
idx = [idx; [r1_idx(a)' r2_idx(b)']];
return;
