function [feat, idx, idx_next, del_idx] = gen_exon_features(CFG, gene, t, reverse_ret)
% GEN_EXON_FEATURES   Generates transcript length feature matrix for PliFs.
%
%   [feat, idx, idx_next, del_idx] = gen_exon_features(CFG, gene, t, reverse_ret)
%
%   -- input --
%   CFG:         configuration struct with
%                   num_plifs: number of supporting points for PLiFs and
%                   max_side_len: maximal number of positions to be
%                   considered at both transcript parts
%   gene:        struct defining a gene with start, stops, exons etc.
%   t:           index of transcript
%   reverse_ret: if true, the positions are considered in reverse
%                direction on the reverse strand
%
%   -- output --
%   feat:        vector of features for P exonic positions
%   idx:         vector of indices to supporting points f
%   idx_next:    vector of indices to supporting points f+1
%   del_idx:     vector of position that do not cover a whole bin
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


if nargin<4  reverse_ret = 0;
end

offset = gene.start-1;
exons = gene.exons{t};

% transcript indices in relative gene coordinates for all exons
eidx = gene.eidx;
eidx = unique(eidx-offset);

% transcript indices in relative gene coordinates
tidx = [];
for e = 1:size(exons,1),
  tidx = [tidx, exons(e,1):exons(e,2)];
  tidx = unique(tidx);
end
tidx = unique(tidx)-offset;
assert(all(tidx>=0));

% transcript indices in relative exonic coordinates
[tmp idx1 idx2] = intersect(tidx, eidx);
assert(isequal(tidx, eidx(idx2)));
tidx = idx2;

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = length(eidx):-1:1;
  tidx = sort(rev_idx(tidx));
end

eidx2tidx = zeros(1, length(eidx));
eidx2tidx(tidx) = 1:length(tidx);

feat = zeros(length(eidx), 1);
idx = zeros(length(eidx), 1);
num_bins = CFG.num_plifs/2 - 1;
lmt = get_limits(CFG.max_side_len, num_bins+1);

% left transcript part
dist = 1:ceil(length(tidx)*0.5);
for b = 1:num_bins,
  fidx = find(lmt(b)<=dist & lmt(b+1)>dist);
  assert(all(feat(tidx(fidx),1)==0));
  if b==num_bins
    feat(tidx(fidx),1) = 0;
  else
    feat(tidx(fidx),1) = (dist(fidx)-lmt(b))./(lmt(b+1)-lmt(b));
  end
  idx(tidx(fidx),1) = b;
end

% right transcript part
dist = 1:(length(tidx)-length(tidx)*0.5);
for b = 1:num_bins,
  fidx = find(lmt(b)<=dist & lmt(b+1)>dist);
  assert(all(feat(tidx(length(tidx)-fidx+1),1)==0));
  if b==num_bins
    feat(tidx(length(tidx)-fidx+1),1) = 0;
  else
    feat(tidx(length(tidx)-fidx+1),1) = (dist(fidx)-lmt(b+1))./(lmt(b)-lmt(b+1));
  end
  idx(tidx(length(tidx)-fidx+1),1) = num_bins*2+1-b+1;
end

% correct features for bins that are only partly covered by transcript positions
midx = ceil(length(tidx)*0.5);
bin1 = idx(tidx(midx));
bin2 = idx(tidx(midx+1));
idx1 = eidx2tidx(find(idx==bin1, 1, 'first'));
idx2 = eidx2tidx(find(idx==bin2, 1, 'last'));
%feat(idx1:idx2) = mean([feat(idx1:idx2), (([idx1:idx2]-idx1)./(idx2-idx1))'], 2);
feat(tidx(idx1:idx2)) = ([idx1:idx2]-idx1)./(idx2-idx1);
del_idx = tidx(idx1:idx2);
idx(del_idx) = bin1;
if nargout>=3,
  profile_weights_adj = get_adj_bins(CFG);
  idx_next = zeros(length(eidx), 1);
  tmp_pw = profile_weights_adj(2,idx(tidx)+CFG.num_plifs*(gene.transcript_len_bin(t)-1));
  idx_next(tidx) = tmp_pw - CFG.num_plifs*(gene.transcript_len_bin(t)-1);
  idx_next(del_idx) = bin2+1;
end

assert(all(sum(feat<=1 & feat>0, 2)==1 | sum(feat,2)==0));
assert(size(feat,1)==length(eidx));