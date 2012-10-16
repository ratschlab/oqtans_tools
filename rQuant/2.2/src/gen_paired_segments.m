function pair_mat = gen_paired_segments(segments, paired_reads)
% GEN_PAIRED_SEGMENTS   Generates matrix of segments connected by paired-end reads.
%
%   pair_mat = gen_paired_segments(segments, paired_reads)
%
%   -- input --
%   segments:     distinguishable segments
%   paired_reads: struct including starts, stops and mates for paired-end reads
%
%   -- output --
%   paired_mat:   matrix (#segments x #segments) of connectivity by paired-end reads
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert, Jonas Behr
%   Copyright (C) 2011 Max Planck Society
%


pair_mat = zeros(size(segments, 1));
for p = 1:size(paired_reads.mates, 2),
  r1 = paired_reads.mates(1,p);
  r2 = paired_reads.mates(2,p);
  if paired_reads.starts(r1)>paired_reads.starts(r2) && paired_reads.stops(r1)>paired_reads.stops(r2)
    r = r1; r1 = r2; r2 = r;
  end
  s0 = find(segments(:,1)<paired_reads.starts(r1) & segments(:,2)>=paired_reads.starts(r1));
  s1 = find(segments(:,1)<paired_reads.stops(r1) & segments(:,2)>=paired_reads.stops(r1));
  s2 = find(segments(:,1)<paired_reads.starts(r2) & segments(:,2)>=paired_reads.starts(r2));
  s3 = find(segments(:,1)<paired_reads.stops(r2) & segments(:,2)>=paired_reads.stops(r2));
  if 1
  % find all segments connected and covered by a read pair
  sl = [];
  if ~isempty(s0), sl = s0; end
  if ~isempty(s1), sl = [sl s1]; end
  if ~isempty(sl), sl = unique(sl); end
  sr = [];
  if ~isempty(s2), sr = s2; end
  if ~isempty(s3), sr = [sr s3]; end
  if ~isempty(sr), sr = unique(sr); end
  for n1 = sl,
    for n2 = sr,
      pair_mat(n1,n2) = pair_mat(n1,n2) + 1;
      if n1~=n2
        pair_mat(n2,n1) = pair_mat(n2,n1) + 1;
      end
    end
  end
  else
  % find only segments connected of the insert size of the read pair
  if length(s1)==1 && length(s2)==1
    pair_mat(s1,s2) = pair_mat(s1,s2) + 1;
    if s1~=s2
      pair_mat(s2,s1) = pair_mat(s2,s1) + 1;
    end
  else
    assert(length(s1)<=1);
    assert(length(s2)<=1);
  end
  end
end