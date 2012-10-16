function ins_sizes = estimate_ins_size_distr(CFG, genes)
% ESTIMATE_INS_SIZE_DISTR   Determines insert size distribution from observed paired-end reads. 
%
%   ins_sizes = estimate_ins_size_distr(CFG, genes)
%
%   -- input --
%   CFG:       configuration struct
%   genes:     struct defining genes with start, stops, exons etc.
%
%   -- output --
%   ins_sizes: vector of insert sizes from observed paired-end reads
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


ins_sizes = [];
for g = 1:length(genes),
  gene = genes(g);
  % get read data
  [coverage, reads_ok, introns, read_starts, paired_reads] = get_coverage_per_read(CFG, gene);
  clear read_starts;
  % estimate insert size distribution from observed reads
  if gene.is_alt==0
    % gene with one transcript: take reads from whole transcript
    offset = gene.start-1;
    gidx = sparse(ones(1,length(gene.eidx)), gene.eidx-offset, 1:length(gene.eidx),  1, gene.stop-gene.start+1);
    r1 = paired_reads.mates(1,:);
    r2 = paired_reads.mates(2,:);
    fidx = find(paired_reads.starts(r1)>paired_reads.starts(r2) & paired_reads.stops(r1)>paired_reads.stops(r2));
    if ~isempty(fidx)
      r = r1(fidx); r1(fidx) = r2(fidx); r2(fidx) = r;
    end
    cand_starts = paired_reads.stops(r1);
    cand_stops = paired_reads.starts(r2);
    idx1 = gidx(cand_starts-offset);
    idx2 = gidx(cand_stops-offset);
    ins_sizes = [ins_sizes idx2-idx1-1];
  else
    % gene with alternative transcripts: take reads only within a distinguishable segment
    segments = define_segments(gene.splicegraph{1}, gene.splicegraph{2}); % distinguishable segments
    for p = 1:size(paired_reads.mates, 2),
      r1 = paired_reads.mates(1,p);
      r2 = paired_reads.mates(2,p);
      if paired_reads.starts(r1)>paired_reads.starts(r2) && paired_reads.stops(r1)>paired_reads.stops(r2)
        r = r1; r1 = r2; r2 = r;
      end
      s1 = find(segments(:,1)<paired_reads.starts(r1) & segments(:,2)>=paired_reads.starts(r1));
      s2 = find(segments(:,1)<paired_reads.stops(r2) & segments(:,2)>=paired_reads.stops(r2));
      if s1==s2
        ins_sizes = [ins_sizes paired_reads.starts(r2)-paired_reads.stops(r1)-1];
      end
    end
  end
end
assert(all(ins_sizes>=0));