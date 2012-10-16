function [coverage reads_ok intron_list read_starts paired_reads] = get_coverage_per_read(CFG, gene, reverse_ret)
% GET_COVERAGE_PER_READ   Gets the reads from the BAM file covering the gene region.
%
%   [coverage reads_ok intron_list read_starts paired_reads] = get_coverage_per_read(CFG, gene, reverse_ret)
%
%   -- input --
%   CFG:            configuration struct
%   gene:           struct defining a gene with start, stops, exons etc.
%   reverse_ret:    if true, the positions are considered in reverse
%                   direction on the reverse strand
%
%   -- output --
%   coverage:       matrix of exonic positions x reads
%   reads_ok:       indicates success of file parsing
%   intron_list:    nx4 list of introns
%                   (intron start, intron stop, confirmation, strand)
%   read_starts:    vector of number of reads starting at the given exonic positions
%   paired_reads:   struct including starts, stops and mates for paired-end reads
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


% initialisation
coverage = [];
reads_ok = 1;
intron_list = zeros(2, 0);
read_starts = zeros(0, 1);
paired_reads.starts = zeros(1, 0);
paired_reads.stops = zeros(1, 0);
paired_reads.mates = zeros(2, 0);

collapse = double(~CFG.paired & nargout<4);

if nargin<3
  reverse_ret = 0;
end

if nargout==5
  assert(CFG.paired==1);
end

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end

if CFG.both_strands
  strand = '0';
else
  strand = gene.strand;
end

% parameters for get_reads
if ~isfield(CFG, 'tracks_max_intron_len')
  CFG.tracks_max_intron_len = 1e9;
end
if ~isfield(CFG, 'tracks_min_exon_len')
  CFG.tracks_min_exon_len = -1;
end
if ~isfield(CFG, 'tracks_max_mismatches')
  CFG.tracks_max_mismatches = CFG.read_len;
end
subsample = 1000; mapped = 1; spliced = 1; maxminlen = 0;


win = CFG.read_len;
eidx = [max(gene.eidx(1)-win,1):gene.eidx(1)-1, gene.eidx, gene.eidx(end)+1:gene.eidx(end)+win];
win_size = length(max(gene.eidx(1)-win,1):gene.eidx(1)-1);

for f = 1:length(CFG.tracks_fn{gene.chr_num}),
  fname = CFG.tracks_fn{gene.chr_num}{f};
  if ~exist(fname, 'file'),
    warning('BAM file %s does not exist', fname);
  end
  try
    if nargout==5
      [coverage_idx_tmp{f}, intron_list_tmp, cov_paired_tmp, mates_tmp] = get_reads(fname, gene.chr, eidx(1), eidx(end), strand, collapse, subsample, CFG.tracks_max_intron_len, CFG.tracks_min_exon_len, CFG.tracks_max_mismatches, mapped, spliced, maxminlen, CFG.paired);
      clear cov_paired_tmp;
    elseif nargout==3 || nargout==4
      [coverage_idx_tmp{f}, intron_list_tmp] = get_reads(fname, gene.chr, eidx(1), eidx(end), strand, collapse, subsample, CFG.tracks_max_intron_len, CFG.tracks_min_exon_len, CFG.tracks_max_mismatches, mapped, spliced, maxminlen, CFG.paired);
    else
      [coverage_idx_tmp{f}] = get_reads(fname, gene.chr, eidx(1), eidx(end), strand, collapse, subsample, CFG.tracks_max_intron_len, CFG.tracks_min_exon_len, CFG.tracks_max_mismatches, mapped, spliced, maxminlen, CFG.paired);
    end
  catch
    warning('get_reads failed');
    intron_list = intron_list';
    reads_ok = 0;
    return;
  end
  if exist('intron_list_tmp', 'var')
    if ~isempty(intron_list_tmp),
      intron_list = [intron_list intron_list_tmp];
    end
  end
  if exist('mates_tmp', 'var')
    if ~isempty(mates_tmp),
      paired_reads.mates = [paired_reads.mates mates_tmp+1];
    end
  end
end

% process coverage: convert to exonic position indices
coverage_idx = [coverage_idx_tmp{:}];
if ~collapse
  coverage_idx = unique(coverage_idx', 'rows')'; % no overlapping reads
  if ~isempty(coverage_idx)
    coverage = sparse(coverage_idx(1,:), coverage_idx(2,:), 1, max(coverage_idx(1,:)), eidx(end)-eidx(1)+1)';
  else
    coverage = sparse([], [], 1, eidx(end)-eidx(1)+1, 0);
  end
  coverage = coverage(eidx(win_size+1:win_size+gene.exonic_len)-eidx(1)+1, :);
  % no overlapping reads
  assert(~any(any(full(coverage>1))));
else
  coverage = sum(coverage_idx, 1);
  coverage = sparse(coverage(eidx(win_size+1:win_size+gene.exonic_len)-eidx(1)+1)');
end
  
% process intron list (1: intron start, 2: intron stop, 3: confirmation, 4: strand)
if nargout>2 && ~isempty(intron_list)
  intron_list = [intron_list', zeros(size(intron_list,2), 1), (gene.strand=='-')*ones(size(intron_list,2), 1)];
  intron_list_unique = unique(intron_list, 'rows');
  for n = 1:size(intron_list_unique,1),
    intron_list_unique(n,3) = sum(intron_list_unique(n,1)==intron_list(:,1) & ...
                                  intron_list_unique(n,2)==intron_list(:,2) & ...
                                  intron_list_unique(n,4)==intron_list(:,4));
  end
  intron_list = intron_list_unique;
  clear intron_list_unique;
end

% process read starts: count reads starting at each exonic position
if nargout>3
  read_starts = zeros(gene.exonic_len, 1);
  for c = 1:size(coverage, 2),
    fidx = find(coverage(:,c)~=0, 1, 'first');
    if ~isempty(fidx),
      if fidx==1 && sum(coverage(:,c),1)<CFG.read_len, continue; end
      read_starts(fidx) = read_starts(fidx) + 1;
    end
  end
end

% process read starts and stops: store indices (eidx based)
if nargout>4
  paired_reads.starts = nan(1, size(coverage, 2));
  paired_reads.stops = nan(1, size(coverage, 2));
  for c = 1:size(coverage, 2),
    fidx = find(coverage(:,c)~=0, 1, 'first');
    lidx = find(coverage(:,c)~=0, 1, 'last');
    if ~isempty(fidx) && ~isempty(lidx)
      if (fidx==1 || lidx==size(coverage, 1)) && sum(coverage(:,c),1)<CFG.read_len, continue; end
      paired_reads.starts(c) = gene.eidx(fidx);
      paired_reads.stops(c) = gene.eidx(lidx);
    end
  end
  idx = find(~isnan(paired_reads.starts) & ~isnan(paired_reads.stops));
  [C fidx1] = intersect(paired_reads.mates(1,:), idx);
  [C fidx2] = intersect(paired_reads.mates(2,:), idx);
  paired_reads.mates = paired_reads.mates(:,intersect(fidx1, fidx2));
end

% collapse coverage as single reads are not required anymore
if CFG.paired || nargout>3
  coverage = sum(coverage, 2);
end

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = size(coverage,1):-1:1;
  coverage = coverage(rev_idx,:);
end
