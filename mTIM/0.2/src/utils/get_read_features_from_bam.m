%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch, Georg Zeller,
% Andre Kahles, Pramod Mudrakarta
% Copyright (C) 2009-2011 Max Planck Society
%

function [exon_cover intron_span intron_list] = get_read_features_from_bam(CFG, region)
% [exon_cover intron_span intron_list] = get_read_features_from_bam(CFG, region)
%
% -- input --
% CFG: configuration struct
% region: struct defining a region with start, stops, exons etc.
%
% -- output --
% exon_cover: vector of read counts aligning at a given position
% intron_span: vector of counts of spliced alignments spanning over a given position
% intron_list: nx5 list of introns (intron start, intron stop, strand,
% number of reads, spliced alignment score)

% check some of the input parameters
check_consistency = 0;
assert(length(region) == 1);
assert(region.start <= region.stop);

fname = CFG.read_map_file;
if ~exist(fname, 'file'),
  error('CFG.read_map_file points to nonexistent file %s', fname);
end
if exist('get_reads') ~= 3,
  warning(['Could not find mex file get_reads.mex (Consider running make all' ...
           'in src/utils/get_reads/)!']);
end

% default return arguments
exon_cover = [];
intron_list = zeros(4, 0);

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end
if CFG.both_strands
  strand = '0';
else
  strand = region.strand;
end

% parameters for generating intron spna features
if ~isfield(CFG, 'read_intron_span_max_intron_len'),
  CFG.read_intron_span_max_intron_len = 1e9;
end

if ~isfield(CFG, 'read_intron_span_min_score'),
  CFG.read_intron_span_min_score = 25;
end
if ~isfield(CFG, 'read_intron_span_max_mismatches'),
  CFG.read_intron_span_max_mismatches = 0;
end

% parameters for extracting splice sites from spliced read alignments
if ~isfield(CFG, 'read_splice site_max_intron_len'),
  CFG.read_splice_site_max_intron_len = 1e9;
  assert(CFG.read_splice_site_max_intron_len >= CFG.read_intron_span_max_intron_len);
end
if ~isfield(CFG, 'read_splice_site_min_score'),
  CFG.read_splice_site_min_score = 8;
  assert(CFG.read_splice_site_min_score <= CFG.read_intron_span_min_score);
end
if ~isfield(CFG, 'read_splice_site_max_mismatches'),
  CFG.read_splice_site_max_mismatches = 0;
  assert(CFG.read_splice_site_max_mismatches == CFG.read_intron_span_max_mismatches);
end

% We need to ask for a larger region to include cropped reads
max_read_len = 250;
win_left  = min(max_read_len, region.start-1);
win_right = min(max_read_len, CFG.chr_lens(region.chr_num)-region.stop);
assert(win_left  >= 0);
assert(win_right >= 0);

c = region.chr_num;
b = region.start - win_left;
e = region.stop  + win_right;

% collapse individual read alignments to position-wise coverage
collapse = 1;
% do not subsample reads
subsample = 1000;
% flags to get_reads function to specify which alignments to include for
% coverage calculations
mapped = 1;
spliced = 1;
% flag to return the length of partial exonic read alignments for spliced
% alignments as part of the intron_list (3rd & 4th rows)
return_align_len = 1;

% get read alignments from bam file
[exon_cover, intron_list] = get_reads(fname, CFG.chr_names{c}, b, e, strand, collapse, subsample, ...
                                      CFG.read_splice_site_max_intron_len, CFG.read_splice_site_min_score, ...
                                      CFG.read_splice_site_max_mismatches, mapped, spliced, return_align_len);

%%% we don't need do process exon coverage
assert(size(exon_cover, 1) == 1);

%%% process intron list (to obtain the following format: 1: intron start, 2:
% intron stop, 3: strand, 4: support, 5:align_score) support will be the
% number of read alignments supporting a splice junction, align_score will
% give the length of the read alignment on either side of the junction
% (precisely the maximum over all reads of the minimum over both sides)
% moreover, extract intron span, a position-wise counts of spliced
% alignments with an intron gap at the given position
if nargout>=2,
  % because of the region extension, it can happen that returned introns
  % are completely outside the region of interest, hence we remove those first 
  keep = find(intron_list(2,:) >= region.start & intron_list(1,:) <= region.stop);
  intron_list = intron_list(:,keep);

  % aggregate scores of individual spliced alignments in terms of unique introns
  algn_overhang = intron_list(3:4,:)';
  intron_list = [intron_list(1:2,:)', ...
                 (region.strand=='-')*ones(size(intron_list,2), 1), ...
                 zeros(size(intron_list,2), 2)];
  intron_list_unique = unique(intron_list, 'rows');

  % count the number of read alignments confirming the splice junction
  for n=1:size(intron_list_unique,1),
    intron_list_unique(n,4) = sum(intron_list_unique(n,1)==intron_list(:,1) & ...
                                  intron_list_unique(n,2)==intron_list(:,2) & ...
                                  intron_list_unique(n,3)==intron_list(:,3));
  end

  % take the maximum over all minimal alignment segments supporting the intron
  for n=1:size(intron_list_unique, 1),
    match_idx = intron_list_unique(n,1)==intron_list(:,1) & ...
        intron_list_unique(n,2)==intron_list(:,2) & ...
        intron_list_unique(n,3)==intron_list(:,3);
    intron_list_unique(n,5) = max(min(algn_overhang(match_idx,:), [], 2));
  end 

  % intron span (position-wise count of the number of spliced alignments
  % spanning over the given position)
  intron_span = zeros(size(exon_cover), 'int32');
  take_idx = find(intron_list_unique(:,2)-intron_list_unique(:,1)+1 ...
                  <= CFG.read_intron_span_max_intron_len ...
                  & intron_list_unique(:,5) >= CFG.read_intron_span_min_score);
  for n=1:length(take_idx),
    i = take_idx(n);
    idx_b = max(1, intron_list_unique(i,1) - b + 1);
    idx_e = min(size(intron_span,2), intron_list_unique(i,2) - b + 1);
    intron_span(idx_b:idx_e) = intron_span(idx_b:idx_e) + intron_list_unique(i,4); 
  end

  if check_consistency,
    ino_blocks = find_blocks(intron_span>0) + b - 1;
    if ~isempty(ino_blocks),
      mem = ismember(ino_blocks(1,:), intron_list_unique(:,1));
      assert(all(mem(2:end-1)));
      mem = ismember(ino_blocks(2,:), intron_list_unique(:,2));
      assert(all(mem(2:end-1)));
    end
  end

  % crop off introns extending over the given region (necessary as we queried
  % a region extended by win_left and win_right)
  crop_idx = find(intron_list_unique(:,1) < region.start ...
                  | intron_list_unique(:,2) > region.stop);
  intron_list_unique(crop_idx,:) = [];
  intron_list = intron_list_unique;
end

% crop off win_left and win_right from exon coverage and intron span
exon_cover(end-win_right+1: end) = [];
exon_cover(1:win_left)           = [];
assert(size(exon_cover,2) == region.stop - region.start + 1);
if nargout>=2,
  intron_span(end-win_right+1: end) = [];
  intron_span(1:win_left)           = [];
  assert(size(intron_span,2) == region.stop - region.start + 1);
end

% eof