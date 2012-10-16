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

function pair_cover = get_pair_features_from_bam(CFG, region)
% pair_cover = get_pair_features_from_bam(CFG, region)
%
% -- input --
% CFG: configuration struct
% region: struct defining a region with start, stops, exons etc.
%
% -- output --
% pair_cover: vector of counts of read pairs aligning left and right at a given position

% check some of the input parameters
check_consistency = 1;
assert(length(region) == 1);
assert(region.start <= region.stop);

fname = CFG.read_map_file;
if ~exist(fname, 'file'),
  error('CFG.read_map_file points to nonexistent file %s', fname);
end
if exist('get_reads') ~= 3,
  error(['Could not find mex file get_reads.mex (Consider running make all' ...
         'in src/utils/get_reads/)!']);
end

% default return arguments
pair_cover = [];

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end
if CFG.both_strands
  strand = '0';
else
  strand = region.strand;
end

% parameters for generating intron spna features
if ~isfield(CFG, 'read_pair_cover_max_intron_len'),
  CFG.read_pair_cover_max_intron_len = 1e9;
end
if ~isfield(CFG, 'read_pair_cover_min_score'),
  CFG.read_pair_cover_min_score = 1;
end
if ~isfield(CFG, 'read_pair_cover_max_mismatches'),
  CFG.read_pair_cover_max_mismatches = 5;
end

% We need to ask for a larger region to include cropped reads
max_read_len = 250;
win_left  = min(max_read_len, region.start-1);
win_right = min(max_read_len, CFG.chr_lens(region.chr_num)-region.stop);
assert(win_left  >= 0);
assert(win_right >= 0);

chr = CFG.chr_names{region.chr_num};
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
return_align_len = 0;
% flag to get read-pair information returned
return_pairs = 1;

% get read-pair information from bam file
[ex, in, pair_cover] = get_reads(fname, chr, b, e, strand, collapse, subsample, ...
                                 CFG.read_pair_cover_max_intron_len, CFG.read_pair_cover_min_score, ...
                                 CFG.read_pair_cover_max_mismatches, mapped, spliced, ...
                                 return_align_len, return_pairs);

assert(size(pair_cover, 1) == 1);

% crop off win_left and win_right from exon coverage and intron span
pair_cover(end-win_right+1: end) = [];
pair_cover(1:win_left)           = [];
assert(size(pair_cover,2) == region.stop - region.start + 1);

% eof