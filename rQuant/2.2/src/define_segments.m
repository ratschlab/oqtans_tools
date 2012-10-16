function [segments, exon_pointer, seg_admat, initial, terminal] = define_segments(exons, admat)
% DEFINES_SEGMENTS   Generates distinguishable segments from splicegraph.
%
%   [segments, exon_pointer, seg_admat, initial, terminal] = define_segments(exons, admat)
%
%   -- input --
%   exons:        exons from splicegraph
%   admat:        adjacency matrix from splicegraph
%
%   -- output --
%   segments:     distinguishable segments
%   exon_pointer: pointer to exons from which segments were generated
%   seg_admat:    adjacency matrix for segments
%   initial:      indicates whether segments is in initial exon
%   terminal:     indicates whether segments is in terminal exon
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Jonas Behr
%   Copyright (C) 2011 Max Planck Society
%


% map of the number of exons covering a position
s = min(exons(1, :));
e = max(exons(2, :));
exons = exons -s+1;

maps = zeros(size(exons, 2), e-s+1);
for j = 1:size(exons, 2)
  exon_start = exons(1, j);
  exon_stop = exons(2, j);
  maps(j, exon_start:exon_stop) = j;
end

% map = sum(maps, 1);

segments = [];
exon_pointer = {};
coverage = [];

num_exons = maps(:,1);
last_start = 1;
k = 0;
while k<size(maps, 2)
  k = k+1;
  if ~isequal(maps(:, k), num_exons)
    segments = [segments, [last_start;k-1]];
    exon_pointer{end+1} = setdiff(unique(maps(:, k-1)), 0);
    %coverage = [coverage, map(k-1)]; 
    while all(maps(:,k)==0) % skip over introns
      k = k+1;
    end
    last_start = k;
    num_exons = maps(:,k);
  end
end
segments = [segments, [last_start;k]];
exon_pointer{end+1} = setdiff(unique(maps(:, k)), 0);
%coverage = [coverage, map(k)];

% built splicegraph based on segments
% there are several things I need to encode here: 
% 1. connection is forbidden: 		-2
% 2. connecton is within one exon: 	-1
% 3. connection is valid intron: 	>=0
num_seg = size(segments, 2);
seg_admat = -2*ones(num_seg);
for j = 1:num_seg
  for k = j+1:num_seg
    seg_end = segments(2, j);
    seg_start = segments(1, k);
    
    if k==j+1 && seg_end == seg_start-1 && ~isempty(intersect(exon_pointer{j}, exon_pointer{k}))
      % if they dont come from the same exon, they come from different transcripts that 
      % exactly touch with their ends
      seg_admat(j, k) = -1;
      seg_admat(k, j) = -1;
      continue
    end
    for e1 = exon_pointer{j}'
      has_intron = 0;
      if exons(2, e1)~=seg_end
        % this segment is not in the end on this exon and therefore it
        % may not share its connectivity
        continue
      end
      for e2 = exon_pointer{k}'
        if exons(1, e2)~=seg_start
          continue
        end
        if admat(e1, e2)>0 || admat(e2, e1)>0 % admat should be symmetric, but this is sometimes violated
          has_intron = 1;
        end
      end
      if has_intron
        seg_admat(j, k) = 0;
        seg_admat(k, j) = 0;
      end
    end
  end	
end

%	exons can be terminal exons, but still oberlapp with other exons. 
%	Then the termination information is lost when reducing the exons to segments 
%
if nargout>3
  terminal = zeros(1, num_seg);
  terminal(end) = 1;
  initial = zeros(1, num_seg);
  initial(1) = 1;
  for j = 1:num_seg
    seg_end = segments(2, j);
    seg_start = segments(1, j);
    %if j==5
    %	keyboard
    %end
    for e1 = exon_pointer{j}'
      if exons(2, e1)==seg_end
        terminal(j) = terminal(j) || is_terminal(exons, admat, e1);
      end
      if exons(1, e1)==seg_start
        initial(j) = initial(j) || is_initial(exons, admat, e1);
      end
    end	
  end
end

segments = segments' + s-1;
return

function ret = is_terminal(exons, admat, idx)
% get all connected nodes
% this in only needed if the exons 
% in the splicegraph are not sorted
cnodes = find(admat(idx,:)>0);	
for j = cnodes
  % check that it is really a child node
  if exons(1, j)>exons(2, idx) 
    ret = 0;
    return
  end
end
ret = 1;
return

function ret = is_initial(exons, admat, idx)
% get all connected nodes
% this in only needed if the exons 
% in the splicegraph are not sorted
cnodes = find(admat(idx,:)>0);	
for j = cnodes
  % check that it is really a parent node
  if exons(2, j)<exons(1, idx) 
    ret = 0;
    return
  end
end
ret = 1;
return

