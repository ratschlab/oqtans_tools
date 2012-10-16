function p_adj = get_adj_bins(CFG)
% GET_ADJ_BINS   Gets adjacent supporting points.
%
%   p_adj = get_adj_bins(CFG)
% 
%   -- input --
%   CFG:   configuration struct
%
%   -- output --
%   p_adj: F*N x 2 matrix of adjacent supporting points (1: previous, 2: next)
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

pw_nnz = get_included_thetas(CFG);
F = CFG.num_plifs;
N = size(CFG.transcript_len_ranges,1);

p_adj = zeros(2, F*N);
p_adj(1,:) = 0:F*N-1;
p_adj(1,1+F*[0:N-1]) = p_adj(1,2+F*[0:N-1]);
fidx = find((~pw_nnz))+1;
for f = 1:length(fidx),
  p_adj(1,fidx(f)) = p_adj(1,fidx(f)-1);
end
p_adj(2,:) = 2:F*N+1;
p_adj(2,F+F*[0:N-1]) = p_adj(2,F-1+F*[0:N-1]);
fidx = find(~pw_nnz)-1;
for f = length(fidx):-1:1,
  p_adj(2,fidx(f)) = p_adj(2,fidx(f)+1);
end