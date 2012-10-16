function perc = prctile(sample, k_percentiles)
% PRCTILE   Computes the k-th percentile of the given sample.
%
%   perc = prctile(sample, k_percentiles)
%
%   -- input --
%   sample:        a vector of which to compute the percentile; if
%                  sample is a matrix, percentiles will be computed 
%                  for each row vector
%   k_percentiles: vector of percentiles
%
%   -- output --
%   perc:          k-th percentiles of the given sample 
%                  (samples x percentiles)
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2007-2011 Georg Zeller, Regina Bohnert
%   Copyright (C) 2007-2011 Max Planck Society
%


perc = zeros(size(sample,1), length(k_percentiles));
for n = 1:size(sample,1),
  for m = 1:length(k_percentiles),
    x = sort(sample(n,:));
    num_x = length(sample(n,:));
    k = k_percentiles(m)/100 * (num_x-1) + 1;
    [k_int, D] = rat(k);
    if isequal(D,1),  % k is an integer, take percentile element
      perc(n, m) = x(k_int);
    else              % take midpoint between two elements
      perc(n, m) = mean(x(floor(k):ceil(k)));
    end
  end
end
