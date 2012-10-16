function limits = get_limits(max_lmt, num_lmt)
% GET_LIMITS   Sets supporting points of PLiF.
%
%   limits = get_limits(max_lmt, num_lmt)
%
%   -- input --
%   max_lmt: maximal supporting point of PLiF
%   num_lmt: number of supporting points
%
%   -- ouput --
%   limits:  supporting points of PLiF
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


limits = round(linspace(0, max_lmt, num_lmt));
%limits = round(linspace(0, sqrt(max_lmt), num_lmt).^2);
limits(1) = 1;
limits(end) = inf;
assert(isequal(limits, unique(limits)));