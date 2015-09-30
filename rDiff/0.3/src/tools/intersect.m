%## Copyright (C) 2000-2012 Paul Kienzle
%## Copyright (C) 2008-2009 Jaroslav Hajek
%##
%## This file is part of Octave.
%##
%## Octave is free software; you can redistribute it and/or modify it
%## under the terms of the GNU General Public License as published by
%## the Free Software Foundation; either version 3 of the License, or (at
%## your option) any later version.
%##
%## Octave is distributed in the hope that it will be useful, but
%## WITHOUT ANY WARRANTY; without even the implied warranty of
%## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ...
%    GNU
%## General Public License for more details.
%##
%## You should have received a copy of the GNU General Public ...
%    License
%## along with Octave; see the file COPYING.  If not, see
%## <http://www.gnu.org/licenses/>.

function [c, ia, ib] = intersect (a, b, varargin)

if (nargin < 2 || nargin > 3)
  print_usage ();
end


if (isempty (a) || isempty (b))
  c = ia = ib = [];
else
  ## form a and b into sets
  if (nargout > 1)
    [a, ja] = unique (a, varargin{:});
    [b, jb] = unique (b, varargin{:});
  else
    a = unique (a, varargin{:});
    b = unique (b, varargin{:});
  end
  
  if (nargin > 2)
    c = [a; b];
    [c, ic] = sortrows (c);
    ii = find (all (c(1:end-1,:) == c(2:end,:), 2));
    c = c(ii,:);
    len_a = rows (a);
  else
    c = [a(:); b(:)];
    [c, ic] = sort (c);               ## [a(:);b(:)](ic) == c
    if (iscellstr (c))
      ii = find (strcmp (c(1:end-1), c(2:end)));
    else
      ii = find (c(1:end-1) == c(2:end));
    end
    c = c(ii);
    len_a = length (a);
  end
  
  if (nargout > 1)
    ia = ja(ic(ii));                  ## a(ia) == c
    ib = jb(ic(ii+1) - len_a);        ## b(ib) == c
  end
  
  if (nargin == 2 && (size (b, 1) == 1 || size (a, 1) == 1))
    c = c.';
  end
end
