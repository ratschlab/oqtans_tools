function f = quad_fun(x, S1, S2, S3)
% QUAD_FUN   Computes value of a quadratic function.
%
%   f = quad_fun(x, S1, S2, S3)
%
%   -- input --
%   x:  variable (scalar)
%   S1: coefficient of quadratic term
%   S2: coefficient of linear term
%   S3: constant
%
%   -- output --
%   f:  function value
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

f = x^2*S1 + x*S2 + S3;
