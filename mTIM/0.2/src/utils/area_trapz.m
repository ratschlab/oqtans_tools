function area = area_trapz(x, y)

% area = area_trapz(x, y)
%
% Calculates the area under the piecewise linear curve specified by a
% vector of x-coordinates and y-coordinates using a trapezoid
% approximation.
%
% x -- a vector of x-coordinates of the piecewise linear function
% y -- a vector of y-coordinates of the piecewise linear function
% returns the area under the curve
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2008

assert(length(x)==length(y));
assert(all(x~=-inf) & all(y~=inf));

d_x = x(2:end)-x(1:end-1);
assert(all(d_x >= 0));

y_max = max([y(2:end); y(1:end-1)]);
y_min = min([y(2:end); y(1:end-1)]);
area = (sum(y_max.*d_x) + sum(y_min.*d_x))./2;
