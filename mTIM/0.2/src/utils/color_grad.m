function c = color_grad(N)

% c = color_grad([N])
%
% Returns a matrix of RGB values which is equal to the default colormap,
% but does not manipulate a figure object as a side-effect.
%
% Written by Georg Zeller, MPI Tuebingen, Germany, 2008

if exist('N', 'var'),
  assert(N <= 64);
end
c = zeros(64,3);
c(40:56,1) = 1;
c(24:40,2) = 1;
c( 8:24,3) = 1;
wedge = [1:16]' * 0.0625;
c(25:40,1) = wedge(1:end);
c(56:64,1) = wedge(end:-1:8);
c( 9:24,2) = wedge(1:end); 
c(40:55,2) = wedge(end:-1:1); 
c( 1: 8,3) = wedge(9:end);
c(24:39,3) = wedge(end:-1:1);

if exist('N', 'var'),
 c = c(round(linspace(1,64,N)),:);
end