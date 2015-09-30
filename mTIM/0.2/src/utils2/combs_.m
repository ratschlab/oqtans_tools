function P = combs(v,m)
%COMBS_  All possible combinations.
%   COMBS_(1:N,M) or COMBS_(V,M) where V is a row vector of length N,
%   creates a matrix with N!/((N-M)! M!) rows and M columns containing
%   all possible combinations of N elements taken M at a time.
%
%   This function is only practical for situations where M is less
%   than about 15.

%   B.A. Jones 2-17-95
%   Copyright (c) 1984-97 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1999/08/06 09:55:17 $
%   also for matlab V4

if nargin~=2, error('Requires 2 input arguments.'); end

v = v(:).'; % Make sure v is a row vector.
n = length(v);
if n == m
   P = v;
elseif m == 1
   P = v.';
else
   P = [];
   if m < n & m > 1
      for k = 1:n-m+1
         Q = combs_(v(k+1:n),m-1);
         P = [P; [v(ones(size(Q,1),1),k) Q]];
      end
   end
end
