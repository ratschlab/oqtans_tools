function x = equal(M1, M2)
% EQUAL(M1, M2) tests M1 and M2 for equality
% we need a special function as the two matrices may be of different
% length and '==' might exit telling that the matrix dimensions don't agree
%
% returns 1 if true and 0 if false

% Copyright (c) 1997  GMD Berlin
% --- All Rights Reserved ---
% THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF 
% GMD Berlin and Lucent Technologies
% The copyright notice above does not evidence any
% actual or intended publication of this work.
%
% Authors: Alex J. Smola, Gunnar R"atsch
% Date   : 12/08/97
%

if all(size(M1)==size(M2)),
	x = all(all(M1==M2)) ;
else
  if isempty(M1)&isempty(M2),
    x=1 ;
  else
	x = 0;
  end ;
end