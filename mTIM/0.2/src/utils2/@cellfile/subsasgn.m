function a = subsasgn(a,s,b)
%    A = SUBSASGN(A,S,B) is called for the syntax A(I)=B, A{I}=B, or
%
%    A.I=B when A is an object.  S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts. 

% File:        @cellfile/subasgn.m
%
% Author:      Gunnar R"atsch
% Created:     04/20/04
% 
% Copyright (c) 2004  Fraunhofer Berlin & MPI for biol. Cybernetics
% Tuebingen - All rights reserved
% THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE of Fraunhofer/MPI
% The copyright notice above does not evidence any
% actual or intended publication of this work.

switch s(1).type,
 case '{}',
  tmp=s(1).subs ;
  if length(tmp)~=1,
    error('bad operation') ;
  end ;
  if length(s)>1,
    [a,b2]=getelem(a,tmp{1}) ;
    b2=subsasgn(b2,s(2:end),b) ;
    a=setelem(a,tmp{1},b2) ;
  else
    a = setelem(a,tmp{1},b) ;
  end ;
 case '()',
  tmp = s(1).subs;
  if length(s)~=1 | length(tmp)~=1 
    error('bad operation') ;
  end ;
  for j=1:length(tmp{1}),
    a=setelem(a,tmp{1}(j),b{j}) ;
  end
 case '.',
  error('inapropriate operation for cellfile object') ;
 otherwise,
  error('???') ;
end ;
