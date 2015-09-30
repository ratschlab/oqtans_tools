function b = subsref(a,s)
%    B = SUBSREF(A,S) is called for the syntax A(I), A{I}, or A.I
%    when A is an object.  S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts.
%

% File:        @cellfile/subsref.m
%
% Author:      Gunnar R"atsch
% Created:     04/20/04
% 
% Copyright (c) 2004  Fraunhofer Berlin & MPI for biol. Cybernetics
% Tuebingen - All rights reserved
% THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE of Fraunhofer/MPI
% The copyright notice above does not evidence any
% actual or intended publication of this work.

switch s(1).type
 case '{}'
  tmp = s(1).subs;
  if length(tmp)~=1
    error('bad operation') ;
  end ;
  [a,b]=getelem(a,tmp{1}) ;
  if length(s)>1
    b=subsref(b,s(2:end)) ;
  end ;
 case '()'
  if length(s)~=1 | length(s(1).subs)~=1,
    error('bad operation') ;
  end ;
  tmp = s(1).subs{1};
  b=cell(1,length(tmp)) ;
  for i=1:length(tmp)
    [a,b{i}]=getelem(a,tmp(i)) ;
  end 
 otherwise
  error('inapropriate operation for cellfile object') ;
end ;

