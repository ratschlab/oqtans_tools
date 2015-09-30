function f = separate(str, delim, merge_multiple)
% SEPARATE  Separates input string by given delimiters.
%
%   f = separate(str, delimiter, merge_multiple_delimiters)
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


if nargin<2, delim=sprintf('\t') ; end ;
if nargin<3, merge_multiple=0 ; end ;
f={} ;
if isempty(delim)
  error('empty delimiter')
end
if isempty(str)
  return
end
idx=[0 find(str==delim) length(str)+1] ;
if merge_multiple,
  idx2=find(idx(1:end-1)~=idx(2:end)-1) ;
  idx=[idx(idx2) idx(end)] ;
end 
for i=1:length(idx)-1
  f{i}=deblank(str(idx(i)+1:idx(i+1)-1)) ;
end ;
