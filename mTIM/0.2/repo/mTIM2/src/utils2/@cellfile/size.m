function [a,b]=size(cf,dim)
% l=length(cf)

if nargout<=1
  a=[1 cf.max_elem] ;
  if nargin==2,
    a=a(dim) ;
  end ;
else
  a=1 ;
  b=cf.max_elem ;
end ;
