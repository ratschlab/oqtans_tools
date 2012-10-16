function [cf,elem]=getelem(cf, id)
% [cf,elem]=getelem(cf, id)

if isempty(id), elem=[] ; return ; end ;
if length(id)~=1,
  error('cannot get multiple elements') ;
end;

if equal(cf.cache_id, id),
  elem=cf.cache_elem ;
  return ;
end ;

fname=getelem_filename(cf, id) ;
if id>cf.max_elem & fexist(fname)
  cf.max_elem=id ;
end ;

if id>cf.max_elem | id<1
  error('index exceeds cellfile dimension') ;
end ;


if fexist(fname),
  L=load(fname, 'elem') ;
  cf.cache_id=id ;
  cf.cache_elem=L.elem ;
  elem=L.elem ;
else
  elem=[] ;
end ;

