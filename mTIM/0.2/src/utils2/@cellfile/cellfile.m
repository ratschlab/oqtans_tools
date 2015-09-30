function cf=cellfile(fdir,depth,width) ;
% cf=cellfile(fdir,depth) ;

cf.cache_id=[] ;
cf.cache_elem=[] ;
cf.fdir=[] ;
cf.depth=[] ;
cf.width=[] ;
cf.max_elem=[] ;

if nargin==0,
  cf=class(cf, 'cellfile') ;
  return ;
end ;

if fdir(end)~='/',
  fdir=[fdir '/'] ;
end ;
if fdir(1)=='/' | fdir(1)=='~'
  cf.fdir=fdir ;
else
  cf.fdir=[pwd '/' fdir] ;
end ;

cf.cache_id=[] ;
cf.cache_elem=[] ;

to_be_copied=0 ;
if isa(depth, 'cellfile')
  to_be_copied=1 ;
  cf.depth=depth.depth ;
  cf.width=depth.width ;
  cf.max_elem=depth.max_elem ;
  cf=class(cf, 'cellfile') ;
else
  if nargin>1,
    cf.depth=depth ;
  else
    cf.depth=1 ;  
  end ;
  if nargin>2,
    cf.width=width ;
  else
    cf.width=100 ;
  end ;
  cf.max_elem=0 ;
  
  cf=class(cf, 'cellfile') ;
end ;
create_subdirs(cf) ;


i=1 ; last_i=1 ;
while fexist(getelem_filename(cf, i))
  last_i=i ;
  i=i+1000 ; 
end ;

i=last_i ;
while fexist(getelem_filename(cf, i))
  i=i+1 ;
end ;
cf.max_elem=i-1 ;
fprintf('cellfile: found %i elements on disk\n', cf.max_elem) ;

if to_be_copied,
  fprintf('copying %i elements\n', depth.max_elem) ;
  for i=1:depth.max_elem,
    [depth,elem]=getelem(depth,i) ;
    cf=setelem(cf,i,elem) ;
    fprintf('%i\r',i) ;
  end ;
  fprintf('\n') ;
end ;
    