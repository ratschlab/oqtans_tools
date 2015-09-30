function cf=setelem(cf, id, elem)
% cf=setelem(cf, id, elem)

if length(id)~=1,
  error('cannot set multiple elements') ;
end;
fname=getelem_filename(cf, id) ;

save(fname, 'elem') ;

cf.cache_id=id ;
cf.cache_elem=elem ;

cf.max_elem=max(cf.max_elem,id) ;
