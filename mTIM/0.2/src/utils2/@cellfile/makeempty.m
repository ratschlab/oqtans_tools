function cf = makeempty(cf) ;

unix(sprintf('find %s -name "[0-9]*.mat" -print -exec rm {} \\;', cf.fdir)) ;

cf.max_elem=0 ;
cf.cache_id=[] ;
cf.cache_elem=[] ;