function create_subdirs(cf)

if ~exist(cf.fdir,'dir')
  unix(['mkdir ' cf.fdir]) ;
end ;

cf.depth=cf.depth-1 ;
if cf.depth<0, return ; end ;
cf2=cf ;
for i=0:cf.width-1,
  cf2.fdir=sprintf('%s%i/', cf.fdir, i) ;
  if cf.depth>0
    fprintf('creating %s\n', cf2.fdir) ;
  end ;
  create_subdirs(cf2) ;
end ;
