function fname=getelem_filename(cf, idx)
% fname=getelem_filename(cf, idx)

i=floor((idx-1)/cf.width) ; fname=cf.fdir ;
for k=1:cf.depth,
  fname=sprintf('%s%i/', fname, mod(i,cf.width)) ;
  i=floor(i/cf.width) ;
end ;
fname=sprintf('%s%i.mat', fname, idx-1) ;
