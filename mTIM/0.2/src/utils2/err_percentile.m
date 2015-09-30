function [m, e]=err_percentile(A, percentile);
% [m,e]=err_percentile(A, percentile);

assert(length(percentile)==1) ;

sA=sort(A) ;
med=median(sA) ;
if size(sA,1)==1,
  sA=sA' ;
end ;
idx=percentile*(size(sA,1)+1) ;
if idx<1,
  idx=1 ;
end ;
if idx>size(sA,1)
  idx=size(sA,1) ;
  warning('wrong parameter value: percentile') ;
end ;

if floor(idx)~=idx
  idx=[floor(idx) ; ceil(idx)] ;
  m=mean(sA(idx,:),1) ;
else
  m=sA(idx,:) ;
end ;

e=abs(med-m) ;