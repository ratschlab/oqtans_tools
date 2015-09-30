function [d,D]=relentropy(x,y)

idx=find( (x>0) ) ;
d=sum(x(idx).*log(x(idx)./y(idx))) ;
if nargout>1,
  D(idx)=x(idx).*log(x(idx)./y(idx)) ;
end ;


