function A=spdiag(diagonal) ;
% A=spdiag(diagonal) ;

if size(diagonal,1)>1
  A=spdiags(diagonal, 0, length(diagonal),length(diagonal)) ;
else
  A=spdiags(diagonal', 0, length(diagonal),length(diagonal)) ;
end ;
