function A=spzeros(m,n,szmax) ;
% A=spzeros(m,n) ;

if nargin<=1
  A=sparse([],[],[],m,m) ;
  return ;
end 
if nargin<=2
  A=sparse([],[],[],m,n) ;
  return ;
end 
A=sparse([],[],[],m,n,szmax) ;
