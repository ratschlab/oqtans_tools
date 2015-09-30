function [x_tls, E, r]= tls(A, b, d, t)
% function [x_tls, E, r]= tls(A, b, d, t)
% 
% This functions computes the Total least sqares solution
% to the system Ax=b.
%
% x_tls=arg min || diag(d) [E r] diag(t) ||_F     s.t. (A+E)x=b+r 
%
% Remark: This algorithm requires about 2mn^2+12n^3 flops

[m,n]=size(A) ;

if nargin<3,
	d=ones(m,1) ;
end ;
if nargin<4,
	t=[1*ones(n,1) ; 1e6];
end ;

if abs(prod(d))<10*eps,
	deb_output(sprintf('product of the diagonal entries of D: %e', prod(d))) ;
	error('Matrix D=diag(d) ist close to singular') ;
end ;
if abs(prod(t))<10*eps,
	deb_output(sprintf('product of the diagonal entries of T: %e', prod(t))) ;
	error('Matrix T=diag(t) ist close to singular') ;
end ;

% comute the svd 
[U,S,V]=svd(diag(d)*[A b]*diag(t),0) ;
clear U ;

diag(log(S))

% determine p
s=diag(S) ;
p=0 ; 
while ( s(n-p) == s(n+1) ),
	p=p+1 
end ;

% compute a householder matrix P such that if VT=VP,
if p>0,
	[h,beta]=house( V(n+1,1:n+1)' ) ;
	VT=V-(beta*(V*h))*h' ;
else
	VT=V ;
end ;
% then must be V(n+1,n-p+1:n)=0

if VT(n+1,n+1)~=0,
	for i=1:n,
		x_tls(i)= - t(i) * VT(i,n+1) / ( t(n+1) * VT(n+1,n+1) ) ;
	end ;
end ;
x_tls


