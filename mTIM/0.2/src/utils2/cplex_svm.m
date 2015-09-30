function [w,b,xis,alphas] = cplex_svm(X, Y, Cp, Cn) ;
% [w,b,xis,alphas] = cplex_svm(X, Y, Cp, Cn) ;

%X=randn(20,100000) ;
%Y=sign(X(1,:)) ;
%C= 1 ;

if nargin<4,
  Cn=Cp ;
end ;

Cs = ones(length(Y),1) ;
Cs(Y==1) = Cp ;
Cs(Y==-1) = Cn ;

lpenv=cplex_license(0) ;

d=size(X,1) ;
N=size(X,2) ;
INF=1e20 ;

% v=[w b xis]
lb = [-INF*ones(d+1,1); zeros(N,1)] ;
ub = INF*ones(d+1+N,1) ;
f = [zeros(d+1,1); Cs] ;
Q=spdiag([ones(d,1); zeros(N+1,1)]) ;

% y_i*X'*w + y_i*b >= 1 - xi_i
% -y_i*X'*w - y_i*b - xi_i <= -1
A = [-spdiag(Y)*X' -Y' -speye(N)] ;
b = -ones(N,1) ;

tic;
  [res,alphas,how]=qp_solve(lpenv, Q, f, A, b, lb, ub, 0, 1, 'bar') ; % 'sift'
toc;
assert(isequal(how,'OK')) ;

w=res(1:d) ;
b=res(d+1) ;
xis=res(d+2:end) ;

