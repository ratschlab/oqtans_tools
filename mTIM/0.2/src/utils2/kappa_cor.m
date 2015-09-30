function kappa=kappa_cor(x1,x2)
%
C(1,1)=sum(x1==1 & x2==1) ;
C(1,2)=sum(x1==1 & x2==-1) ;
C(2,1)=sum(x1==-1 & x2==1) ;
C(2,2)=sum(x1==-1 & x2==-1) ;

l=length(x1) ;
theta1=sum(diag(C))/l ;

theta2=(C(1,1)+C(1,2))*(C(1,1)+C(2,1))/(l^2)+ ...
       (C(2,1)+C(2,2))*(C(1,2)+C(2,2))/(l^2) ;

kappa=(theta1-theta2)/(1-theta2) ;