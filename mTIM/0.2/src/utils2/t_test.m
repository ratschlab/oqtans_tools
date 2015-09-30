function [sig,t]=t_test(a,b) ;
%

ma=mean(a) ;
la=length(a) ;
mb=mean(b) ;
lb=length(b) ;

sd=sqrt( ( sum( (a-ma).^2 ) + sum( (b-mb).^2 ) )*(1/la+1/lb)/(la+lb-2)) ;
t=abs((ma-mb)/sd) ;

nu=la+lb-2 ;
sig=betainc(nu/(nu+t^2), nu/2, 1/2) ;