function [sig,t]=t_test(ma, std_a, la, mb, std_b, lb) ;
%


sd=sqrt( ( std_a^2*(la-1) + std_b^2*(lb-1) )*(1/la+1/lb)/(la+lb-2)) ;
t=abs((ma-mb)/sd) ;

nu=la+lb-2 ;
sig=betainc(nu/(nu+t^2), nu/2, 1/2) ;