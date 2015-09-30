function MSE=mse(Out, DestOut) 
% MSE=mse(Out, DestOut) 
% 
% Computes the mean-squared-error

%#realonly
%#inbounds

[dim,p]=size(DestOut) ;
MSE=sum(sum((Out-DestOut).^2))/p ;