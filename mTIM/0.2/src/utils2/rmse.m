function RMSE=rmse(O, L) 
% MSE=rmse(O, L) 
% 
% Computes the root of the mean-squared-error
% 
% O    is the output of the predictor
% L    is the label, i.e. desired output

[dim,p]=size(L) ;
if ~isempty(O) | ~isempty(L),
  RMSE=sqrt(sum(sum((O-L).^2))/p) ;
else
  RMSE=nan ;
end ;

