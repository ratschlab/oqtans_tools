function [q2,q2e]=q2err(O, L) 
% [q2,q2e]=q2err(O, L) 
% 
% Computes Kristin Bennet's error measure
% 
% O    is the output of the predictor
% L    is the label, i.e. desired output

mp = L-mean(L) ;
q2 = sum( (L-O).^2 ) / sum(mp.^2) ;
q2e= 1 - sum( (O-mean(O)) .* mp ) / sum(mp.^2) ;
