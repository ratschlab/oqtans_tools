function E=kldiv(p, q) 
% E=kldiv(Out, DestOut) 
% 

%#realonly
%#inbounds

E=sum(p.*log(q./p)) ;
