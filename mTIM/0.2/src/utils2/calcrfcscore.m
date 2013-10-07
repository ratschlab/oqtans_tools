function [score,tp,tdr]=calcrfcscore(output,LT)
% score=calcrfcscore(output,LT)
map=(LT==1) ;
[s,idx]=sort(-output) ; s=-s ;
map=map(idx) ;

tp=cumsum(map)/sum(LT==1) ;
tdr=cumsum(map)./[1:length(map)] ;

t=tp(2:end)-tp(1:end-1) ;
score = sum(0.5*[tdr(1:end-1)+tdr(2:end)].*t) ;


