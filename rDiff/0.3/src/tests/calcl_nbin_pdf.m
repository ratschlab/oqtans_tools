function [RET]=calcl_nbin_pdf(VAL,R,P)
%uses formula used by matlab for compatability

RET=zeros(size(VAL,1),size(VAL,2));
IDX=VAL==round(VAL);
%RET(IDX) = exp(gammaln(R+VAL(IDX))-gammaln(VAL(IDX)+1)-gammaln(R)+R*log(P)+VAL(IDX)*log1p(-P));


nk = (gammaln(R + VAL(IDX)) - gammaln(VAL(IDX) + 1) - gammaln(R));
RET(IDX) = nk + (R.*log(P)) + (VAL(IDX).*log(1 - P));
RET(IDX) = exp(RET(IDX));
