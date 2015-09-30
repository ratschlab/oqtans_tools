function [VARIANCE]=predict_variance(MEANS,VARIANCE_FCT)
%predict the variances for the means
if isstruct(VARIANCE_FCT)
    %use the estimated variance function
    VARIANCE=exp(predict(VARIANCE_FCT,log(full(MEANS))));
else
    %use parmeters a,b,c to compute the variance function
    a=VARIANCE_FCT(1);
    b=VARIANCE_FCT(2);
    c=VARIANCE_FCT(3);
    VARIANCE=a+MEANS*b+(MEANS.^2)*c;
end

%Make sure the variance is bigger than the poisson
%variance. Otherwise the NB cannot be definied

VARIANCE=max([VARIANCE,MEANS*(1e-8)],[],2);

return
