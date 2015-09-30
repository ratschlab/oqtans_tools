function [VARIANCE]=estimate_variance(MEANS,VARS) 
%estimates the variance using the locfit package

IDX_ZERO=and(VARS>=0,MEANS>0);
VARIANCE =locfit(log(MEANS(IDX_ZERO)),VARS(IDX_ZERO),'family','gamma','acri','cp','deg',2,'maxk',1000);

return
