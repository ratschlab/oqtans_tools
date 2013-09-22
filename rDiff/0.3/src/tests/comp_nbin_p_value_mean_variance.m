function [P_VALUE,FLAG_ERR]=comp_nbin_p_value_mean_variance(MEAN1,MEAN2,SIGMA1,SIGMA2,COUNTS1,COUNTS2)
%This function computes the p_value for a negative binomial
%hypothesis test given the mean and variance of two nbin
%distributions and counts from the respective distributions

  FLAG_ERR=0;
  r1=(MEAN1^2)/(SIGMA1-MEAN1);
  p1=1-(MEAN1/(SIGMA1));
  r2=(MEAN2^2)/(SIGMA2-MEAN2);
  p2=1-(MEAN2/(SIGMA2));
  
  P_OBS=calcl_nbin_pdf(COUNTS1,r1,1-p1)*calcl_nbin_pdf(COUNTS2,r2,1-p2);
  
  
  P_A=calcl_nbin_pdf(0:(COUNTS1+COUNTS2),r1,1-p1);
  P_B=calcl_nbin_pdf((COUNTS1+COUNTS2):(-1):0,r2,1-p2);
  P_COMB=sort(P_A.*P_B);

  if sum(P_COMB)>0
    POS= find(P_COMB<=P_OBS,1,'last');
    if isempty(POS)
      P_VALUE=min(P_COMB);
    else
      P_VALUE=sum(P_COMB(1:POS))/sum(P_COMB);
    end
  else
    P_VALUE=1;
  end
  
  
  