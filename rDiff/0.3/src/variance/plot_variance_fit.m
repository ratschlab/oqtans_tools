function []=plot_variance_fit(VARIANCE,MEANS,VARS,MAX_X,MAX_Y)

if nargin<=3
    MAX_X=1000;
    MAX_Y=1000;
end

%values for which to plot the variance function
pp=0.1:1:MAX_X;

%Define the legend entries
LEG={};
FIG_HANDLES=[];
figure;
if and(not(isempty(MEANS)),not(isempty(VARS)))
    LEG{end+1}='Observation';
    FIG_HANDLES(end+1)=plot(MEANS,VARS,'.');
end
hold on
%plot the predicted variances
FIG_HANDLES(end+1)=plot(pp',predict_variance(pp',VARIANCE),'r');
LEG{end+1}='Variance fit';

plot(pp,pp,'g')

%Plot sliding window
RR=zeros(1,MAX_X);
for i=1:MAX_X
    TIDX=and(MEANS>i,MEANS<i*1.1+1);
    RR(i)= mean(VARS(TIDX));
end 
FIG_HANDLES(end+1)=plot(1:MAX_X,RR,'k');
LEG{end+1}='Sliding window';

%Plot poisson
FIG_HANDLES(end+1)=plot(pp,pp,'g');
LEG{end+1}='Poisson variance';

legend(FIG_HANDLES,LEG)
xlabel('mean')
ylabel('variance')

% Change window
xlim([0,MAX_X])
ylim([0,MAX_Y])

return