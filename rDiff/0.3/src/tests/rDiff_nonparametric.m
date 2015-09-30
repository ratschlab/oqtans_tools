function [pval,info] = rDiff_nonparametric(CFG,READS1,READS2,variance_function_nonparametric_1, variance_function_nonparametric_2)

bootstraps=CFG.bootstraps;

%Initizialize the reads
SUMS1=[]; %the readcoverage for sample 1
SUMS2=[]; %the readcoverage for sample 2
reads1=[]; %the total of reads in sample 1
reads2=[]; %the total of reads in sample 2
NON_EMPTY=[];
for i=1:length(READS1)
  if not(isempty(READS1{i}))
    NON_EMPTY=[NON_EMPTY,i];
  end
end
READS1={READS1{ NON_EMPTY}};
NON_EMPTY=[];
for i=1:length(READS2)
  if not(isempty(READS2{i}))
    NON_EMPTY=[NON_EMPTY,i];
  end
end
READS2={READS2{NON_EMPTY}};
for i=1:length(READS1)
  SUMS1=[SUMS1;sum(READS1{i},1)];
  reads1=[reads1;READS1{i}];
end
for i=1:length(READS2)
  SUMS2=[SUMS2;sum(READS2{i},1)];
  reads2=[reads2;READS2{i}];
end

if size(reads1,1)==0 || size(reads2,1)==0
  pval=1;
  info=1;
  bootstrap_results=1;
  statistic=1;
  return
end


%Get masks for the highly covered positions
[MASKS]=get_nonparametric_masks(CFG,reads1,reads2);


%Determine number of masks
NR_OF_MASKS=size(MASKS,1);

% Predict variance
VARIANCE1=(1e-8)*ones(size(reads1,2),1);
for i=1:length(READS1)
    TEMP_VARIANCE1=predict_variance(sum(READS1{i},1)',variance_function_nonparametric_1);
    TEMP_VARIANCE1(isnan(TEMP_VARIANCE1))=0;
    TEMP_VARIANCE1(isinf(TEMP_VARIANCE1))=0;
    VARIANCE1=VARIANCE1+TEMP_VARIANCE1;
end

VARIANCE2=(1e-8)*ones(size(reads1,2),1);
for i=1:length(READS2)
    TEMP_VARIANCE2=predict_variance(sum(READS2{i},1)',variance_function_nonparametric_2);
    TEMP_VARIANCE2(isnan(TEMP_VARIANCE2))=0;
    TEMP_VARIANCE2(isinf(TEMP_VARIANCE2))=0;
    VARIANCE2=VARIANCE2+TEMP_VARIANCE2;
end

%Get the mean variance
VARIANCE1=VARIANCE1/length(READS1);
VARIANCE2=VARIANCE2/length(READS2);

%Concatenation of all reads
allreads = [reads1;reads2];


%Determine the subsampling rate
p=sum(allreads,1)/size(allreads,1);
R=size(reads1,1);
c=(p.*(1-p))/(size(allreads,1)-1);
N1s=size(allreads,1)*c./(c+(VARIANCE1'/(R^2)));

p=sum(allreads,1)/size(allreads,1);
R=size(reads2,1);
c=(p.*(1-p))/(size(allreads,1)-1);
N2s=size(allreads,1)*c./(c+(VARIANCE2'/(R^2)));

%Round to get integer values
N1s=ceil(full(N1s));
N2s=ceil(full(N2s));


%Total number of reads in each replicate
N1 = size(reads1,1);
N2 = size(reads2,1);
N = N1 + N2;


bootstrap_results=ones(NR_OF_MASKS,bootstraps);
COUNTER=0;
R1s=[];
R2s=[];
statistic=ones(NR_OF_MASKS,bootstraps);
nr_of_slices=CFG.nr_of_slices;

TOTAL_SUM1=sum(reads1,1);
TOTAL_SUM2=sum(reads2,1);
TOTAL_SUM=TOTAL_SUM1+TOTAL_SUM2;

clear reads1
clear reads2
%Use the transpose to make the selection of columms faster
all_reads_trans=allreads';
clear allreads
if length(unique(TOTAL_SUM))<nr_of_slices+5
    pval=1;
    bootstrap_results=[];
    statistic=[];
    pval=1;
    info = 1;
    return
end

[SUM_SORTED,SUM_POS]=sort(TOTAL_SUM);
NR_OF_ZEROS=sum(TOTAL_SUM==0);

%precompute this sum once 
all_reads_trans_row_sum=sum(all_reads_trans,1);

%These field contain
STAT_DIST=zeros(NR_OF_MASKS,nr_of_slices);
TEMP_DIST=zeros(1,nr_of_slices);

%Precompute normalizing constants
SUM_TOTAL_SUM1=sum(TOTAL_SUM1);
SUM_TOTAL_SUM2=sum(TOTAL_SUM2);
SUM_SUM_TOTAL_SUM=SUM_TOTAL_SUM1+SUM_TOTAL_SUM2;

%Precompute the some of the values or the slices
SLICES=zeros(nr_of_slices,size(TOTAL_SUM,2))==1;
N1_arr=zeros(nr_of_slices,1);
N2_arr=zeros(nr_of_slices,1);
FACT_arr=zeros(nr_of_slices,1);
V1_cell=cell(nr_of_slices,1);
V2_cell=cell(nr_of_slices,1);

%Detemine regions wher to match the variances
for s=nr_of_slices:(-1):1
        LOWER_POS=min(NR_OF_ZEROS+ceil(((nr_of_slices-s)/nr_of_slices)*(length(TOTAL_SUM)-NR_OF_ZEROS))+1,length(TOTAL_SUM));
        UPPER_POS=min(NR_OF_ZEROS+ceil(((nr_of_slices-s+1)/nr_of_slices)*(length(TOTAL_SUM)-NR_OF_ZEROS))+1,length(TOTAL_SUM));
        
        SLICE_LOWER_BOUND=SUM_SORTED(LOWER_POS);
        SLICE_UPPER_BOUND=SUM_SORTED(UPPER_POS);
        if s==nr_of_slices
	  SLICE=and(TOTAL_SUM>=SLICE_LOWER_BOUND,TOTAL_SUM<=SLICE_UPPER_BOUND);
	else
	  SLICE=and(TOTAL_SUM>SLICE_LOWER_BOUND,TOTAL_SUM<=SLICE_UPPER_BOUND);
        end
	SLICES(s,:)=SLICE;
        
        N1s_temp=ceil(median(N1s(SLICE)));
        N2s_temp=ceil(median(N2s(SLICE)));
        N1s_temp=min(N1s_temp,N1);
        N2s_temp=min(N2s_temp,N2);
        
        N1_arr(s)=N1s_temp;
        N2_arr(s)=N2s_temp;
        
        FACT_arr(s)=sum(TOTAL_SUM1(SLICE)+TOTAL_SUM2(SLICE))/(SUM_SUM_TOTAL_SUM);
        
        V1_cell{s}=TOTAL_SUM1(SLICE)/SUM_TOTAL_SUM1;%temporary variable to safe time
        V2_cell{s}=TOTAL_SUM2(SLICE)/SUM_TOTAL_SUM2;%temporary variable to safe time
        for mask_ix=1:NR_OF_MASKS
            STAT_DIST(mask_ix,s)=eucl_dist_weigthed(V1_cell{s},V2_cell{s},MASKS(mask_ix,SLICES(s,:)));
        end
end
SUM_SLICES=sum(SLICES,2);
STAT_DIST(isnan(TEMP_DIST))=0;
all_reads_trans_slices=cell(nr_of_slices,1);
for s=nr_of_slices:(-1):1
    all_reads_trans_slices{s}=all_reads_trans(SLICES(s,:),:);
end

for i = 1:bootstraps
    % permutation of the reads
    read_per = randperm(N);
     
    %Peform the computation for each region where the variances are matched
    for s=nr_of_slices:(-1):1
        if SUM_SLICES(s)==0;
            continue;
        end 
	%Create random samples 1 and 2
        sample1 = sum(all_reads_trans_slices{s}(:,read_per(1:N1_arr(s))),2);
        sample2 = sum(all_reads_trans_slices{s}(:,read_per((N1_arr(s)+1):(N1_arr(s)+N2_arr(s)))),2);
        
        W1=sample1/sum(sample1)*FACT_arr(s);%temporary variable to safe time
        W2=sample2/sum(sample2)*FACT_arr(s);%temporary variable to safe time
       
        for mask_ix=1:NR_OF_MASKS
            TEMP_DIST(mask_ix,s)=eucl_dist_weigthed(W1,W2,MASKS(mask_ix,SLICES(s,:))');
        end
        
    end 
    
    %make sure the normalisation doe not intruduces nan's
    TEMP_DIST(isnan(TEMP_DIST))=0;
   
    COUNTER=COUNTER+1;  
    %Compute the average from the different matching regionon
    statistic(:,COUNTER)=mean(STAT_DIST,2);
    bootstrap_results(:,COUNTER)=mean(TEMP_DIST,2);
end

bootstrap_results=bootstrap_results(:,1:COUNTER);
statistic=statistic(:,1:COUNTER);

pval=double(sum(bootstrap_results>=statistic,2)) / COUNTER;

info = {bootstrap_results,statistic,pval};
pval=min(pval)*10;

function result = eucl_dist_weigthed(A,B,W)

result = sqrt(sum( W.*((A - B) .^ 2) ));
