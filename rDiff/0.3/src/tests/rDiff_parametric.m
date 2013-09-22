function [P_VALUE, RET_STRUCT]= rDiff_parametric(CFG,gene,Counts_rDiff_parametric,Gene_expression,variance_function_parametric_1, variance_function_parametric_2)
 
% Calculates the p-Values of a negative binomial test on each
% alternative regions and combines the p-values using Bonferroni's correction


%Initialize gene.name
NR_OF_TRANS=size(gene.transcripts,2);
if NR_OF_TRANS<=1
  RET_STRUCT='NR_OF_TRANS too small';
  P_VALUE=1;
  return
end


%get the samples that are expressed (have more than 10 reads)
TEMP_SAMPLE1=and(CFG.SAMPLES==1,Gene_expression>=10);
TEMP_SAMPLE2=and(CFG.SAMPLES==2,Gene_expression>=10);
SAMPLE1=find(TEMP_SAMPLE1);
SAMPLE2=find(TEMP_SAMPLE2);

%Check wether Counts_rDiff_parametric is nonempty
for j=1:length(TEMP_SAMPLE1)
  TEMP_SAMPLE1(j)=and(not(isempty(Counts_rDiff_parametric{j})),TEMP_SAMPLE1(j));
end
for j=1:length(TEMP_SAMPLE2)
  TEMP_SAMPLE2(j)=and(not(isempty(Counts_rDiff_parametric{j})),TEMP_SAMPLE2(j));
end

SAMPLE1=find(TEMP_SAMPLE1);
SAMPLE2=find(TEMP_SAMPLE2);


SAMPLE_LENGTH1=length(SAMPLE1);
SAMPLE_LENGTH2=length(SAMPLE2);

if min(SAMPLE_LENGTH1,SAMPLE_LENGTH2)==0
    RET_STRUCT='SAMPLE_LENGTH too small';
    P_VALUE=1;
    return
end
 
% Get the region counts
region_counts_1=zeros(SAMPLE_LENGTH1,length(Counts_rDiff_parametric{1,1}));
for j=1:SAMPLE_LENGTH1
    region_counts_1(j,:)=Counts_rDiff_parametric{1,SAMPLE1(j)};
end
region_counts_2=zeros(SAMPLE_LENGTH2,length(Counts_rDiff_parametric{1,1}));
for j=1:SAMPLE_LENGTH2
    region_counts_2(j,:)=Counts_rDiff_parametric{1,SAMPLE2(j)};
end

% Get the gene expression
gene_expression_1=Gene_expression(SAMPLE1);
gene_expression_2=Gene_expression(SAMPLE2);

% compute the expected mean and the variance under the null
% hypothesis
[EXPECTED_MEAN,EXPECTED_VARIANCE]=get_mean_variance_seg(gene_expression_1,gene_expression_2,region_counts_1,region_counts_2,variance_function_parametric_1, variance_function_parametric_2);

%compute the p-values
OBSERVED_COUNTS=round([mean(region_counts_1,1);mean(region_counts_2,1)]);

P_LIST=ones(1,size(EXPECTED_MEAN,2));
%Iterate over the regions
SKIPPED_TESTS=0;
for i=1:length(P_LIST)
    if sum(OBSERVED_COUNTS(:,i))==0
        SKIPPED_TESTS=SKIPPED_TESTS+1;
        continue
    end
    
    [P_VALUE,FL]=comp_nbin_p_value_mean_variance(EXPECTED_MEAN(1,i),EXPECTED_MEAN(2,i),EXPECTED_VARIANCE(1,i),EXPECTED_VARIANCE(2,i),OBSERVED_COUNTS(1,i),OBSERVED_COUNTS(2,i));
    P_LIST(i)=P_VALUE;
end
if length(P_LIST)-SKIPPED_TESTS<=0
  P_VALUE=1;
else
  P_VALUE=min(P_LIST)*(length(P_LIST)-SKIPPED_TESTS);
end
RET_STRUCT={};
return