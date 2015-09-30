function [P_VALUE, RET_STRUCT]= rDiff_poisson(CFG,gene,Counts_rDiff_parametric,Gene_expression)

% Calculates the p-Values of a poisson test on each
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
  TEMP_SAMPLE1(j)=and(not(isempty(Counts_rDiff_parametric{j})), ...
		      TEMP_SAMPLE1(j));
end
for j=1:length(TEMP_SAMPLE2)
  TEMP_SAMPLE2(j)=and(not(isempty(Counts_rDiff_parametric{j})), ...
		      TEMP_SAMPLE2(j));
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
gene_expression_1=sum(Gene_expression(SAMPLE1));
gene_expression_2=sum(Gene_expression(SAMPLE2));

%compute the p-values
OBSERVED_COUNTS=[sum(region_counts_1,1);sum(region_counts_2,1)];
TOTAL_COUNTS=sum(OBSERVED_COUNTS,1);

%calculate gene expression ratio
LAMBDA=gene_expression_1/(gene_expression_2+gene_expression_1);


P_LIST=ones(1,size(TOTAL_COUNTS,2));
%Iterate over the regions
SKIPPED_TESTS=0;
for i=1:length(P_LIST)
    if sum(TOTAL_COUNTS(i))==0
        SKIPPED_TESTS=SKIPPED_TESTS+1;
        continue
    end
    
    MEAN=LAMBDA*TOTAL_COUNTS(i);
    VARIANCE=(TOTAL_COUNTS(i)*LAMBDA*(1-LAMBDA)).^0.5;      
    Y=sum(((MEAN-OBSERVED_COUNTS(1,i))./VARIANCE).^2).^0.5;
    
    %Calculate the p-value
    P_VALUE=1-gammainc((Y^2)/2,1/2);
    
    P_LIST(i)=P_VALUE;
end
if length(P_LIST)-SKIPPED_TESTS<=0
  P_VALUE=1;
else
  P_VALUE=min(P_LIST)*(length(P_LIST)-SKIPPED_TESTS);
end

RET_STRUCT={};
return



