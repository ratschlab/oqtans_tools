function [VARIANCE]=estimate_variance_helper(CFG,SAMPLE_IX,COUNTS_PER_GENE,GENE_EXPRESSION_MATRIX)


%compute means and variances
%Get the number of of positions to be used in MEANS and VARS
LENGTH=0;
for i=1:size(COUNTS_PER_GENE,1)
    try
        TEMP_SUM=COUNTS_PER_GENE{i, SAMPLE_IX(1)};
        for j=2:length(SAMPLE_IX)
            TEMP_SUM=TEMP_SUM+COUNTS_PER_GENE{i,SAMPLE_IX(j)};
        end
        LENGTH=LENGTH+sum(TEMP_SUM>0);
    catch
        LENGTH=LENGTH;
    end
end

%Compute the means and variances, normalizing for gene expression
MEANS=zeros(1,LENGTH);
VARS=zeros(1,LENGTH);
COUNTER=1;
for i=1:size(COUNTS_PER_GENE,1)
    try
        TS=GENE_EXPRESSION_MATRIX(i, SAMPLE_IX);
        S=length(TS)*TS/sum(TS);
        
        TEMP_SUM=sparse(zeros(length(SAMPLE_IX),size(COUNTS_PER_GENE{i, SAMPLE_IX(1)},2)));
	%if S(1)==0
	%  TEMP_SUM(1,:)=0;
	%else
	%  TEMP_SUM(1,:)=COUNTS_PER_GENE{i, SAMPLE_IX(1)}/S(1);
        %end
	for j=1:length(SAMPLE_IX)
	  if S(j)==0
	    TEMP_SUM(j,:)=0;
	  else
            TEMP_SUM(j,:)=COUNTS_PER_GENE{i,SAMPLE_IX(j)}/S(j);
	  end
	end
        TEMP_SUM=TEMP_SUM(:,sum(TEMP_SUM,1)>0);
        MEANS(COUNTER:(COUNTER+size(TEMP_SUM,2)-1))=mean(TEMP_SUM(:,sum(TEMP_SUM,1)>0),1);
        VARS(COUNTER:(COUNTER+size(TEMP_SUM,2)-1))=var(TEMP_SUM(:,sum(TEMP_SUM,1)>0));
        COUNTER=COUNTER+size(TEMP_SUM,2);
    catch
        LENGTH=LENGTH;
    end
end

%filter those which have a zeros mean
NONZERO_IDX=MEANS>0;
MEANS=MEANS(NONZERO_IDX);
VARS=VARS(NONZERO_IDX);
% Subsample the mean variance pairs to reduce the number of
% samples which have a low mean
[MEANS,POS]=sort(MEANS);
MAX_VAL=max(MEANS);
MIN_VAL=min(MEANS);
BOUNDS=exp(linspace(log(MIN_VAL),log(MAX_VAL), CFG.variance_samplebins+1));
SAMPLE=[];
for i=1:CFG.variance_samplebins
    NR_IN_BIN=length((find(MEANS>=BOUNDS(i),1,'first'):find(MEANS<BOUNDS(i+1),1,'last')));
    if NR_IN_BIN==0
        continue;
    elseif NR_IN_BIN<= CFG.variance_samples_per_bin
        SAMPLE=[SAMPLE find(MEANS>=BOUNDS(i),1,'first'):find(MEANS<BOUNDS(i+1),1,'last')];
    else
        CHUNK_SAMPLE=find(MEANS>=BOUNDS(i),1,'first'):find(MEANS<BOUNDS(i+1),1,'last');
        TEMP_SAMPLE=randperm(NR_IN_BIN);
        SAMPLE=[SAMPLE CHUNK_SAMPLE(TEMP_SAMPLE(1:CFG.variance_samples_per_bin))];
    end
end
VARS=VARS(POS);

%estimate variance function
[VARIANCE]=estimate_variance(full(MEANS(SAMPLE))',full(VARS(SAMPLE)'));
return
