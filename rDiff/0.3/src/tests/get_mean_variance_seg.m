function [COMB_MEAN,COMB_VARIANCE]=get_mean_variance_seg(gene_expression_1,gene_expression_2,region_counts_1,region_counts_2,variance_function_parametric_1, variance_function_parametric_2)


COMB_READS_PER_EXON=[region_counts_1;region_counts_2];
GENE_EXPRESSION=[gene_expression_1';gene_expression_2'];
IX_SAMPLE1=1:length(gene_expression_1);
IX_SAMPLE2=(1+length(gene_expression_1)):(length(gene_expression_1)+length(gene_expression_2));
  
COMB_VARIANCE=[];
COMB_MEAN=[];
for i=1:size(region_counts_1,2)
    INTEN1=region_counts_1(:,i);
    INTEN2=region_counts_2(:,i);
    INTENSITY=[INTEN1;INTEN2];
    
    SR_VECT=GENE_EXPRESSION;
    
    %Are there any counts at all in the region?
    if sum(INTENSITY)>0
        %Compute the means under the null hypothesis
        Q=(INTENSITY./SR_VECT)/sum(SR_VECT>0);
        
        if sum(isnan(SR_VECT))
            MEAN1=0;
            MEAN2=0;
        else
            MEAN1=mean(sum(Q)*SR_VECT(IX_SAMPLE1));
            MEAN2=mean(sum(Q)*SR_VECT(IX_SAMPLE2));
        end
        
        COMB_MEAN=[COMB_MEAN,[MEAN1;MEAN2]];
    else
        COMB_MEAN=[COMB_MEAN,[0;0]];
    end
    
end

VARIANCE1= predict_variance(COMB_MEAN(1,:)',variance_function_parametric_1)';
VARIANCE2= predict_variance(COMB_MEAN(2,:)',variance_function_parametric_2)';

COMB_VARIANCE=[VARIANCE1;VARIANCE2];