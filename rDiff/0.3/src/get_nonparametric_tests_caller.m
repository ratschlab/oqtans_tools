function []=get_nonparametric_tests_caller(PAR)


CFG = PAR.CFG;
genes = PAR.genes;
OUT_STR='';

% add paths
addpath(CFG.paths);
%load local variables
data_dir=CFG.data_dir;
OUT_STR=[];

variance_function_nonparametric_1=PAR.variance_function_nonparametric_1;
variance_function_nonparametric_2=PAR.variance_function_nonparametric_2;

Counts_rDiff_nonparametric=PAR.Counts_rDiff_nonparametric;
Gene_expression=PAR.Gene_expression;

%clear variabe PAR
clear PAR;

NUMBER_OF_TESTS_PER_GENE=(CFG.perform_nonparametric+CFG.perform_mmd);
if NUMBER_OF_TESTS_PER_GENE==0
    return
end
P_VALS=cell(size(genes,2),NUMBER_OF_TESTS_PER_GENE+2);


%iterate over genes
for i=1:size(genes,2)
    
    %TEMP_COUNT contains the counts for the current gene
    TEMP_COUNT=cell(1,3);
    gene = genes(i);
    
    
    OLD_OUT_STR=OUT_STR;
    OUT_STR=['Current gene: ' gene.name ];
    %print progress
    if CFG.use_rproc
        fprintf([OUT_STR '\n'])
    else
        % Erase old progress
        fprintf(repmat('\b',1,length(OLD_OUT_STR)));
        fprintf([OUT_STR])
    end
    
    %set default return values
    P_VALS{i,1}=gene.name;
    
    %check that the gene has exons defined
    if isempty(gene.exons)
        P_VALS{i,4}='Exons field empty in gene structure';
        continue;
    end
    
    %check that the gene is longer than the Reads. Otherwise the
    %definition of regions does not makes sense
    if gene.stop-gene.start<CFG.sequenced_length+3
        continue;
    end
    %perform 
    
    
    %Get the reads
    READ_SET={};
    for bam_ix=1:length(CFG.BAM_FILES)
        CFG.curr_bamfile=CFG.BAM_FILES{bam_ix};
        READ_SET{end+1} = get_reads_for_gene(CFG,gene);
    end
    %merge the reads per sample
    SAMPLE1=find(CFG.SAMPLES==1);
    SAMPLE2=find(CFG.SAMPLES==2);
    
       
    COUNTER=2;
    
    if CFG.perform_mmd
        reads1=[];
        reads2=[];
        SAMPLE_LENGTH1=length(SAMPLE1);
        SAMPLE_LENGTH2=length(SAMPLE2);
        for bam_ix=1:SAMPLE_LENGTH1
            reads1=[reads1;READ_SET{SAMPLE1(bam_ix)}];
        end
        
        for bam_ix=1:SAMPLE_LENGTH2
            reads2=[reads2;READ_SET{SAMPLE2(bam_ix)}];
        end
        [PV, INFO]= rDiff_mmd(CFG,reads1,reads2);
        P_VALS{i,COUNTER}={PV, INFO};
        clear reads1
        clear reads2
        COUNTER=COUNTER+1;
    end
    
    if CFG.perform_nonparametric
        [PV, INFO]= rDiff_nonparametric(CFG,READ_SET(SAMPLE1),READ_SET(SAMPLE2),variance_function_nonparametric_1, variance_function_nonparametric_2);
        P_VALS{i,COUNTER}={PV, INFO};
        COUNTER=COUNTER+1;
    end
    
end
fprintf('\n')
%Save the p-values
OUT_FILENAME=[CFG.outfile_prefix '.mat'];
save(OUT_FILENAME,'P_VALS')

