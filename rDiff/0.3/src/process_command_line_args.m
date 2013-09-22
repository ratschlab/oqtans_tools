function [CFG]=process_command_line_args(CFG,ARGS)
% [CFG]=process_command_line_args(CFG,ARGS)

%This function adds the command line arguments to the CFG config variable

%Parse the ARGS 

%Check wether octave or matlab is used
if size(ver('Octave'),1)
  INTERPR = 1;
else
  INTERPR = 0;
end

%turn of warning
if INTERPR
    warning('off', 'Octave:precedence-change');
    warning('off', 'Octave:function-name-clash');
    warning('off', '');
    warning('off', 'Octave:num-to-str');
    warning('off', 'Octave:function-name-clash');
    warning('off', 'Octave:divide-by-zero');
    warning('off', 'Octave:future-time-stamp');
    warning('off', 'solve_qp:constraints');
    warning('off', 'Octave:assign-as-truth-value');
    warning('off', 'Octave:matlab-incompatible');
else
   warning('off', 'MATLAB:typeaheadBufferOverflow');
end

%Seperate the Variable Field
if INTERPR
    ARGS = strsplit(ARGS,';');
else
    ARGS=regexp(ARGS,';','split');
end

%Assign the local variables

%iterate over Fields
for i=1:length(ARGS)
    if isempty(ARGS{i})
        continue
    end
    if INTERPR
        VALS = strsplit(ARGS{i},':');
    else
        VALS=regexp(ARGS{i},':','split');
    end
    
    
    if length(VALS)>2
        error([" more than one field for variable: " VALS{1} ":" VALS{:} "\n  Maybe there are colons in the input argument?"])
    end
    
    if strcmp(VALS{1},"RDIFF_RES_DIR"),RDIFF_RES_DIR=VALS{2};continue,end
    if strcmp(VALS{1},"RDIFF_INPUT_DIR"),RDIFF_INPUT_DIR=VALS{2};continue,end
    if strcmp(VALS{1},"BAM_INPUT1"),BAM_INPUT1=VALS{2};continue,end
    if strcmp(VALS{1},"BAM_INPUT2"),BAM_INPUT2=VALS{2};continue,end
    if strcmp(VALS{1},"GFF_INPUT"),GFF_INPUT=VALS{2};continue,end
    if strcmp(VALS{1},"READ_LENGTH"),READ_LENGTH=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"MIN_READ_LENGTH"),MIN_READ_LENGTH=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"EST_GENE_EXPR"),EST_GENE_EXPR=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"ONLY_GENE_EXPR"),ONLY_GENE_EXPR=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"VAR_PATH1"),VAR_PATH1=VALS{2};continue,end
    if strcmp(VALS{1},"VAR_PATH2"),VAR_PATH2=VALS{2};continue,end
    if strcmp(VALS{1},"SAVE_VAR1"),SAVE_VAR1=VALS{2};continue,end
    if strcmp(VALS{1},"SAVE_VAR2"),SAVE_VAR2=VALS{2};continue,end
    if strcmp(VALS{1},"PRED_VAR1"),PRED_VAR1=VALS{2};continue,end
    if strcmp(VALS{1},"PRED_VAR2"),PRED_VAR2=VALS{2};continue,end
    if strcmp(VALS{1},"ONLY_GENE_START"),ONLY_GENE_START=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"SUBSAMPLE"),SUBSAMPLE=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"CLIP"),CLIP=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"BOOTSTRAP"),BOOTSTRAP=str2num(VALS{2});continue,end
    if strcmp(VALS{1},"TEST_METH_NAME"),TEST_METH_NAME=VALS{2};continue,end
    if strcmp(VALS{1},"MERGE_SAMPLE"),MERGE_SAMPLE=str2num(VALS{2});continue,end
end

%Process Bamfiles
if INTERPR
    BAMS1 = strsplit(BAM_INPUT1,',');
    BAMS2 = strsplit(BAM_INPUT2,',');
else
    BAMS1=regexp(BAM_INPUT1,',','split');
    BAMS2=regexp(BAM_INPUT2,',','split');
end

CFG.BAM_FILES={BAMS1{:},BAMS2{:}};
    
%Name of the experiment. Use the FILENAMES if the entries are empty.
CFG.NAMES=CFG.BAM_FILES;
for i=1:length(CFG.NAMES)
	CFG.NAMES{i}=strrep(CFG.NAMES{i},"/","_");
end

% Give the directory where the input-files are
CFG.data_dir = [RDIFF_INPUT_DIR '/'];

% Indicate to which sample the bam-files belong
CFG.SAMPLES=[repmat(1,1,size(BAMS1,2)),repmat(2,1,size(BAMS2,2))];

%Process directories



% Location of the gene structure
CFG.genes_path=GFF_INPUT;

% Output directory

CFG.out_base =  [RDIFF_RES_DIR '/'];

% Output directory for temporary files
CFG.out_base_temp =   [CFG.out_base '/temp/'];
mkdir(CFG.out_base_temp);


%Check which method to perform
if strcmp(TEST_METH_NAME,'poisson')
    CFG.perform_poisson=1;
end
if strcmp(TEST_METH_NAME,'param')
    CFG.perform_parametric=1;
end
if strcmp(TEST_METH_NAME,'nonparam')
    CFG.perform_nonparametric=1;
end
if strcmp(TEST_METH_NAME,'mmd')
    CFG.perform_mmd=1;
end


%Process arguments for gene expression estimation
CFG.estimate_gene_expression=EST_GENE_EXPR;
CFG.only_gene_expression=ONLY_GENE_EXPR;


%Set options for the variance function
CFG.merge_sample1=MERGE_SAMPLE;
CFG.merge_sample2=MERGE_SAMPLE;

%If samples contains leass than one sample, merge replicates
if length(BAMS1)<2
    CFG.merge_sample1=1;
end
if length(BAMS2)<2
    CFG.merge_sample2=1;
end

%Use predefined  parameters
CFG.predefined_variance_function1=[];
if not(isempty(PRED_VAR1))
    if INTERPR
        VALS = strsplit(PRED_VAR1,',');
    else
        VALS=regexp(PRED_VAR1,',','split');
    end
    for i=1:length(VALS)
        CFG.predefined_variance_function1(end+1)=str2num(VALS{i});
    end
end
CFG.predefined_variance_function2=[];
if not(isempty(PRED_VAR2))
    if INTERPR
        VALS = strsplit(PRED_VAR2,',');
    else
        VALS=regexp(PRED_VAR2,',','split');
    end
    for i=1:length(VALS)
        CFG.predefined_variance_function2(end+1)=str2num(VALS{i});
    end
end


if not(isempty(SAVE_VAR1))
    CFG.save_variance_function_1=SAVE_VAR1;
else
    CFG.save_variance_function_1='variance_function_1.mat';
end
if not(isempty(SAVE_VAR2))
    CFG.save_variance_function_2=SAVE_VAR2;
else
    CFG.save_variance_function_2='variance_function_2.mat';
end

if not(isempty(VAR_PATH1))
    CFG.variance_function_1=VAR_PATH1;
end
if not(isempty(VAR_PATH2))
    CFG.variance_function_2=VAR_PATH2;
end

%Option not implemented yet
CFG.compute_variance_function_1=1;
CFG.compute_variance_function_2=1;

%use only gene starts and stops for rDiff.nonparametric variance
%function esitmation
CFG.only_gene_start=ONLY_GENE_START;

%Process read arguments
CFG.sequenced_length=READ_LENGTH;
CFG.min_read_length=min(CFG.sequenced_length,MIN_READ_LENGTH);

CFG.rDiff_subsample=SUBSAMPLE;
CFG.rDiff_nonparametric_subsample_variance_estimation=CFG.rDiff_subsample;	

CFG.bases_to_clip=CLIP;

%Process arguments for rDiff.nonparametric
CFG.bootstraps=BOOTSTRAP;

return
