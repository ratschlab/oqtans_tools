function CFG = configure_rdiff(CFG)
% configure_rdiff(CFG)

%%% rDiff parameters %%%

% Give the filenames of the bam-files to be considered
CFG.BAM_FILES={'condition_A_replicate_1.bam','condition_A_replicate_2.bam','condition_B_replicate_1.bam','condition_B_replicate_2.bam'};

%Name of the experiment. Use the FILENAMES if the entries are empty.
CFG.NAMES={'A1','A2','B1','B2'};


% Give the directory where the bam-files are
CFG.data_dir = '' ;

% Indicate to which sample the bam-files belong
CFG.SAMPLES=[1,1,2,2];

% Location of the gene structure
CFG.genes_path='' ;

% Output directory
CFG.out_base =  '' ;

% Output directory for temporary files
CFG.out_base_temp =  '' ;

%Length of the reads
CFG.sequenced_length=75;

% Prefix for the chromosome name when getting geetting reads from
% the bam-files
CFG.chr_prefix='';

%%% Read filters %%%

% Minimal read length
CFG.min_read_length=30;



%%% Parameters for gene expression estimation
%Count the number of reads ( CFG.estimate_gene_expression=1 for yes
%give the Files for the expresison in CFG.GENE_EXPR_FILES
CFG.estimate_gene_expression=1;

% Use the following files in CFG.GENE_EXPR_FILES for the
% gene_expression. Those must be Tab-delimitered files where each
% line contains the gene name folowed by the expressiob
CFG.Counts_gene_expression='';
CFG.Counts_rDiff_parametric='';
CFG.Counts_rDiff_nonparametric='';



%%% Parameters for variance function

% Use a parametric form for the variance function for sample 1: sigma= a + bx + cx^2
% (CFG.predefined_variance_function1=[] if not; CFG.predefined_variance_function1=[a,b,c] otherwise)
% If CFG.predefined_variance_function1=[a,b,c] is given, the other
% parameters for the variance function estimations are ignored for
% sample 1
CFG.predefined_variance_function1=[];

% Use a parametric form for the variance function for sample 2: sigma= a + bx + cx^2
% (CFG.predefined_variance_function2=[] if not; CFG.predefined_variance_function2=[a,b,c] otherwise)
% If CFG.predefined_variance_function2=[a,b,c] is given, the other
% parameters for the variance function estimations are ignored 
% for sample 2
CFG.predefined_variance_function2=[];

% compute variance function for sample 1 ( 1 = yes , 0 = use precomputed
% variance function saved under CFG.variance_function_1) 
CFG.compute_variance_function_1=1;
CFG.variance_function_1='';
CFG.save_variance_function_1='variance_function_1.mat';

% compute variance function for sample 2 ( 1 = yes , 0 = use precomputed
% variance function saved under CFG.variance_function2) 
CFG.compute_variance_function_2=1;
CFG.variance_function_2='';
CFG.save_variance_function_2='variance_function_2.mat';

% subsample points for the variance function estimate for rDiff.nonparametric
CFG.rDiff_nonparametric_subsample_variance_estimation=10000;	

% Subsample the mean-variance pairs to increas the speed of the
% local regression.CFG.variance_samplebins is the number of bins to
% use and CFG.variance_samples_per_bin is how many samples should
% be drwan per bin
CFG.variance_samplebins=100;
CFG.variance_samples_per_bin=500;



%%% Testing parameters %%%

% subsample reads down to rDiff.subsample to increase speed ( If no
% subsampling shall be done set CFG.rDiff_subsample to 0
CFG.rDiff_subsample=10000;

% Clib the first CFG.bases_to_clip bases at the end of the reads
CFG.bases_to_clip=3;

%Number of bootraps for nonparametric test
CFG.bootstraps=1000;

%Number of bins for variance matching
CFG.nr_of_slices=10;

% Tests to perform
CFG.perform_nonparametric=0;
CFG.perform_parametric=0;
CFG.perform_mmd=0;
CFG.perform_poisson=0;


%%%%% rproc settings %%%%%
CFG.use_rproc = 0; % 1: cluster submission or 0: locally
if CFG.use_rproc,
    CFG.rproc_num_jobs              = 100;
    CFG.rproc_memreq                = 8000;
    CFG.rproc_par.priority          = 55;
    CFG.rproc_par.resubmit          = 3;
    CFG.rproc_par.mem_req_resubmit  = [ 24000 40000 60000];
    CFG.rproc_par.time_req_resubmit = [60*60 100*60 90*60];
    CFG.rproc_par.express           = 0;
    CFG.rproc_par.immediately_bg    = 0;
    CFG.rproc_par.immediately       = 0;
    CFG.rproc_par.arch              = 64;
    CFG.rproc_par.identifier        = '';
    CFG.rproc_par.verbosity         = 0;
    CFG.rproc_time                  = 15*60; % mins
else
    CFG.rproc_num_jobs              = 1;
end





