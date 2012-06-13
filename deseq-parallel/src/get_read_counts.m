function get_read_counts(anno_dir, outfile, varargin)
% get_read_counts(ANNO, BAM_FILES)
%
% -- input --
% anno_dir: directory of genes
% outfile: output file 
% varargin: list of BAM files (at least two)

% DESeq paths
global DESEQ_PATH DESEQ_SRC_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% SAMTools path
global SAMTOOLS_DIR

% configuration 
global CFG

%%%% paths
addpath(sprintf('%s/tools', DESEQ_PATH));
addpath(sprintf('%s/mex', DESEQ_PATH));
addpath(sprintf('%s', DESEQ_SRC_PATH));
addpath('/mnt/galaxyTools/tools/mtimpred/src/rproc');
addpath('/mnt/galaxyTools/tools/mtimpred/src/utils');

CFG.paths = {sprintf('%s/tools', DESEQ_PATH), sprintf('%s/mex', DESEQ_PATH), sprintf('%s', DESEQ_SRC_PATH), '/mnt/galaxyTools/tools/mtimpred/src/utils'} ;
CFG.out_base = outfile;

chr1_flag = 1;

%%% read list of bam files from variable length argument list
file_list = cell(1,size(varargin, 2));
for idx = 1:size(varargin, 2)
    file_list(idx) = varargin(idx);
end
file_list=sort(file_list) ;
idx=strmatch('', file_list, 'exact') ;
file_list(idx)=[] ;
clear idx;
BAM_FILES = file_list;

%%%% adapt to number of input arguments
RUNS=1:length(BAM_FILES);
NR_OF_TASK=size(RUNS,2);

RESULTS=cell(1,NR_OF_TASK);

%%% get annotation file
load(sprintf('%s', anno_dir));


%%%% get overlapping genes -> later not counted
[genes]=mask_dubl(genes,0);

%%% remove genes with no annotated exons or where no 
idx = find(arrayfun(@(x)(~isempty(x.exons)*~isempty(x.start)*~isempty(x.stop)), genes));
fprintf('removed %i of %i genes, which had either no exons annotated or lacked a start or stop position', size(genes, 2) - size(idx, 2), size(genes, 2))
genes = genes(idx);
clear idx;

JB_NR = 1;

%%% iterate over all given bam files
for RUN=1:size(RUNS,2)  
    CFG.expr1_bam = fullfile('', BAM_FILES{RUN});
    STAT=cell(size(genes, 2),1);

    CFG.IDX1 = RUNS(1,RUN);

  	%%%% load genes
  	nb_of_chunks=CFG.nb_of_chunks;
  	idx=[(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/size(genes,2))'];

  	%%% submit jobs to cluster
  	for i = 1:nb_of_chunks
    	%%% setting parameters
    	PAR.genes = genes(idx(idx(:,2)==i,1));
    	CFG.rproc_memreq = 10000;
    	CFG.rproc_par.mem_req_resubmit = [ones(1,7)*10000 12000 20000 32000];
    	CFG.rproc_par.identifier = sprintf('DSq-%i-',JB_NR);
    	CFG.JB_NR = JB_NR;
		CFG.RUN = RUN;
    	PAR.CFG = CFG;
		fprintf(1, 'Submitting job %i to cluster\n', JB_NR);
    	JOB_INFO(JB_NR) = rproc('diff_expr_caller_counts', PAR, CFG.rproc_memreq,CFG.rproc_par, 20);
    	JB_NR = JB_NR + 1;
    end 
end

pause(100)
[JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1, 0);
if num_crashed>0, pause(60); end; % some delay to wait until results are written

JB_NR = 1;
TASK_NR = 1;
nb_of_chunks = CFG.nb_of_chunks;
ERR=[];

idx=[(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/ size(genes,2))'];

for RUN=1:size(RUNS,2)
  CFG.IDX1 = RUNS(1,RUN);
  PROBS = cell(size(genes));
  % collect results
%  fprintf(1, 'collecting results for RUN %i', RUN) ;
  for i = 1:nb_of_chunks
    LOAD_DIR=[CFG.out_base, '_exon_counts_' int2str(JB_NR)   '.mat'];
%    fprintf(1, 'collecting results from %s\n', LOAD_DIR) ;
    try
      load (LOAD_DIR,'STAT')
      temp_idx=idx(idx(:,2)==i,1);
      for j=1:size(temp_idx)
        PROBS{temp_idx(j)} = STAT{j};
      end
    catch
      ERR=[ERR;JB_NR];
    end
    JB_NR=JB_NR+1;
  end
  RESULTS{TASK_NR}=PROBS;
  TASK_NR=TASK_NR+1;
end


S=size(genes,2);

READCOUNTS_ALL=zeros(S,NR_OF_TASK);
READCOUNTS_EXON=zeros(S,NR_OF_TASK);
LENGTH_ALL=zeros(S,NR_OF_TASK);
LEN_EXON=zeros(S,NR_OF_TASK);

for j=1:NR_OF_TASK
  for i=1:S
    T=RESULTS{j}{i};
    if isempty(T)
        continue
    else
        READCOUNTS_ALL(i,j)=T{2};
        READCOUNTS_EXON(i,j)=T{3};          
        LENGTH_ALL(i,j)=T{5}(1);
        LEN_EXON(i,j)=T{5}(2);
    end
  end
end

%%%% write results for all bam files

%MIN_RPKM=30 ;
MIN_RPKM=0 ;

READSUM=(sum(READCOUNTS_EXON, 1)./1000000);

fid_counts=fopen(sprintf('%s_COUNTS.tab', outfile) ,'w');
%fid_rpkm=fopen(sprintf('%s_RPKM.tab', outfile) ,'w');

fprintf(fid