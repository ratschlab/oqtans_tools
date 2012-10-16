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

deseq_config;

CFG.paths = {sprintf('%s/tools', DESEQ_PATH), sprintf('%s/mex', DESEQ_PATH), sprintf('%s', DESEQ_SRC_PATH)} ;
CFG.out_base = outfile;

%%% read list of replicate groups from variable length argument list
rg_list = cell(1,size(varargin, 2));
file_list = cell();
file_cond_ids = [];
file_rep_ids = [];
varargin
for idx = 1:size(varargin, 2)
    rg_list(idx) = varargin(idx);
end
idx = strmatch('', rg_list, 'exact');
rg_list(idx) = [];
for idx = 1:length(rg_list),
    items = separate(rg_list{idx}, ':');
    for idx2 = 1:length(items)
		if isempty(deblank(items{idx2})),
			continue;
		end;
        file_list{end + 1} = items{idx2};
        file_cond_ids(end + 1) = idx; 
        file_rep_ids(end + 1) = idx2; 
    end;
end;
clear idx idx2;

%%%% adapt to number of input arguments
file_num = length(file_list);
RESULTS = cell(1, file_num);

%%% get annotation file
load(sprintf('%s', anno_dir));

%%%% mask overlapping gene regions -> later not counted
[genes] = mask_dubl(genes,0);

%%% remove genes with no annotated exons or where no 
idx = find(arrayfun(@(x)(~isempty(x.exons)*~isempty(x.start)*~isempty(x.stop)), genes));
fprintf('removed %i of %i genes, which had either no exons annotated or lacked a start or stop position\n', size(genes, 2) - size(idx, 2), size(genes, 2))
genes = genes(idx);
clear idx;

%%% check if genes have field chr_num
if ~isfield(genes, 'chr_num')
	chrms = unique({genes(:).chr});
	for i = 1:length(genes)
		genes(i).chr_num = strmatch(genes(i).chr, chrms, 'exact');
	end;
end;

JB_NR = 1;

%%% iterate over all given bam files
for f_idx = 1:file_num  
    CFG.expr1_bam = fullfile('', file_list{f_idx});
    STAT = cell(size(genes, 2),1);

    CFG.RUN = f_idx;

    %%% build up a chromosome mapping list between read files and annotation
	[stat, res] = unix(sprintf('%s/samtools view -H %s', SAMTOOLS_DIR, CFG.expr1_bam));
	chr_names = separate(res, sprintf('\n'));
	for it_idx = 1:length(chr_names),
	    chr_names{it_idx} = regexprep(chr_names{it_idx}, '.*SN:(.+)\s.*', '$1');
	end;
    CFG.chr_names = chr_names;

  	%%%% configuration
  	%CFG = configure_difftest(CFG);

  	%%%% load genes
  	nb_of_chunks = CFG.nb_of_chunks;
  	idx = [(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/size(genes,2))'];

  	%%% submit jobs to cluster
  	for i = 1:nb_of_chunks
    	%%% setting parameters
    	PAR.genes = genes(idx(idx(:,2)==i,1));
    	CFG.rproc_memreq = 10000;
    	CFG.rproc_par.mem_req_resubmit = [ones(1,7)*10000 12000 20000 32000];
    	CFG.rproc_par.identifier = sprintf('DSq-%i-',JB_NR);
    	CFG.JB_NR = JB_NR;
    	PAR.CFG = CFG;
        fprintf(1, 'Submitting job %i to cluster\n', JB_NR);
    	JOB_INFO(JB_NR) = rproc('diff_expr_caller_counts', PAR, CFG.rproc_memreq,CFG.rproc_par, 20);
	    %diff_expr_caller_counts(PAR);
    	JB_NR = JB_NR + 1;
    end;
end;
fprintf(1, '\nWaiting for cluster jobs to finish.\n');

pause(10);

[JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1, 1);
if num_crashed>0, pause(60); end; % some delay to wait until results are written

fprintf(1, '\nCluster jobs finished.\nCollect results.\n');

JB_NR = 1;
TASK_NR = 1;
nb_of_chunks = CFG.nb_of_chunks;
ERR=[];

idx = [(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/ size(genes,2))'];

for f_idx = 1:file_num,
  CFG.IDX1 = f_idx;
  PROBS = cell(size(genes));
  % collect results
  for i = 1:nb_of_chunks
    LOAD_DIR=[CFG.out_base, '_exon_counts_' int2str(JB_NR)   '.mat'];
%try
      while ~fexist(LOAD_DIR), fprintf('waiting for file %s to appear\n', LOAD_DIR); pause(1);  end ;

      load (LOAD_DIR,'STAT')
      temp_idx=idx(idx(:,2)==i,1);
      for j=1:size(temp_idx)
        PROBS{temp_idx(j)} = STAT{j};
      end
%    catch
%      ERR=[ERR;JB_NR];
%    end
    JB_NR=JB_NR+1;
  end
  RESULTS{TASK_NR}=PROBS;
  TASK_NR=TASK_NR+1;
end

ERR

S=size(genes,2);

READCOUNTS_ALL=zeros(S, file_num);
READCOUNTS_EXON=zeros(S, file_num);
LENGTH_ALL=zeros(S,file_num);
LEN_EXON=zeros(S, file_num);

for j=1:file_num,
  for i=1:S
    T=RESULTS{j}{i};
    if isempty(T)
        warning('result %i %j empty', j, i) ;
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
fid_conditions = fopen(sprintf('%s_CONDITIONS.tab', outfile), 'w');
fid_counts = fopen(sprintf('%s_COUNTS.tab', outfile) ,'w');

fprintf(fid_counts,'gene');
fprintf(fid_conditions, 'file\tcondition\treplicate\n');

for j = 1:length(file_list)
    fname = file_list{j} ;
    fname = separate(fname, '/');
    fname = fname{end};
    fname = strrep(fname, '.bam', '') ;
    fprintf(fid_counts,'\t%s', fname);
    
    fprintf(fid_conditions, '%s\t%i\t%i\n', fname, file_cond_ids(j), file_rep_ids(j));
end;
fprintf(fid_counts,'\n') ;

for i = 1:size(genes,2)
    fprintf(fid_counts,'%s',genes(i).name);
    for j = 1:length(file_list),
      fprintf(fid_counts,'\t%i', READCOUNTS_EXON(i,j));
    end
    fprintf(fid_counts,'\n');
end
fclose(fid_counts);
fclose(fid_conditions);

% clean up
JB_NR = 1;
for f_idx = 1:file_num,
    for i = 1:nb_of_chunks
        delete([CFG.out_base, '_exon_counts_' int2str(JB_NR)   '.mat']);
		JB_NR = JB_NR + 1;
    end
end

exit;
