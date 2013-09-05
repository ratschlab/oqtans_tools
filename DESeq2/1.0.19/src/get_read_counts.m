function get_read_counts(anno_dir, outfile, varargin)
%
% -- input --
% anno_dir: directory of genes
% outfile: output file 
% varargin: list of BAM files (at least two)

% DESeq paths
global DESEQ2_PATH DESEQ2_SRC_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% SAMTools path
global SAMTOOLS_DIR

%%%% paths
addpath(sprintf('%s/tools', DESEQ2_PATH));
addpath(sprintf('%s/mex', DESEQ2_PATH));
addpath(sprintf('%s', DESEQ2_SRC_PATH));

deseq_config;

%%% read list of replicate groups from variable length argument list
rg_list = cell(1,size(varargin, 2));
file_list = cell();
file_cond_ids = [];
file_rep_ids = [];
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

%%%%% adapt to number of input arguments
file_num = length(file_list);
RESULTS = cell(1, file_num);

%%%% get annotation file
load(sprintf('%s', anno_dir));

%%%%% mask overlapping gene regions -> later not counted
[genes] = mask_dubl(genes,0);

%%%% remove genes with no annotated exons or where no 
idx = find(arrayfun(@(x)(~isempty(x.exons)*~isempty(x.start)*~isempty(x.stop)), genes));
fprintf('removed %i of %i genes, which had either no exons annotated or lacked a start or stop position\n', size(genes, 2) - size(idx, 2), size(genes, 2))
genes = genes(idx);
clear idx;

%%%% check if genes have field chr_num
if ~isfield(genes, 'chr_num')
	chrms = unique({genes(:).chr});
	for i = 1:length(genes)
		genes(i).chr_num = strmatch(genes(i).chr, chrms, 'exact');
	end;
end;

%%%% iterate over all given bam files
for f_idx = 1:file_num  
    expr1_bam = fullfile('', file_list{f_idx});
    STAT = cell(size(genes, 2),1);
    for i=1:size(genes,2)
        RESULT = cell(1,7);
        gene = genes(i);
        RESULT{4} = f_idx; 
        RESULT{1} = gene.name;
        if isempty(gene.exons)
            RESULT{2} = inf;
            RESULT{3} = inf;
            RESULT{5} = [inf,inf];
            STAT{i} = RESULT;
            continue;
        elseif or(isempty(gene.start),isempty(gene.stop))
		    RESULT{2} = inf;
		    RESULT{3} = inf;
		    RESULT{5} = [inf,inf];
		    STAT{i} = RESULT;
		    continue;
        end
        if ~isempty(gene.chr_num),
		    [mask1, read_intron_list] = get_reads(expr1_bam, gene.chr, gene.start, gene.stop, '0');
		    clear read_intron_list; 
	    else
		    mask1 = []; 
        end;

        if isempty(mask1)
            reads1 = zeros(0,gene.stop-gene.start+1);
        else
            reads1 = sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
        end
	    if ~isempty(reads1);
		    [reads1,FLAG] = remove_reads_from_other_genes(reads1,gene);
	    end
        L = size(reads1);
        RESULT{2}=[size(reads1,1)]; % number of all reads falling in that gene
        EXON_IDX=zeros(1,gene.stop-gene.start+1);
        for t=1:size(gene.transcripts,2)
            for e=1:size(gene.exons{t},1)
                EXON_IDX((gene.exons{t}(e,1)-gene.start+1):(gene.exons{t}(e,2)-gene.start+1))=1;
            end
        end
        reads1 = reads1(sum(reads1(:,find(EXON_IDX)),2)>0,:); 
        L1 = sum(EXON_IDX);
        RESULT{3}=[size(reads1,1)]; % number of reads overlapping to exons
        RESULT{5}=[L, L1]; % size of reads1, number of exonic positions
        % old and weighted poisson new ,weighted regions reads and
        % unexplained reads
        clear reads1;
        STAT{i} = RESULT;
    end;
    RESULTS{f_idx} = STAT;
end;

S=size(genes,2);
READCOUNTS_ALL=zeros(S, file_num);
READCOUNTS_EXON=zeros(S, file_num);
LENGTH_ALL=zeros(S,file_num);
LEN_EXON=zeros(S, file_num);

for j=1:file_num,
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

%%%%% write results for all bam files
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
exit;
