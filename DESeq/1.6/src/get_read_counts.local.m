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

addpath(sprintf('%s/tools', DESEQ_PATH));
addpath(sprintf('%s/mex', DESEQ_PATH));
addpath(sprintf('%s', DESEQ_SRC_PATH));

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
%load(sprintf('%s', anno_dir), '-mat');
load(sprintf('%s', anno_dir));

%%% restrict genes to only one chromosome, if neccessary
%if chr1_flag;
%    idx = find([genes(:).chr_num] == 1);
%    genes = genes(idx);
%    clear idx;
%end
%%%% get overlapping genes -> later not counted
[genes]=mask_dubl(genes,0);

%%% remove genes with no annotated exons or where no 
%idx = find(~isempty(genes(:).exons) && ~isempty(genes(:).start) && ~isempty(genes(:).stop));
idx = find(arrayfun(@(x)(~isempty(x.exons)*~isempty(x.start)*~isempty(x.stop)), genes));
fprintf('removed %i of %i genes, which had either no exons annotated or lacked a start or stop position', size(genes, 2) - size(idx, 2), size(genes, 2))
genes = genes(idx);
clear idx;

%%% iterate over all given bam files
for RUN=1:size(RUNS,2)  
    expr1_bam = fullfile('', BAM_FILES{RUN});
    STAT=cell(size(genes, 2),1);

    for i=1:size(genes,2)
        RESULT=cell(1,7);
        gene = genes(i);
        
        RESULT{4}=RUN;
        RESULT{1}=gene.name;
        load_only = false;
        if isempty(gene.exons)
            %RESULT{2}=[inf,inf];
            %RESULT{3}=[inf,inf];
            %RESULT{5}=inf;
            RESULT{2}=inf;
            RESULT{3}=inf;
            RESULT{5}=[inf,inf];
            STAT{i}=RESULT;
            continue;
        elseif or(isempty(gene.start),isempty(gene.stop))
            %RESULT{2}=[inf,inf];
            %RESULT{3}=[inf,inf];
            %RESULT{5}=inf;
            RESULT{2}=inf;
            RESULT{3}=inf;
            RESULT{5}=[inf,inf];
            STAT{i}=RESULT;
            continue;
        end
         
        [mask1, read_intron_list] = get_reads(expr1_bam, gene.chr, gene.start, gene.stop, '0');
        clear read_intron_list; 
        if isempty(mask1)
            reads1=zeros(0,gene.stop-gene.start+1);
        else
            reads1=sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
        end
		if ~isempty(reads1);
			[reads1,FLAG] = remove_reads_from_other_genes(reads1,gene);
		end
        L=size(reads1);
        RESULT{2}=[size(reads1,1)]; % number of all reads falling in that gene
        
        EXON_IDX=zeros(1,gene.stop-gene.start+1);
        for t=1:size(gene.transcripts,2)
          for e=1:size(gene.exons{t},1)
            EXON_IDX((gene.exons{t}(e,1)-gene.start+1):(gene.exons{t}(e,2)-gene.start+1))=1;
          end
        end
        reads1=reads1(sum(reads1(:,find(EXON_IDX)),2)>0,:); 
        L1=sum(EXON_IDX);
        
        RESULT{3}=[size(reads1,1)]; % number of reads overlapping to exons
        RESULT{5}=[L,L1]; % size of reads1, number of exonic positions
        % old and weighted poisson new ,weighted regions reads and
        % unexplained reads
        clear reads1;
        STAT{i}=RESULT;
    end;
	RESULTS{RUN} = STAT;
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

fprintf(fid_counts,'gene\t');
%fprintf(fid_rpkm,'gene\t');

%maximum=zeros(size(genes,2), 1);

for j=1:length(BAM_FILES)
    fname=BAM_FILES{j} ;
    fname=separate(fname, '/');
    fname=fname{end};
    fname=strrep(fname, '.bam', '') ;
    fprintf(fid_counts,'%s\t', fname);
    %fprintf(fid_rpkm,'%s\t', fname);
end
fprintf(fid_counts,'\n') ;
%fprintf(fid_rpkm,'\n') ;
for i=1:size(genes,2)
    %maximum(i) = max(READCOUNTS_EXON(i,:)./(READSUM.*LEN_EXON(i,:)./1000));
    RPKM = READCOUNTS_EXON(i,:)./(READSUM.*LEN_EXON(i,:)./1000);
    if any(RPKM > MIN_RPKM);
        fprintf(fid_counts,'%s\t',genes(i).name);
        %fprintf(fid_rpkm,'%s\t',genes(i).name);
        for j=1:length(BAM_FILES),
          fprintf(fid_counts,'%i\t', READCOUNTS_EXON(i,j));
          %fprintf(fid_rpkm,'%i\t', RPKM(j));
        end
        fprintf(fid_counts,'\n');
        %fprintf(fid_rpkm,'\n');
    end
end
fclose(fid_counts);
%fclose(fid_rpkm);
exit;
