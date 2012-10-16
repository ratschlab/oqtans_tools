function [STAT, dummy] = diff_expr_caller_counts(PAR)

CFG = PAR.CFG;
JB_NR = CFG.JB_NR;

dummy=[] ;

%%%%% define globals

%% rDiff paths
%global CFG.DESEQ_PATH CFG.DESEQ_SRC_PATH

%% interpreter paths
%global CFG.INTERPRETER CFG.MATLAB_BIN_PATH CFG.OCTAVE_BIN_PATH

%% SAMTools path
%global CFG.SAMTOOLS_DIR

%These are the non-missing regions
genes = PAR.genes;
for i=1:length(genes),
	genes(i).start=min(genes(i).exons{1}(:)) ;
	genes(i).stop=max(genes(i).exons{1}(:)) ;
        for j=2:length(genes(i).exons)
		if genes(i).start>min(genes(i).exons{j}(:)),
		     genes(i).start=min(genes(i).exons{j}(:)) ;
                end
		if genes(i).stop<max(genes(i).exons{j}(:)),
		     genes(i).stop=max(genes(i).exons{j}(:)) ;
                end
	end 
end
clear PAR;

%%%% paths
for p = 1:length(CFG.paths),
	addpath(CFG.paths{p});
end

expr1_bam = CFG.expr1_bam;

STAT = cell(size(genes));
do_save = 0;

for i=1:size(genes,2)
    RESULT = cell(1,7);
    gene = genes(i);
        
    RESULT{4} = CFG.RUN;
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

    if ~isempty(gene.chr_num) && (gene.chr_num <= length(CFG.chr_names)),
		CHR = CFG.chr_names{gene.chr_num};
		[mask1, read_intron_list] = get_reads(expr1_bam, CHR, gene.start, gene.stop, '0');
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

OUT_FILENAME=[CFG.out_base, '_exon_counts_' int2str(JB_NR)   '.mat'];
save(OUT_FILENAME,'STAT')

