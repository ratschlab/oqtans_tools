function STAT = diff_expr_caller_counts(PAR)
% genes = opt_transcripts_caller(PAR)
%
% -- input --
% PAR contains
%     CFG: configuration struct
%     genes: struct defining genes with start, stops, exons etc.
%     profile_weights: weights of profile functions
%     intron_dists: distances to closest intron
%
% -- output --
% genes: struct with additional fields of eg. estimated transcript weights

CFG = PAR.CFG;
JB_NR = CFG.JB_NR;

%%%%% define globals

%These are the non-missing regions
genes = PAR.genes;
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
    load_only = false;
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
         
   [mask1] = get_reads(expr1_bam, gene.chr,gene.start, gene.stop, '0');
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
    reads1=reads1(sum(reads1(:,find(EXON_IDX)),2)>0,:); 
    L1=sum(EXON_IDX);
        
    RESULT{3}=[size(reads1,1)]; % number of reads overlapping to exons
    RESULT{5}=[L,L1]; % size of reads1, number of exonic positions
    % old and weighted poisson new ,weighted regions reads and
    % unexplained reads
    clear reads1;
    STAT{i}=RESULT;
end;

OUT_FILENAME=[CFG.out_base, '_exon_counts_' int2str(JB_NR)   '.mat']
save(OUT_FILENAME,'STAT')

