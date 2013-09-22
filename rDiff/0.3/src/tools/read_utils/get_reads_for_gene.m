function [reads1] = get_reads_for_gene(CFG,gene)
  

% Get the reads from the bam-file
if strcmp(gene.strand,'-')
    [mask1, read_intron_list] = get_reads(CFG.curr_bamfile, [CFG.chr_prefix gene.chr],gene.start+1,gene.stop+1, '0');
else 
    [mask1, read_intron_list] = get_reads(CFG.curr_bamfile, [CFG.chr_prefix gene.chr],gene.start,gene.stop, '0');
end

% Bring the reads into a matrix form
if isempty(mask1)
    reads1=zeros(0,gene.stop-gene.start+1);
    return
else
    reads1=sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
end

% remove reads which are shorter than CFG.min_read_length
reads1=reads1(sum(reads1,2)>CFG.min_read_length,:);

if isempty(reads1)
    return
end
%remove reads which could stem from other genes
[reads1,FLAG]=remove_reads_from_other_genes(reads1,gene); 

%Subsample reads
if isempty(reads1)
    return
end
if CFG.rDiff_subsample>0
    if size(reads1,1)>CFG.rDiff_subsample
        RP=randperm(size(reads1,1));
        reads1=reads1(RP(1:CFG.rDiff_subsample),:);
    end
end

%Clip reads
if isempty(reads1)
    return
end
if CFG.bases_to_clip>0
    CFG.bases_to_clip=3;
    [reads1]=clip_reads(reads1,CFG.bases_to_clip); 
end

