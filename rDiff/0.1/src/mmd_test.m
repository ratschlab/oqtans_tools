function [p_values]=mmd_test(genes,bam1,bam2)

min_read_length = 0;

p_values=ones(size(genes,2),1);

for j=1:size(genes,2)
  fprintf(1, 'gene %i: %s\n', j, genes(j).name);
  gene=genes(j);
  eidx = [gene.start,gene.stop];
 
    [mask1, read_intron_list] = get_reads(bam1, gene.chr,gene.start,gene.stop, '0');
    [mask2, read_intron_list] = get_reads(bam2, gene.chr,gene.start,gene.stop, '0');  

	%% read coverage is empty set p-value as 1 and continue
	if isempty(mask1) || isempty(mask2)
		p_values(j)=[1];
		continue
	end

    %bring both readmask to the same size and filter for read length
    max11=max(mask1(1,:));
    max12=max(mask1(2,:));
    max21=max(mask2(1,:));
    max22=max(mask2(2,:));
    if max22<max12     
      reads1=(spconvert([[mask1',ones(size(mask1,2),1)];max11+1,1,-1;max11+2,gene.stop-gene.start+1,-1]));
      reads1=reads1(sum(reads1,2)>min_read_length,:);
      reads2=(spconvert([[mask2',ones(size(mask2,2),1)];max21+1,gene.stop-gene.start+1,-1;max21+2,1,-1]));
      reads2=reads2(sum(reads2,2)>min_read_length,:);    
    end
    if max22>max12
      reads1=(spconvert([[mask1',ones(size(mask1,2),1)];max11+1,gene.stop-gene.start+1,-1;max11+2,1,-1]));
      reads1=reads1(sum(reads1,2)>min_read_length,:);
      reads2=(spconvert([[mask2',ones(size(mask2,2),1)];max21+1,1,-1;max21+2,gene.stop-gene.start+1,-1]));
      reads2=reads2(sum(reads2,2)>min_read_length,:);    
    end
    if max22==max12
      reads1=(spconvert([[mask1',ones(size(mask1,2),1)];max11+1,gene.stop-gene.start+1,-1;max11+2,1,-1]));
      reads1=reads1(sum(reads1,2)>min_read_length,:);
      reads2=(spconvert([[mask2',ones(size(mask2,2),1)];max21+1,1,-1;max21+2,gene.stop-gene.start+1,-1]));
      reads2=reads2(sum(reads2,2)>min_read_length,:);    
    end
 
       
    [P_VALUE, STRUCT]= diff_mmd(reads1,reads2,gene);
    p_values(j)=[P_VALUE];
end
