function [READS_OUT,FLAG]=remove_reads_from_other_genes(READS,GENE)
%This funtion removes the reads in READS which could ome from other
%annotated genes. FLAG is 1 if this was sucsesfull and 0 otherwise
READS_IN=READS;
if isfield(GENE,'non_unique_regions')
  EXONS=GENE.non_unique_regions;
  IDX=zeros(1,GENE.stop-GENE.start+1);
  
  for i=1:size(EXONS,1)
    START=max(EXONS(i,1),GENE.start)-GENE.start+1;
    STOP=min(EXONS(i,2),GENE.stop)-GENE.start+1;
    IDX(START:STOP)=1;
  end  
  READS=READS(not(sum(READS(:,IDX>0),2)==sum(READS,2)),:); 
  FLAG=1;
  READS_OUT=READS;
else
  READS_OUT=READS_IN;
  FLAG=0;  
end

