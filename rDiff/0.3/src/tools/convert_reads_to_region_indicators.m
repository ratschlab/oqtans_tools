function [NEW_READS,UNEXPLAINED_READS,UNEXPLAINED_INDEX] = convert_reads_to_region_indicators(READS, gene)
% Convert the reads into counts of EIRS 

%UNEXPLAINED_REGIONS
UNEXPLAINED_REGIONS = 1:(length(gene.splicingevents) - 1);
UNEXPLAINED_REGIONS(gene.sequence == 0);

%Extend GRAPHNODES to also include introns
EXT_GRAPHNODES = zeros(size(gene.graphnodes,1), length(gene.splicingevents) - 1);
EXT_GRAPHNODES(:, gene.unique_new_exons) = gene.graphnodes;
% This puts gene.graphnodes into the exonic positions of EXT_GRAPHNODES

%Mark for each read into which region it falls
TEMP_READS = sparse(zeros(size(READS,1), length(gene.splicingevents) - 1));
for i = 1:(length(gene.splicingevents) - 1)
  TEMP_READS(:,i) = sum(READS(:, gene.splicingevents(i):(gene.splicingevents(i + 1) - 1)), 2) > 0;
end

%UNEXPLAINED_INDEX = (((sum(EXT_GRAPHNODES,1) == 0) * TEMP_READS1') > 0);
UNEXPLAINED_INDEX = ((( gene.sequence== 0) * TEMP_READS') > 0);
UNEXPLAINED_READS = TEMP_READS(UNEXPLAINED_INDEX,:);

TEMP_READS = TEMP_READS(not(UNEXPLAINED_INDEX),:);

%find the row of EXT_GRAPHNODES minimizing the mismatch between the rows of EXT_GRAPHNODES
%and the reads
%[MAX_VAL,MAX_NODE] = max( (1+diag(1./sum(EXT_GRAPHNODES,2))) * (EXT_GRAPHNODES * TEMP_READS'),[],1);
[MAX_VAL,MAX_NODE] = max(EXT_GRAPHNODES*TEMP_READS'+diag(1./sum(EXT_GRAPHNODES,2))*EXT_GRAPHNODES*TEMP_READS',[],1);

%Create sparse read matrix
NEW_READS = spconvert([[1,1,1; size(EXT_GRAPHNODES,1),2,1]; [MAX_NODE',(3:(length(MAX_NODE)+2))', ones(size(MAX_VAL))']])';
NEW_READS = NEW_READS(3:end,:);


return

