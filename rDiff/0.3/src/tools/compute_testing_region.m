function [genes]=compute_testing_region(CFG,genes)

%Iterate over gene

for i=1:size(genes,2)
    [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(genes(i));
    genes(i).splicingevents=SPLICINGEVENTS;
    genes(i).sequence=SEQUENCE;
    genes(i).exonsequence=EXONSEQUENCE;
    
    [UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ]=transform_single_end_reads(genes(i),CFG.sequenced_length-CFG.bases_to_clip*2);
    genes(i).unique_new_exons=UNIQUE_NEW_EXONS;
    genes(i).graphnodes=GRAPHNODES;
    genes(i).order_of_graphnodes=ORDER_OF_GRAPHNODE;
    genes(i).eirs_in_seq=EIRS_IN_SEQ;
 
end
