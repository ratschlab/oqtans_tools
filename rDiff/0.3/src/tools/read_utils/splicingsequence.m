function [SPLICINGEVENTS, SEQUENCE, EXONSEQUENCE, IDENTIFICATIONLENGTH] = splicingsequence(GENE)
% This function generates all sequence of all splicesites
% SPLICINGEVENTS, a SEQUENCE marking which elements in
% SPLICINGEVENTS are introns and exons and a sequence for each
% transcript which indicates which Exons are in cluded in this transcript
% SEQUENCE contains indices of SPLICINGEVENTS, containing 1, if the position right to the idx in the transcript
% are exonic, 0 otherwise

  EXONS = GENE.exons;
  START = GENE.start;
  STOP = GENE.stop;
  
  NB_OF_TRANSCR = size(EXONS,2);
  
  SPLICINGEVENTS = [];
  for i = 1:NB_OF_TRANSCR
    SPLICINGEVENTS = [SPLICINGEVENTS, EXONS{i}(:,1)', EXONS{i}(:,2)' ];
  end
  SPLICINGEVENTS = SPLICINGEVENTS - START + 1;
  SPLICINGEVENTS = unique(SPLICINGEVENTS);
  
  EXONSEQUENCE = zeros(NB_OF_TRANSCR, length(SPLICINGEVENTS) - 1);
  
  for i = 1:NB_OF_TRANSCR
    %%% for every exon in transcript i
    for j = 1:size(EXONS{i}, 1)
      POS_START = find(SPLICINGEVENTS == EXONS{i}(j,1) - START + 1, 1, 'first');
      POS_STOP = find(SPLICINGEVENTS < EXONS{i}(j,2) - START + 1, 1, 'last');
      %%% set splicing events active in transcript i to 1
      EXONSEQUENCE(i, POS_START:POS_STOP) = 1;
    end 
  end

  %%% merge active splicing events of all transcripts
  SEQUENCE = max(EXONSEQUENCE, [], 1);
  IDENTIFICATIONLENGTH = [];

