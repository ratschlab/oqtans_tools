function [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE,IDENTIFICATIONLENGTH]= ...
    splicingsequence(GENE)
%This function generates all sequence of all splicesites
%SPLICINGEVENTS, a SEQUENCE marking which elements in
%SPLICINGEVENTS are introns and exons and a sequence for each
%transcript which indicates which Exons are in cluded in this transcript

  EXONS=GENE.exons;
  START=GENE.start;
  STOP=GENE.stop;
  %assert(START<STOP,'Start>=Stop')
  
  NB_OF_TRANSCR=size(EXONS,2);
  
  
  SPLICINGEVENTS=[];
  for i=1:NB_OF_TRANSCR
    SPLICINGEVENTS=[SPLICINGEVENTS,EXONS{i}(:,1)'-START+1,EXONS{i}(:,2)'-START+1];
  end
  SPLICINGEVENTS=unique(SPLICINGEVENTS);
  
  
  EXONSEQUENCE=zeros(NB_OF_TRANSCR,length(SPLICINGEVENTS)-1);
  
  for i=1:NB_OF_TRANSCR
    for j=1:size(EXONS{i},1)
      POS_START=find(SPLICINGEVENTS==EXONS{i}(j,1)-START+1,1,'first');
      POS_STOP=find(SPLICINGEVENTS<EXONS{i}(j,2)-START+1,1,'last');
      EXONSEQUENCE(i,POS_START:POS_STOP)=1;
      
    end 
  end

  SEQUENCE=max(EXONSEQUENCE,[],1);