function [P_VALUE, STRUCT]= diff_poisson_bonf_2(READS1,READS2,GENE,START_OR_COVERAGE)

% Calculates the p-Values of a poisson test on each of the intervalls
% and combines them to one p-values. SEGMENTS is a n x 2 array where
% each row is the start and the end of an exon
%START_OR_COVERAGE tells whether the test should be based on the
%startpoints or the coverage. START_OR_COVERAGE=0 means start =1
%means coverage
if nargin<4
  START_OR_COVERAGE = 1;
end;

if(and(size(READS1,1)>0,size(READS2,1)>0))

r=size(READS1,1);
S=spconvert([1:r;1:r;(1./full(sum(READS1,2)))']');
READS1=S*READS1;

r=size(READS2,1);
S=spconvert([1:r;1:r;(1./full(sum(READS2,2)))']');
READS2=S*READS2;
   


STRUCT=[];

%Determine the relevant segments
[SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(GENE);
k=size(EXONSEQUENCE,1);

if (k>1)
%Marker for the potential differentialy expressed exons
  IDX=find(not(or(sum(EXONSEQUENCE,1)==0,sum(EXONSEQUENCE,1)==k))==1);
  
  
  
  %Get segments
  SEGMENTS=[SPLICINGEVENTS(IDX)',SPLICINGEVENTS(IDX+1)'];
  
  READS_PER_EXON=zeros(2,size(SEGMENTS,1));
  P_VALUES=zeros(1,size(SEGMENTS,1));
  
  if START_OR_COVERAGE==0
    for j=1:size(SEGMENTS,1)
      N1=sum(and((sum(READS1(:,1:(SEGMENTS(j,1)-1)),2)==0),(sum(READS1(:,SEGMENTS(j,1):(SEGMENTS(j,2)-1)),2)>0)));
      N2=sum(and((sum(READS2(:,1:(SEGMENTS(j,1)-1)),2)==0),(sum(READS2(:,SEGMENTS(j,1):(SEGMENTS(j,2)-1)),2)>0)));
      READS_PER_EXON(1,j)=N1;
      READS_PER_EXON(2,j)=N2;	
    end
  else
    for j=1:size(SEGMENTS,1)
      j;
      N1=sum(sum(READS1(:,SEGMENTS(j,1):(SEGMENTS(j,2)-1))));
      N2=sum(sum(READS2(:,SEGMENTS(j,1):(SEGMENTS(j,2)-1))));
      READS_PER_EXON(1,j)=N1;
      READS_PER_EXON(2,j)=N2;
    end  
  end
  
  EXS_SEQ=EXONSEQUENCE(:,IDX);
  [NEWCOLS,IDX2,POS]=unique(EXS_SEQ','rows');
  NEWCOLS=NEWCOLS';
  
  COMB_READS_PER_EXON=[];
  
  for i=1:max(POS)
    TEMP_IDX=(POS==i);
    if (sum(TEMP_IDX)>0)
      TEMP=sum(READS_PER_EXON(:,TEMP_IDX),2);
      if sum(TEMP)>0
	COMB_READS_PER_EXON=[COMB_READS_PER_EXON,TEMP];
      end
    end
  end
  READS_PER_EXON;
  COMB_READS_PER_EXON;
  
  %renormalization
  MIN_P=1;
  
  READS=sum(COMB_READS_PER_EXON,1);
  NONZERO=[];
  if sum(READS)>0
    NONZERO=READS>0;
    READS=READS(:,NONZERO);
    MINREADS=min(COMB_READS_PER_EXON(:,NONZERO),[],1);
    COMB_READS_PER_EXON(:,NONZERO);
    
    for i=1:size(READS,2)
      MEAN=0.5*READS(i);
      VARIANCE=(READS(i)*0.25).^0.5;      
      Y=sum(((MEAN-MINREADS(i))./(VARIANCE)).^2).^0.5;
      
      %Calculate the p-value
      P_VALUE=1-gammainc((Y^2)/2,1/2);
      if P_VALUE<MIN_P
	MIN_P=P_VALUE;
      end   
    end
    
    MIN_P=MIN_P*size(READS,2);
    P_VALUE=MIN_P;
  else  
    P_VALUE=1;
  end
  
  STRUCT=cell(1,4);
  
  STRUCT{1}=size(READS,2);
  STRUCT{2}=POS;
  STRUCT{3}=NONZERO;
  STRUCT{4}=IDX;
else
  STRUCT=[];
  P_VALUE=1;
end
else
  STRUCT=[];
  P_VALUE=1;
end
