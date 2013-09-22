function [reads1]=clip_reads(reads1,CLIP_NUCL)
% This function clips the first CLIP_NUCL bases from reads1
if or(size(reads1,1)==0,size(reads1,2)==0)
  return
end

for C_ITER=1:CLIP_NUCL
  %Clip the first bases
  [TEMP,FIRST_VECT]=max(reads1,[],2);
  FIRST_MASK=sparse((1:length(FIRST_VECT))',FIRST_VECT,ones(length(FIRST_VECT),1),size(reads1,1),size(reads1,2));
  reads1(FIRST_MASK>0)=0;
  
  %Clip the last bases 
  [TEMP,FIRST_VECT]=max(reads1(:,end:(-1):1),[],2);
  FIRST_VECT=size(reads1,2)-FIRST_VECT+1;
  FIRST_MASK=sparse((1:length(FIRST_VECT))',FIRST_VECT,ones(length(FIRST_VECT),1),size(reads1,1),size(reads1,2));
  reads1(FIRST_MASK>0)=0;
  
end
