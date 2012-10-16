function [seqs,names]=read_fasta_file(fname, name_only)
% [seqs,names]=read_fasta_file(fname)
% 
% reads all entries in a fasta file
  
if nargin<2, name_only=0; end ; 

fd=fopen(fname,'r') ;
str=fread(fd) ;

% fix new lines
num13=sum(str==13) ;
num10=sum(str==10) ;
if num10==0 && num13>0,
  str(str==13)=10 ;
end ;
if num10>0 && num13>0,
  assert(num10==num13) ;
  str(str==13)=[] ;
end ;

newseq=find(str=='>') ;
newseq(end+1) = length(str)+1 ;

seqs={} ;
names={} ;

% find those that are not at the beginning of a line
idx_rm = [] ;
for i=1:length(newseq)-1,
  if newseq(i)>1 && str(newseq(i)-1)~=10
    idx_rm(end+1)=i ;
  end ;
end ;
newseq(idx_rm)=[] ;

for i=1:length(newseq)-1,
  seq=char(str(newseq(i):newseq(i+1)-1)) ;
  idx=find(seq==10, 1, 'first') ;
  assert(length(idx)==1) ;
  name = seq(2:idx-1)' ;
  names{end+1}=name ;
  if ~name_only,
    seq  = seq(idx+1:end)' ;
    seq(seq==10) = [] ;
    seqs{end+1} = seq ;
  end ;
end ;

fclose(fd) ;
