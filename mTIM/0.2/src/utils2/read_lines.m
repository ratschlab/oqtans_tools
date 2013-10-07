function lines=read_lines(fname) ;

fd=fopen(fname,'r+') ;
str=fread(fd) ;
fclose(fd) ;
str=setstr(str)' ;
lines=[] ;
if isempty(str),
  return ;
end ;
id=find(str==10) ;
id=[0 id] ; 
for i=1:length(id)-1 ;
  lines{i}=str(id(i)+1:id(i+1)-1) ;
end ;
