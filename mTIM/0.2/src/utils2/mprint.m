function str=mprint(M) ;
% str=mprint(M) ;

str='' ;
[m,n]=size(M) ;
str=[str '['] ;
for i=1:m,
  str=[str '['] ;
  for j=1:n,
    if ceil(M(i,j))==M(i,j),
      str=[str sprintf('%i', M(i,j))] ;
    elseif abs(M(i,j))<1e-1,
      str=[str sprintf('%1.1e', M(i,j))] ;
    else
      str=[str sprintf('%1.1f', M(i,j))] ;
    end ;
    if j~=n
      str=[str ', '] ;
    end ;
  end ;
  str=[str ']'] ;
  if i~=m,
    str=[str ', '] ;
  end ;
end ;

str=[str ']'] ;