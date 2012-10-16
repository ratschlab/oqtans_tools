function array = load_cell(infile_name, arrayname)
% array = load_cell(infile_name, arrayname)
  
if nargin<2,
  arrayname='array' ;
end ;
 
L=load(infile_name) ;
if isfield(L,'numvars')
  eval(sprintf('array = L.%s1;', arrayname)) ;
  for ix = 2:L.numvars
    if size(array,1) == 1
      array = [array, eval(sprintf('L.%s%d', arrayname, ix))];
    elseif size(array,2) == 1
      array = [array; eval(sprintf('L.%s%d', arrayname, ix))];
    end
  end
else
  eval(sprintf('array=L.%s;', arrayname)) ;
end
