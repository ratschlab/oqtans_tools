function array = arrayload(infile_name, arrayname)
% array = arrayload(infile_name, arrayname)
% 
% obsolete
  
if nargin<2,
  arrayname='array' ;
end ;

warning('this function is obsolete and will be removed soon, please use "load_cell" instead') ;

array = load_cell(infile_name, arrayname) ;

  
  