function array = arrayload2(infile_name, arrayname)
% array = arrayload2(infile_name, arrayname)
% 
% obsolete
  
if nargin<2,
  arrayname='array' ;
end ;

warning('this function is obsolete and will be removed soon, please use "load_struct" instead') ;

array = load_struct(infile_name, arrayname) ;

  
  