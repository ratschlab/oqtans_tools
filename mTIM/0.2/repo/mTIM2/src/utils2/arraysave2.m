function arraysave2(outfile_name, array, arrayname)
% arraysave2(outfile_name, array[, arrayname])
% 
% obsolete
  
if nargin<3
  arrayname='array' ;
end ;

warning('this function is obsolete and will be removed soon, please use "save_struct" instead') ;

save_struct(outfile_name, array, arrayname) ;
