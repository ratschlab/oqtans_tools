function arraysave(outfile_name, array, arrayname)
% arraysave(outfile_name, array[, arrayname])
% 
% obsolete
  
if nargin<3
  arrayname='array' ;
end ;

warning('this function is obsolete and will be removed soon, please use "save_cell" instead') ;

save_cell(outfile_name, array, arrayname) ;

