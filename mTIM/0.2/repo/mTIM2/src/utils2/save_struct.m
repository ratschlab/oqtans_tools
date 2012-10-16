function save_struct(outfile_name, array, arrayname, MAXSIZE)
% save_struct(outfile_name, array[, arrayname])

if nargin<3
  arrayname='array' ;
end ;

if nargin<4, 
  MAXSIZE = 1e9; % about 1Gb
end ;

whosinfo = whos('array');
if whosinfo.bytes > MAXSIZE
  % linearize first, store original size
  dimensions = size(array) ;
  array = reshape(array, 1, prod(dimensions)) ;

  elementsize = whosinfo.bytes/length(array);
  chunklength = round(MAXSIZE/elementsize);
  numvars = ceil(length(array)/chunklength);
  %fprintf('arraysave2 processes %i chunks:', numvars) ;
  for ix = 1:numvars-1
    eval(sprintf('%s%d = array(%ld:%ld);',arrayname, ix, (ix-1)*chunklength+1, ix*chunklength));
    %fprintf('.') ;
  end
  %fprintf('\n') ;
  eval(sprintf('%s%d = array(%ld:end);',arrayname, numvars, (numvars-1)*chunklength+1));
  clear array ix ;
  fprintf('save_struct saves %i chunks:', numvars) ;
  save(outfile_name, '-V7', 'chunklength', 'elementsize', 'numvars', 'dimensions') ;
  for ix = 1:numvars
    save(sprintf('%s_%i.mat', outfile_name, ix), '-V7', sprintf('%s%d', arrayname, ix)) ;
    fprintf('.') ;
  end ;
  fprintf('\n') ;
else
  %eval(sprintf('%s=array;', arrayname)) ;
  save_append(outfile_name, 0, arrayname, array) ;
end
