function arraysave(outfile_name, array, arrayname)
% arraysave(outfile_name, array[, arrayname])

assert(size(array,1)==1 || size(array,2)==1)

if nargin<3
  arrayname='array' ;
end ;

MAXSIZE = 1e9; % about 1Gb

whosinfo = whos('array');
if whosinfo.bytes > MAXSIZE
  elementsize = whosinfo.bytes/length(array);
  chunklength = round(MAXSIZE/elementsize);
  numvars = ceil(length(array)/chunklength);
  for ix = 1:numvars-1
    eval(sprintf('%s%d = array(%d:%d);',arrayname, ix,(ix-1)*chunklength+1, ...
                 ix*chunklength));
  end
  eval(sprintf('%s%d = array(%d:end);',arrayname, numvars,(numvars-1)*chunklength+1));
  clear array ix ;
  save(outfile_name);
else
  eval(sprintf('%s=array;', arrayname)) ;
  save(outfile_name,arrayname);
end
