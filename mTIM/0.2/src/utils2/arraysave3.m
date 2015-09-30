function arraysave3(outfile_name, array)
% arraysave(outfile_name, array[, arrayname])

assert(size(array,1)==1 || size(array,2)==1)

cf = cellfile(outfile_name, 1) ;
cf = makeempty(cf) ;

fprintf('creating cellfile of size %i\n', length(array)) ;
if iscell(array)
  for i=1:length(array)
    cf{i} = array{i} ;
    if mod(i,100)==0,
      fprintf('%i \r', i) ;
    end ;
  end ;
elseif isstruct(array)
  for i=1:length(array)
    cf{i} = array(i) ;
    if mod(i,100)==0,
      fprintf('%i \r', i) ;
    end ;
  end ;
else
  error('dont know how to handle type') ;
end ;
