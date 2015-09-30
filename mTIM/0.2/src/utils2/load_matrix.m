function var=load_matrix(filename, varname)
% var=load_matrix(filename, varname)

load(filename, 'meta')
if ~exist('meta')
  load(filename, varname);
  eval(sprintf('var = %s;',varname))
  return
end

var = [] ;
col = 1;
for j=1:meta.num_splits
  %if col>1
  %  assert(~all(var(:,col-1)==-1* ones(meta.size(1),1)))
  %end
  %assert(all(var(:,col)==-1*ones(meta.size(1),1)))
   
  load(filename, sprintf('%s_split_%i',varname,j));
  eval(sprintf('temp = %s_split_%i;',varname,j));

  if j==1,
    if isnumeric(temp)
      var = zeros(meta.size(1), meta.size(2), class(temp));
    else
      var = zeros(meta.size(1), meta.size(2),'uint8');
      var = char(var);
    end
  end ;

  var(:,col:col + size(temp, 2)-1)=temp;
  col = col + size(temp, 2);
end
%assert(~all(var(:,end)==-1* zeros(meta.size(1),1)))
