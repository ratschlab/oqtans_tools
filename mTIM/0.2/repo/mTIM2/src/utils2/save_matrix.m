function save_matrix(filen, M, varname, split_size)
%save_matrix(filen, M, varname, split_size)
%split matrix into parts with size of about split_size 
%bytes (1e+8 default) and save them in one file

if nargin<4
  split_size=1e+8;
end

info = whos('M');
num_splits = ceil(info.bytes/split_size);
num_cols = ceil(size(M, 2)/num_splits);

  
  meta.size=size(M);
  meta.num_splits = num_splits;
  mata.name = varname;
  if fexist(filen)
    save_append(filen, 1, 'meta', meta);
  else
    save(filen, 'meta');
  end
  for j=1:num_splits-1
    %eval(sprintf('%s_split_%i = M(:, 1:num_cols);',varname, j));
    var_split__ = M(:, 1:num_cols);
    M(:, 1:num_cols)=[];
    %save(filen, '-append', sprintf('%s_split_%i',varname,j));
    save_append(filen, 1, sprintf('%s_split_%i',varname,j), var_split__);
    clear var_split__
  end
  if isempty(j), j=0; end;
  var_split__ = M; 
  save_append(filen, 1, sprintf('%s_split_%i',varname,j+1), var_split__);
