function array = load_struct(infile_name, arrayname)
% array = load_struct(infile_name, arrayname)
  
if nargin<2,
  arrayname='array' ;
end ;

infile_name
L_=load(infile_name);

if isequal(infile_name(end-3:end), '.mat'),
  infile_name = infile_name(1:end-4);
end


if isfield(L_,'numvars')
  fprintf('arrayload2 loads %i chunks:', L_.numvars) ;
  ix=1;
  L=load(sprintf('%s_%i.mat', infile_name, ix),'-mat', sprintf('%s%d', arrayname, ix)) ;
  eval(sprintf('array = L.%s1;', arrayname)) ;
  clear L;
  pos = length(array) ;
  array(prod(L_.dimensions))=array(1) ;
  %if keyboard_allowed(), keyboard ; end ;
  fprintf('.') ;
  for ix = 2:L_.numvars
    L=load(sprintf('%s_%i.mat', infile_name, ix), '-mat', sprintf('%s%d', arrayname, ix)) ;
    %if size(array,1) == 1
    %  array = [array, eval(sprintf('L.%s%d', arrayname, ix))];
    %elseif size(array,2) == 1
    %  array = [array; eval(sprintf('L.%s%d', arrayname, ix))];
    %end
    tmp=eval(sprintf('L.%s%d', arrayname, ix)) ; clear L ;
    assert(pos+length(tmp)<=length(array)) ;
    array(pos+1:pos+length(tmp)) = tmp ;
    pos = pos+length(tmp) ;
    clear tmp 
    fprintf('.') ;
  end
  %if keyboard_allowed(), keyboard ; end ;
  fprintf('\n') ;
  if isfield(L_, 'dimensions'),
    array = reshape(array, L_.dimensions) ;
  end ;
else
  eval(sprintf('array=L_.%s;', arrayname)) ;
end
