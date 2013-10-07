function idx = find4structarrays(A,val,fieldname)
  
% idx = find4structarrays(A,val,fieldname)
% This function searches for val in the structure array A, where
% fieldname specifies the field that is searched. 
%
% INPUT   A           structure array, 
%         val         number or string 
%         fieldname   
% 
% The values in A.(fieldname) can be of class 'cell', 'char', or any number,
% Yet all entries must be of same class.  (e.i. genes(1).name = 'na1';
% genes(2)=[] results in an error.)
 
Len = zeros(size(A));
idx_all = [];
field_type=class(A(1).(fieldname));
for i=1:length(A), 
  a=A(i).(fieldname);
  len(i) = length(a);
  if ~isequal(field_type,'char')
    idx_all = [idx_all repmat(i,1,len(i))] ;
  end
  assert(isequal(field_type,class(a)));
end
assert(sum(len)==length(idx_all));



if isequal(field_type,'char')
  idx = find(strcmp({A.(fieldname)},val));
elseif isequal(field_type,'cell')
  idx = find4struct(A.(fieldname),val);
else
  if any(len>1) | any(len==1)
    idx = find([A.(fieldname)]==val);
    idx = idx_all(idx);
  else
    idx = find([A.(fieldname)]==val);  
  end  
end


function idx = find4struct(varargin); 


A = varargin(1:end-1);
val = varargin{end};
idx = [];

if isequal(class(val),'char')
  for i=1:length(A), 
    if any(strcmp(val,A{i})), 
      idx = [idx,i];
    end, 
  end
else
  for i=1:length(A), 
    if ismember(val,[A{i}{:}]), 
      idx = [idx,i];
    end, 
  end
end
