function blocks = find_blocks(v)

% blocks = find_blocks(v)
%
% Finds blocks of consecutive '1's.
%
% v -- a 0/1 logical vector
% returns blocks of consecutive ones as [block_starts; block_ends]
%
% written by Georg Zeller & Regina Bohnert, MPI Tuebingen, Germany 2007-2008

% trivial cases
if isempty(v),
  blocks = zeros(0,2);
  return
end
assert(size(v,1)==1 || size(v,2)==1);
if length(v)==1 
  if v==1,
    blocks = [1; 1];
  else
    blocks = zeros(2,0);
  end
  return
end

blocks(1,:) = find(v~=0 & [0 v(1:end-1)]==0);
blocks(2,:) = find(v~=0 & [v(2:end)   0]==0);

% tests for correctness
assert(all(blocks(1,:)<=blocks(2,:)));
check = zeros(size(v));
for i=1:size(blocks,2),
  check(blocks(1,i):blocks(2,i)) = 1;
end
assert(isequal(check~=0, v~=0));

% eof