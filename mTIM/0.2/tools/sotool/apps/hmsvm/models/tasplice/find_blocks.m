function blocks = find_blocks(v)
% blocks = find_blocks(v)
% expects a 0/1 vector v as input and returns blocks of consecutive
% '1's as [block_starts; block_ends]

assert(size(v,1)==1);
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