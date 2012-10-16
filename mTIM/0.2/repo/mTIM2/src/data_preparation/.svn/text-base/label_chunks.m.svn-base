function [label id_intervals] = label_chunks(chunks, label_map)

% [label id_intervals] = label_chunks(chunks, label_map)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

L = sum(chunks(:,3)-chunks(:,2)+1);
LABELS = get_label_set_mTIM();
undef = intmin('int8');
label = undef * ones(1, L, 'int8');
fn = fieldnames(LABELS);
for f=1:length(fn)
  assert(getfield(LABELS, fn{f}) ~= undef);
end
% these will contain chunk id and corresponding
% start and end positions in label
id_intervals = zeros(size(chunks,1), 3);

b = 1;
e = 1;
for c=1:size(chunks,1),
  idx = chunks(c,2):chunks(c,3);
  e = b + length(idx) - 1;
  label(b:e) = label_map{chunks(c,1)}(idx);
  id_intervals(c,:) = [chunks(c,4), b, e];
  b = e + 1;
end
assert(length(label) == L);
assert(e == L);
assert(~any(label == undef));
assert(~any(any(id_intervals == 0)));

% eof