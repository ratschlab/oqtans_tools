function [idx_1 idx_2 len_1 len_2] = compare_intervals_sorted(intervals_1, intervals_2, mode, t)

% [idx_1 idx_2 len_1 len_2] = compare_intervals(I, J, mode, t)
%
% I is expected to be an nx2 matrix of (positive integer) intervals [I(i,1),
% I(i,2)] for i=1...n. which is sorted such that either I(i,1) < I(i+1,1) or
% I(i,1) = I(i+1,1) and I(i,2) <= I(i+1,2) for all i=1...n.  Similarly, J is
% expected to be an mx2 matrix of sorted (poitive integer) intervals
%  
% Returns two cell arrays idx_1 and idx_2 of length n and m
% respectively where the i-th entry in idx_1
% contains a vector of intervals from J similar/overlapping 
% to interval i (from I) and vice versa for idx_2.
% Optionally also returns for each interval the length of overlap to
% other intervals. Note that length calculations assume closed integer
% intervals, i.e. the overlap length between [1,4] and [2,8] is 
% 4-2+1 = 3 !
%
% Mode must be one of 'similar' (intervals with boundaries that
% differ each by at most t are returned) or 'overlap' (intervals
% that overlap by at least t are returned)
% 
% t is a threshold value for similar boundaries/overlap
%
% Written by Georg Zeller, MPI Tuebingen, Germany, 2008

assert(all(all(intervals_1 >= 0)));
assert(all(all(intervals_2 >= 0)));

switch mode,
 case 'similar', % find intervals with boundaries that differ 
                 % each by at most t
  assert(t>=0);
  idx_1 = cell(1,size(intervals_1,1));
  idx_2 = cell(1,size(intervals_2,1));
  q = size(intervals_2,1);
  for p=1:size(intervals_1,1),
    % look to the left in intervals_2
    while q >= 1 && intervals_2(q,1)+t >= intervals_1(p,1),
      if abs(intervals_1(p,1) - intervals_2(q,1)) <= t ...
        && abs(intervals_1(p,2) - intervals_2(q,2)) <= t,
        idx_1{p}(end+1) = q;
        idx_2{q}(end+1) = p;
      end
      q = q - 1;
    end
    if q==0, 
      q = 1; 
    end
    
    % look to the right in intervals_1
    while q <= size(intervals_2,1) ...
             && intervals_2(q,1)-t <= intervals_1(p,1),
      if abs(intervals_1(p,1) - intervals_2(q,1)) <= t ...
            && abs(intervals_1(p,2) - intervals_2(q,2)) <= t,
        idx_1{p}(end+1) = q;
        idx_2{q}(end+1) = p;
      end
      q = q + 1;
    end
    if q>size(intervals_2,1), 
      q = size(intervals_2,1); 
    end
  end
 case 'overlap' % find intervals that overlap by at least t
  assert(t>=1);
  idx_1 = cell(1,size(intervals_1,1));
  idx_2 = cell(1,size(intervals_2,1));
  
  % find all cases where intervals_1(p,1) <= intervals_2(q,1)
  q = 1;
  for p=1:size(intervals_1,1),
    % init q
    while q <= size(intervals_2,1) ...
             && intervals_2(q,1) <= intervals_1(p,2),
      q = q + 1;
    end
    if q>size(intervals_2,1), 
      q = size(intervals_2,1); 
    end
    
    % look to the left in intervals_2    
    while q >= 1 && intervals_2(q,1) >= intervals_1(p,1),
      if min(intervals_1(p,2), intervals_2(q,2)) ...
            - max(intervals_1(p,1), intervals_2(q,1)) + 1 >= t, 
        idx_1{p}(end+1) = q;
        idx_2{q}(end+1) = p;
      end      
      q = q - 1;      
    end
    if q==0, 
      q = 1; 
    end
  end
  
  % find all cases where intervals_1(p,1) >= intervals_2(q,1)
  p = 1;
  for q=1:size(intervals_2,1),
    % init p
    while p <= size(intervals_1,1) ...
             && intervals_1(p,1) <= intervals_2(q,2),
      p = p + 1;
    end
    if p>size(intervals_1,1), 
      p = size(intervals_1,1); 
    end
    
    % look to the left in intervals_2    
    while p >= 1 && intervals_1(p,1) >= intervals_2(q,1),
      if min(intervals_1(p,2), intervals_2(q,2)) ...
            - max(intervals_1(p,1), intervals_2(q,1)) + 1 >= t, 
        idx_1{p}(end+1) = q;
        idx_2{q}(end+1) = p;
      end      
      p = p - 1;
    end
    if p==0, 
      p = 1; 
    end
  end
end

switch nargout,
  case 1,
   for i=1:length(idx_1),
     idx_1{i} = unique(idx_1{i});
   end
 case 2,
  for i=1:length(idx_1),
    idx_1{i} = unique(idx_1{i});
  end
  for i=1:length(idx_2),
    idx_2{i} = unique(idx_2{i});
  end
 case 3,
  for i=1:length(idx_1),
    idx_1{i} = unique(idx_1{i});
    len_1(i) = 0;
    for j=1:length(idx_1{i}),
      ovl = min(intervals_1(i,2), intervals_2(idx_1{i}(j),2)) ...
            - max(intervals_1(i,1), intervals_2(idx_1{i}(j),1)) + 1;
      len_1(i) = len_1(i) + ovl;
    end
  end
  for i=1:length(idx_2),
    idx_2{i} = unique(idx_2{i});
  end
 case 4,
  for i=1:length(idx_1),
    idx_1{i} = unique(idx_1{i});
    len_1(i) = 0;
    for j=1:length(idx_1{i}),
      ovl = min(intervals_1(i,2), intervals_2(idx_1{i}(j),2)) ...
            - max(intervals_1(i,1), intervals_2(idx_1{i}(j),1)) + 1;
      len_1(i) = len_1(i) + ovl;
    end
  end
  for i=1:length(idx_2),
    idx_2{i} = unique(idx_2{i});
    len_2(i) = 0;
    for j=1:length(idx_2{i}),
      ovl = min(intervals_2(i,2), intervals_1(idx_2{i}(j),2)) ...
            - max(intervals_2(i,1), intervals_1(idx_2{i}(j),1)) + 1;
      len_2(i) = len_2(i) + ovl;
    end
  end
 otherwise,
  error('expected 1 to 4 output arguments');
end
