function [TP FP TN FN thresholds] = eval_separation(y, label, max_points)

% [TP FP TN FN thresholds] = eval_separation(y, label, max_points)
%
% Calculates the number of true positives (TP), false positives (FP),
% true negatives (TN) and false negatives (FN) given a (continuous)
% prediction score and a label.
%
% y -- (continuous) prediction scores, each row is treated as result of a
%   different prediction method (y can hence have dimension k x n, where
%   k is the number of prediction methods and n the number of data
%   points)
% label -- corresponding labels of size 1 x n (only position where label
%   is either -1 or +1 are considered)
% max_points -- the maximum number of points where TP, FP, TN and FN are
%   assessed 
% returns true positives (TP), false positives (FP), true negatives (TN)
%   and false negatives (FN) each is an k x m matrix where k is equal to
%   the number of rows in y and m is the number of reasonable thresholds
%   (at most max_points if specified). Also returns the thresholds
%   between y's used to trade-off TP and FN (or TN and FP respectively).
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2009

num_curves = size(y,1);

assert(size(y,2) == size(label,2))

if exist('max_points', 'var'),
  % interpolated curves
  num_points = max_points;
else
  % exact curves
  num_points = size(y,2) + 2;
end

TP = zeros(num_curves, num_points);
FN = zeros(num_curves, num_points);
FP = zeros(num_curves, num_points);
TN = zeros(num_curves, num_points);

%tic
for i=1:num_curves,
  assert(~any(isnan(y(i,:))));
  y_min = min(y(i,:));
  y_max = max(y(i,:));
  y_len = size(y,2);
  if exist('max_points', 'var'),
    % interpolated curves
    % try to select thresholds such that they cover a broad range of
    % sensitivity values in an approximately uniform manner
    epsilon = 10^-4;
    t    = nan(1, num_points);
    tp   = nan(1, num_points);
    fn   = nan(1, num_points);
    sn   = nan(1, num_points);
    n    = nan(1, num_points);
    t(1) = y_min-epsilon*(y_min+y_max);
    t(2) = y_max+epsilon*(y_min+y_max);
    t(3) = (y_min+y_max) / 2;
    for j=1:3,
      is_positive = y(i,:) >= t(j);
      tp(j) = sum(label == +1 & is_positive);
      fn(j) = sum(label == +1 & ~is_positive);
      sn(j) = tp(j)/(tp(j) + fn(j));
      n(j) = sum(is_positive);
    end
    for j=4:num_points,
      [t perm] = sort(t);
      sn = sn(perm);
      % keep track of the number of points between consecutive threshold
      % values in order to avoid attempts to split large sensitivity
      % intervals although there are not enough points to further resolve it
      n = n(perm);
      % sort sensitivity intervals in descending order as a candidate
      % list for further resolution
      [tmp p] = sort(sn(1:j-2) - sn(2:j-1), 'descend');
      for k=1:length(p),
        idx = p(k);
        t_try = (t(idx) + t(idx+1)) / 2;
        is_positive = y(i,:) >= t_try;
        n_try = sum(is_positive);
        % take this interval if there are enough data points for further
        % resolution
        if n(idx) - n_try > epsilon*y_len ...
              && n_try - n(idx+1) > epsilon*y_len,
          break
        end
        % otherwise take the next largest interval
      end
      
      t(j) = (t(idx) + t(idx+1)) / 2;
      is_positive = y(i,:) >= t(j);
      n(j) = sum(is_positive);
      tp(j) = sum(label == +1 & is_positive);
      fn(j) = sum(label == +1 & ~is_positive);
      sn(j) = tp(j)/(tp(j) + fn(j));
    end
    [t perm] = sort(t);
    thresholds(i,:) = t;

    clear t tp fn sn n p idx perm tmp t_try n_try is_positive
  else
    % exact curves
    thresholds(i,:) = [y_min-1 sort(y(i,:)) y_max+1];
  end
end
%toc

%tic
for i=1:num_curves,
    for j=1:num_points,
      is_positive = y(i,:) >= thresholds(i,j);
      TP(i,j) = sum(label == +1 & is_positive);
      FN(i,j) = sum(label == +1 & ~is_positive);
      FP(i,j) = sum(label == -1 & is_positive);
      TN(i,j) = sum(label == -1 & ~is_positive);
    end
  assert(all(TP(i,:)+FN(i,:) == TP(i,1)+FN(i,1)));
  assert(all(TN(i,:)+FP(i,:) == TN(i,1)+FP(i,1)));
end
%toc

% eof