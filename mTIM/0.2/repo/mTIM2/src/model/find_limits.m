function limits = find_limits(num_plif_nodes, feat_values)

% limits = find_limits(num_plif_nodes, feat_values)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2011

feat_values = sort(feat_values);
unq = unique(feat_values);
if length(unq) < num_plif_nodes,
  % for discrete features with only very few possible values
  % leave some limits = 0 and fill the unique values into the right-hand
  % side of the limits vector (only works, if all feature values >= 0) 
  unq(unq==0) = [];
  limits = zeros(1,num_plif_nodes);
  if any(unq ~= 0),
    assert(all(unq>0));
    limits(end-length(unq)+1:end) = unq;
  end
elseif length(unq) < 2*num_plif_nodes
  % for discrete features with relatively few possible values
  % start from the number of unique feature values and eliminate the most
  % rarely seen ones until left with the desired number of limits
  limits = unq;
  n = hist(feat_values, limits);
  while length(limits) > num_plif_nodes;
    [mn mn_idx] = min(n);
    if mn_idx == 1,
      n(mn_idx+1) = n(mn_idx+1) + mn;
      limits(mn_idx) = [];
      n(mn_idx) = [];
    elseif mn_idx == length(limits),
      n(mn_idx-1) = n(mn_idx-1) + mn;
      limits(mn_idx) = [];
      n(mn_idx) = [];      
    else
      n(mn_idx-1) = n(mn_idx-1) + round(mn/2);
      n(mn_idx+1) = n(mn_idx+1) + round(mn/2);
      limits(mn_idx) = [];
      n(mn_idx) = [];
    end
  end
else
  % for features with many possible discrete or continuous values
  % first, reduce the number of 0-valued features
  if ~any(feat_values<0) && mean(feat_values==0)>1/num_plif_nodes, 
    feat_values(feat_values==0) = [];
    feat_values = [zeros(1,round(length(feat_values)/num_plif_nodes)) feat_values];
  end
  % second, place limits such that there are equally many feature values
  % "around" each limit (i.e. limits are centered on equally populated bins)
  N = num_plif_nodes;
  limits = [];
  while length(limits) < num_plif_nodes;
    limits = linspace(1, length(feat_values), N);
    limits = round((limits(1:end-1)+limits(2:end))/2);
    % the following unique statement can result in a reduced number of
    % limits; that's why we increase N in the loop (but in most cases we
    % will just pass through the while loop once) 
    limits = unique(feat_values(limits));
    N = N + 1;
%    if mod(N,num_plif_nodes) == 0,
%      fprintf('%i trials, %i limits found\r', N-num_plif_nodes, length(limits));
%    end
    % give up the hope that this strategy is successful
    if N > 10000*num_plif_nodes,
      fprintf('N=%i\n',N);
      %  keyboard
      break
    end
  end
  while length(limits) > num_plif_nodes,
    % this can (in rare cases) happen due to the above unique statement
    bins = [-inf mean([limits(1:end-1); limits(2:end)]) inf];
    cnt = zeros(1,length(limits));
    for b=1:length(bins)-1,
      cnt(b) = sum(feat_values>=bins(b) & feat_values<bins(b+1));
    end
    [tmp del_idx] = min(cnt);
    limits(del_idx) = [];
  end

%  fprintf('%i trials, %i limits found\n', N-num_plif_nodes, length(limits));
end

assert(length(limits) == num_plif_nodes);
%{
bins = [-inf mean([limits(1:end-1); limits(2:end)]) inf];
cnt = zeros(1,length(limits));
tl = {};
for b=1:length(bins)-1,
  cnt(b) = sum(feat_values>=bins(b) & feat_values<bins(b+1));
  tl{b} = num2str(limits(b));
end
figure
bar(cnt)
set(gca, 'XTickLabel', tl);
%}

% eof
