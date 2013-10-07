function [normalizer] = get_scaling_factor(trainset)
% Because values are sometimes too small or too big
% a (scalar) scaling factor is sometimes needed during training.
% If not needed at alljust return 1.0.
%
% Multiplying an optimization problem with a constant
% (scalar) does change the function value but not 
% the optimimal solution w.

exm_id_intervals = trainset.exm_id_intervals;
train_exm_ids = trainset.train_exm_ids;

idx = find(ismember(exm_id_intervals(:,1), train_exm_ids));
L = sum(exm_id_intervals(idx,3) - exm_id_intervals(idx,2) + length(idx));
fprintf('\n%i training sequences with a total length of %i\n', length(train_exm_ids), L);

% Normalization is needed because otherwise the values of
% the subgradients and the losses will grow with the length
% of the training examples.
% If the values become too big the optimization problem
% might contain too big values as well (>10^9).
normalizer = 1.0;
for i=1:length(train_exm_ids),
    idx = find(exm_id_intervals(:,1)==train_exm_ids(i));
    len = exm_id_intervals(idx,3)-exm_id_intervals(idx,2) + 1;
    if (normalizer<len), normalizer=len; end
end
% normalizer contains now the size of the longest
% training sequence. It's sufficient that it scales
% with the longest sequence:
normalizer=normalizer/10;
