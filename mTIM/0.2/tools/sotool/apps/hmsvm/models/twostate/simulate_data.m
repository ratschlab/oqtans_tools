function data_file = simulate_data()
% data_file = simulate_data()
%
% Generates simulated data, where features are generated from the label
% sequence with some labels swapped and Gaussian noise added.
%
% returns the name of the file to which data was written
%
% written by Georg Zeller, MPI Tuebingen, Germany

num_exm = 1000;         % number of examples
exm_len = 500;          % length of each example sequence
num_features = 8;       % total number of features
num_noise_features = 4; % number features to be pure noise
block_len = [10, 100];   % min an max lentgh of positive block
num_blocks = [0, 6];    % min and max number of positive block per example
num_subsets = 1;        % number of subsets for crossvalidation

prop_distort = 0.2;     % proportion of wrong labels
noise_std = 4;          % standard deviation of Gaussian noise


exm_id = [];
pos_id = [];
label = [];
subset_id = [];
i = 1;
while (i<=num_exm),
  % generate label sequence randomly
  % containing num_blocks(1) to num_blocks(2) blocks of positive labels
  % each of length between block_len(1) and block_len(2)
  l = -ones(1,exm_len);
  rnb = num_blocks(1) + ceil((num_blocks(2)-num_blocks(1)).*rand(1)) - 1;
  for j=1:rnb,
    rl = block_len(1) + ceil((block_len(2)-block_len(1)).*rand(1)) - 1;
    rp = ceil((exm_len-rl).*rand(1));
    l(rp:rp+rl) = 1;
  end
  
  % sanity check: every end of a example should
  % be a stop state
  if (l(1)~=-1 || l(end)~=-1),
    disp(sprintf('Invalid start/end label (%i)! Retry.',i));
  else 
    exm_id = [exm_id i*ones(1,exm_len)];
    pos_id = [pos_id 1000*i + (1:exm_len)];

    label = [label l];
    rs = ceil((num_subsets).*rand(1));
    subset_id = [subset_id rs*ones(1,exm_len)];
    
    i = i+1;
  end
end

% generate features by i) introducing label noise, i.e. flipping a
% proportion prop_distort of labels and ii) adding gaussian noise to the
% (distorted) label sequence
for i=1:num_features,
  distort = randperm(length(label));
  d1 = distort(1:round(length(label)*prop_distort));
  d2 = distort(end-round(length(label)*prop_distort)+1:end);
  l = label;
  l(d1) = l(d2);
  signal(i,:) = l+noise_std*randn(size(label));
end
% substitute some features by pure noise
ridx = randperm(num_features);
ridx = ridx(1:num_noise_features);
signal(ridx,:) = noise_std*randn(length(ridx), size(label,2));
fprintf('noise features: %i\n', ridx);

data_dir = '../../../../data/';
if ~exist(data_dir, 'dir'),
  mkdir(data_dir);
end
data_file = [data_dir 'twostate_toy.mat'];
save(data_file, 'label', 'signal', 'pos_id', 'exm_id', 'subset_id');
