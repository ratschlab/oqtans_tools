function [PAR, label, signal, state_label, exm_id_intervals] = ...
  load_multiple_organisms(PAR, base_dir, in_dirs, num_train, num_val, organism)
% Load multiple organisms and glue them in memory together.
%
% Uses 'load_data_rnd.m'.
%
% IN
%   PAR
%   base_dir
%   in_dirs
%   num_train
%   num_val
%   (optional) organism : same structure as returned PAR.organism
%
% OUT
%   PAR
%   label
%   signal
%   state_label
%   exm_id_intervals
%
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

num_orgs = length(in_dirs);

PAR.organism = [];

% return values
label = [];
signal = [];
state_label = [];
exm_id_intervals = [];


train_exm_ids = [];
vald_exm_ids = [];
test_exm_ids = [];

% increment id by this value for each 
% organism -> assume that 
% no organism has more genes than 'id_inc'.
id_inc = 10000;

% for each organism
for i=1:num_orgs,
  % name of the current organism
  PAR.organism{i}.name = in_dirs{i}; 
  
  PAR.num_train_exm = num_train;
  PAR.num_vald_exm = num_val;
  PAR.max_num_vald_exms = PAR.num_vald_exm;
  
  if (exist('organism','var')),
    % find corresponding organism
    train_inds = [];
    vald_inds = [];
    for j=1:length(organism),
      if (strcmpi(organism{j}.name,PAR.organism{i}.name)),
        train_inds = organism{j}.org_train_exms;
        vald_inds = organism{j}.org_vald_exms;
      end
    end
    assert(~isempty(train_inds));
    assert(~isempty(vald_inds));
    
    [n_PAR, n_label, n_signal, n_state_label, n_exm_id_intervals] = ...
      load_data_rnd(PAR, [base_dir '/' in_dirs{i}], train_inds, vald_inds);
  else
    [n_PAR, n_label, n_signal, n_state_label, n_exm_id_intervals] = ...
      load_data_rnd(PAR, [base_dir '/' in_dirs{i}]);
  end

  % glue them together
  offset = length(signal);
    
  % only for sanity check
  check_len = 4;
  test_sequences = zeros(size(n_exm_id_intervals,1),1);
  test_start_sequences = zeros(size(n_exm_id_intervals,1),check_len);
  test_stop_sequences = zeros(size(n_exm_id_intervals,1),check_len);
  
  % set new offsets taking the new position in the
  % overal sequence into account
  for j=1:size(n_exm_id_intervals,1),
    assert(n_exm_id_intervals(j,1)<id_inc);

    % store a small test sequence
    start = n_exm_id_intervals(j,2);
    stop =  n_exm_id_intervals(j,3);
    if ((stop+check_len-1)<=length(n_signal)),
      test_sequences(j) = 1;
      test_start_sequences(j,:) = n_signal(start:start+check_len-1);
      test_stop_sequences(j,:) = n_signal(stop:stop+check_len-1);
    end
    
    % set new id    
    n_exm_id_intervals(j,1) = (i-1)*id_inc + n_exm_id_intervals(j,1);
    % add offset to bounds
    n_exm_id_intervals(j,2) = offset + n_exm_id_intervals(j,2);
    n_exm_id_intervals(j,3) = offset + n_exm_id_intervals(j,3);    
  end

  % collect new training/validation and test example ids
  new_test_exm_ids = n_PAR.test_exms+(i-1)*id_inc;
  new_vald_exm_ids = n_PAR.vald_exms+(i-1)*id_inc;
  new_train_exm_ids = n_PAR.train_exms+(i-1)*id_inc;

  train_exm_ids = [train_exm_ids n_PAR.train_exms+(i-1)*id_inc];
  vald_exm_ids = [vald_exm_ids n_PAR.vald_exms+(i-1)*id_inc];
  test_exm_ids = [test_exm_ids n_PAR.test_exms+(i-1)*id_inc];
  
  % glue all together
  signal = [signal n_signal];
  label = [label n_label];
  state_label = [state_label n_state_label];
  exm_id_intervals = [exm_id_intervals; n_exm_id_intervals];
  
  % test for start and stop sequences in the new environment
  for j=1:size(n_exm_id_intervals,1),
    if (test_sequences(j)),
      % store a small test sequence
      start = n_exm_id_intervals(j,2);
      stop =  n_exm_id_intervals(j,3);
      
      assert(all(test_start_sequences(j,:)==signal(start:start+check_len-1)));
      assert(all(test_stop_sequences(j,:)==signal(stop:stop+check_len-1)));
    end
  end
  
  % store a footprint of the current organism
  PAR.organism{i}.from = offset;
  PAR.organism{i}.to = offset+1 + length(n_signal);
  PAR.organism{i}.test_exms = new_test_exm_ids;
  PAR.organism{i}.vald_exms = new_vald_exm_ids;
  PAR.organism{i}.train_exms = new_train_exm_ids;
  PAR.organism{i}.org_train_exms = n_PAR.train_exms;
  PAR.organism{i}.org_vald_exms = n_PAR.vald_exms;
end

% the final 'super'-organism
PAR.num_train_exm = num_train*num_orgs;
PAR.num_vald_exm = num_val*num_orgs;
PAR.max_num_vald_exms = PAR.num_vald_exm;

PAR.train_exms = train_exm_ids;
PAR.vald_exms = vald_exm_ids;
PAR.test_exms = test_exm_ids;
