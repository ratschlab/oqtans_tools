function [PAR, label, signal, state_label, exm_id_intervals] = ...
  load_data_rnd(PAR, dirname, train_inds, vald_inds)
% Loads prepared data
%   [ 25% train_data | 25% vald_data | 50% test_data ]
%
% ARGs
%		PAR         : general parameter structure containing at least
%                 num_train_exm/num_vald_exm/..
%   train_inds  : (optional) if train_inds is present so has vald_inds.
%   vald_inds     These are the indices of the data to use. If not present
%                 then the data is shuffled.
%
% adapted by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011


% partition data for cross-validation
data_file = [dirname '/data.mat'];
if (exist(data_file,'file')),
  fprintf('Loading data file %s...',data_file);
  load(data_file);
  fprintf('Done!.\n');
else
  error('Could not find data file %s.\n',data_file);
end

% get example ids
all_exm_ids = unique(exm_id);

% last half of the data is used for testing purposes
parts = ceil(length(all_exm_ids)/2);
% the remaining half is split into two parts: first for testing
% and the second for validation
% ASSUMES that there are enough examples available for validation and
% training
vald_start = ceil(parts/2);
test_start = parts;


fprintf('Converting label to state sequences..');
state_model = make_model(PAR);
state_label = nan(size(label));
for i=1:length(all_exm_ids),
  idx = exm_id_intervals(i,2):exm_id_intervals(i,3);
  true_label_seq = label(idx);
  obs_seq = signal(:,idx);
  [true_state_seq, problems] = labels_to_states(true_label_seq, state_model, obs_seq, PAR);
  assert(isempty(problems));
  state_label(idx) = true_state_seq;
end
fprintf('Done!\n');

% for train and validation 
exm_id = all_exm_ids(1:parts-1);

% check size
assert(length(exm_id)>=(PAR.num_train_exm+PAR.num_vald_exm));
assert(vald_start>PAR.num_train_exm);
assert((test_start-vald_start)>=PAR.num_train_exm);

% assign training / validation and test examples
PAR.train_exms = [exm_id(1:vald_start-1)];
PAR.vald_exms  = [exm_id(vald_start:test_start-1)];
PAR.test_exms  = [all_exm_ids(test_start:end)];

% use pre-defined indices for train and vald if defined otherwise
% shuffle train/vald
if (exist('train_inds','var') && exist('vald_inds','var')),
  % don't shuffle data instead use train_inds/vald_inds
  PAR.train_exms = train_inds;
  PAR.vald_exms  = vald_inds;  
else
  % shuffle data
  PAR.train_exms = PAR.train_exms(randperm(length(PAR.train_exms)));
  PAR.vald_exms = PAR.vald_exms(randperm(length(PAR.vald_exms)));

  % cut to fit the num_train_exm/num_val_exm 
  % assign training and validation data
  PAR.train_exms = PAR.train_exms(1:PAR.num_train_exm);
  PAR.vald_exms  = PAR.vald_exms(1:PAR.num_vald_exm);
end

fprintf('There are %i test examples and %i available in total.\n',...
  length(PAR.test_exms),length(all_exm_ids));
fprintf('Currently %i are used for training and %i for validation.\n',...
  length(PAR.train_exms),length(PAR.vald_exms));
fprintf('For training and validation %i can be used.\n',...
  length(all_exm_ids)-length(PAR.test_exms));

assert(isempty(intersect(PAR.train_exms, PAR.vald_exms)));
assert(isempty(intersect(PAR.test_exms, [PAR.train_exms, PAR.vald_exms])));
assert(length(PAR.train_exms) >= PAR.num_train_exm);

% eof

