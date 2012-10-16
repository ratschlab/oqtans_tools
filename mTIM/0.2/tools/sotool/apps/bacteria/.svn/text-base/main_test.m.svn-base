function main_test
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008
% adapted by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

% add mosek path
% addpath('/usr/local/mosek/5/toolbox/r2007a');

addpath('model');

addpath('../hmsvm');
addpath('../hmsvm/native');

addpath('../../src');
addpath('../../src/linesearch');
addpath('../../src/losses');
addpath('../../src/solver/prloqo');

% number of subsets used for cross-validation
num_xval_subsets = 3;

model = 'bacteria';
organism = 'Escherichia_coli_BW2952_uid59391';
%organism = 'Enterobacter_638_uid58727';

% names of parameters to be independently specified (possibly differently
% between training runs) 
param_names = {'C_small', ...
               'num_train_exm', ...
               'train_subsets', ...
              };

% parameter combinations to be used for independent training
parameters = { ...
    [1], [10], [1]; ...
             };
assert(size(parameters,2) == length(param_names));

% basic data directory
dr_base = ['../../tmp/' model '_' organism '/result_' datestr(now,'yyyy-mm-dd_HHhMM')]

% seed for random number generation
rand('seed', 11081979);

% partition data for cross-validation
data_file = ['../../data/' organism '/data.mat'];
load(data_file, 'exm_id');
exm_id = unique(exm_id);
exm_id = exm_id(randperm(length(exm_id)));
subset_ends = round(linspace(0,length(exm_id),num_xval_subsets+1));
for i=1:num_xval_subsets,
  exm_subsets{i} = sort(exm_id(subset_ends(i)+1:subset_ends(i+1)));
end
assert(isequal(sort([exm_subsets{:}]), sort(exm_id)));

for i=1:size(parameters,1),
  PAR = [];
  % constant parameters
  PAR.out_dir = [dr_base '/model' num2str(i) '/']; % output directory
  PAR.model_name = [model '_' organism];           % name of the learning model
  PAR.model_dir = ['models/' PAR.model_name '/'];  % model directory
  PAR.data_file = data_file;                       % name of the training data file
  PAR.num_plif_nodes = 64;                         % number of supporting points
                                                   % for each scoring function
  PAR.constraint_margin = 10;                      % use heuristic training procedure
  PAR.submit_vald = 0;
  PAR.include_paths = {};
  PAR.optimizer = 'prloqo';
  PAR.classifier = @train_primal_sosvm;

  PAR.C_coupling = 0.0;
  PAR.C_smooth = 0.0

% TODO: Nico, das geht leider nicht mehr...
%  PAR.extra_checks = 1;

  PAR.perm_feature_ranges = set_permitted_feature_ranges();
  
  % parameters which vary across HM-SVM training runs
  fprintf('Training model %i...\n', i);
  for j=1:length(param_names),
    if length(parameters{i,j}) == 1,
      fprintf('  %s = %f\n', param_names{j}, parameters{i,j});
    else
      p_str = [];
      for k=1:length(parameters{i,j}),
        p_str = [p_str sprintf('%f ', parameters{i,j}(k))];
      end
      fprintf('  %s = %s\n', param_names{j}, p_str);
    end
    PAR = setfield(PAR, param_names{j}, parameters{i,j});
  end
  fprintf('\n\n');

  % assign example sequences to training, validation and test set
  assert(isfield(PAR, 'train_subsets'));
  PAR.vald_subsets = PAR.train_subsets(end)+1;
  if PAR.vald_subsets>num_xval_subsets,
    PAR.vald_subsets = mod(PAR.vald_subsets,num_xval_subsets);
  end
  assert(all(ismember(PAR.train_subsets, [1:num_xval_subsets])));
  assert(all(ismember(PAR.vald_subsets,  [1:num_xval_subsets])));
  PAR.test_subsets = setdiff(1:num_xval_subsets, ...
                             [PAR.train_subsets PAR.vald_subsets]);
  PAR.train_exms = [exm_subsets{PAR.train_subsets}];
  PAR.vald_exms  = [exm_subsets{PAR.vald_subsets}];
  PAR.test_exms  = [exm_subsets{PAR.test_subsets}];
  assert(isempty(intersect(PAR.train_exms, PAR.vald_exms)));
  assert(isempty(intersect(PAR.test_exms, [PAR.train_exms, PAR.vald_exms])));
  assert(length(PAR.train_exms) >= PAR.num_train_exm);
  
  disp(PAR)
  time = clock();
  train_hmsvm(PAR);
  fprintf('Complete training took %3.2f sec\n',etime(clock(), time));
end

% eof

