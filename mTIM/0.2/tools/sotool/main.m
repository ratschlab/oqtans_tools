function main
% Starts HM-SVM training (calling train_hmsvm.m) for different
% combination of hyperparameters. Can thus be used for cross-validation
% and model selection.
%
% see train_hmsvm.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008
% adapted by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2010

% add mosek path
%addpath('/usr/local/mosek/5/toolbox/r2007a');
addpath('~/Documents/mosek/6/toolbox/r2009b');

addpath('native');
addpath('lai');

addpath('../../src');
addpath('../../src/linesearch');
addpath('../../src/losses');
addpath('../../src/solver/prloqo');

rmpath('models/tasplice');
rmpath('models/twostate');

addpath('models/twostate');
addpath('models/tasplice');

% which model to use
model = 'twostate';
model = 'tasplice';

% number of subsets used for cross-validation
num_xval_subsets = 2;

% names of parameters to be independently specified (possibly differently
% between training runs) 
param_names = {'C_small', ...
               'C_smooth', ...
               'C_coupling', ...
               'num_train_exm', ...
               'reg_type'...
               'train_subsets', ...
              };

% parameter combinations to be used for independent training
parameters = { ...
%    [0.1],  [5], [5], [200], 'QP', [1]; ...
    [0.1],  [0], [1.5], [100], 'QP', [1]; ...
%    [0.1],  [1], [0.1], [100], 'QP', [1]; ...
%    [0.0001],  [0.001], [0.001], [200], 'QP', [1]; ...
%    [0.1] ,  [1], [0.1], [200], 'QP', [3]; ...
%    [0.01],  [0.1], [0.05], [100], 'QP', [1]; ...
             };
assert(size(parameters,2) == length(param_names));

% basic data directory
dr_base = ['../../tmp/' model '/result_' datestr(now,'yymmdd_HHMM')]

% seed for random number generation
rand('seed', 11081979);

% partition data for cross-validation
data_file = ['~/Data/lsltoy/' model '_toy.mat'];

data_file = ['../../data/' model '_toy.mat'];

load(data_file, 'exm_id');
exm_id = unique(exm_id);
exm_id = exm_id(randperm(length(exm_id)));
subset_ends = round(linspace(0,length(exm_id),num_xval_subsets+1));
for i=1:num_xval_subsets,
  exm_subsets{i} = sort(exm_id(subset_ends(i)+1:subset_ends(i+1)));
end
assert(isequal(sort([exm_subsets{:}]), sort(exm_id)));

JOB_INFO = [];
for i=1:size(parameters,1),
  PAR = [];

  %PAR.loss.fct = @loss_logistic_native;
  %PAR.loss.sg_fct = @loss_logistic_sg;
  %PAR.loss.fct = @loss_sqrhinge_native;
  %PAR.loss.sg_fct = @loss_sqrhinge_sg;
  PAR.loss.fct = @loss_hinge_native;
  PAR.loss.sg_fct = @loss_hinge_sg;
  
  % constant parameters
  PAR.out_dir = [dr_base '/model' num2str(i) '/']; % output directory
  PAR.model_name = model;                    % name of the learning model
  PAR.model_dir = ['models/' PAR.model_name '/'];  % model directory
  PAR.data_file = data_file;                       % name of the training data file
  PAR.num_plif_nodes = 5;                          % number of supporting points
                                                   % for each scoring function
  PAR.constraint_margin = 10;                      % use heuristic training procedure
  PAR.submit_vald = 0;
  PAR.include_paths = {};
  PAR.optimizer = 'prloqo';

  % enable slack rescaling?
  PAR.slack_rescaling = 1;
  PAR.classifier = @train_primal_sosvm_def;
%  PAR.classifier = @train_primal_sosvm;
%  PAR.classifier = @train_hmm;

  % enable mtl
  PAR.mtl.mtl_enable = 0;
  PAR.mtl.mtl_b = 0.5;
  PAR.mtl.mtl_wstar = ones(324,1);

%  PAR.num_features = 8;
  
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

