function [pred_label_seqs pred_state_seqs] = predict_hmsvm(data_file, predictor_dirs)

% [pred_label_seqs pred_state_seqs] = predict_hmsvm(data_fn, predictor_dirs)
%
% Predicts label sequences given a trained HM-SVM model and a data_file
% containing features and example ids.
%
% data_file -- a data file containing a signal matrix of feature values
%   and a vector exm_id containing sequence identifiers (same format as
%   required for train_hmsvm.m)
% predictor_dirs -- a cell array in which each entry specifies a
%   directory where a trained HM-SVM model has been saved; if more than
%   one predictor is specified, a cross-validation test set will be
%   reconstructed by concatenation of test examples from individual
%   predictors
% returns concatenated predicted label sequences (same format as label
%   vector required for train_hmsvm.m), optionally also returns the
%   corresponding predicted state sequences (vector of the same dimension)
%
% see train_hmsvm.m
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008

% adjust set_hmsvm_paths.m to point to the correct directories
init_paths();

% switch to evaluate on test examples (if 1) 
% or on validation examples (if 0)
USE_TEST_DATA = 1

% control of the amount of output
VERBOSE = 0

% maximum number of examples evaluated (just to save time)
% set this to inf, if you want to evaluate all available examples
MAX_NUM_EXM = inf;

% number of iterations used for training the mSTAD model 
% (if empty, the model from the final iteration will be loaded)
iteration = [];

% seed for random number generation
rand('seed', 11081979);


load(data_file, 'signal', 'exm_id');
pred_label_seqs = nan(size(exm_id));
pred_state_seqs = nan(size(exm_id));
exms = unique(exm_id);
num_exms = length(exms);

% identify test examples of each of the predictors
pred_exms = zeros(length(predictor_dirs), num_exms);
for p=1:length(predictor_dirs),
  if predictor_dirs{p}(end) ~= '/',
    predictor_dirs{p}(end+1) = '/';
  end
  if isempty(iteration)
    fn = sprintf('%slsl_final.mat', predictor_dirs{p});
  else
    fn = sprintf('%slsl_iter%i.mat', predictor_dirs{p}, iteration);
  end
  predictors(p) = load(fn, 'PAR', 'score_plifs', 'transition_scores');
  
  if VERBOSE>=1,
    fprintf(out_fid, 'Loaded predictor from %s...\n\n', fn);
    disp(predictors(p).PAR);
  end

  if USE_TEST_DATA,
    pred_exms(p, ismember(exms, predictors(p).PAR.test_exms)) = 1;
  else
    pred_exms(p, ismember(exms, predictors(p).PAR.vald_exms)) = 1;    
  end
end

if length(predictors) > 1,
  % if more than one predictor is given, re-assemble cross-validation test
  % data, by assigning a predictor at random where test examples overlap
  idx = find(sum(pred_exms)>1);
  for j=1:length(idx),
    p = find(pred_exms(:,idx(j)));
    p = p(randperm(length(p)));
    p = p(1);
    single_p = zeros(length(predictors),1);
    single_p(p) = 1;
    pred_exms(:,idx(j)) = single_p;
  end
  assert(all(sum(pred_exms)==1));
end
assert(max(max(pred_exms)==1));

% make predictions
fprintf('Predicting on data from %s\n', data_file);
for p=1:length(predictors),
  if VERBOSE>=2,
    fprintf('  Using predictor %i for chunks ...\n', p);
  end
  transition_scores = predictors(p).transition_scores;
  score_plifs = predictors(p).score_plifs;
  PAR = predictors(p).PAR;
  if exist(PAR.model_name)~=7,
    addpath(PAR.model_dir);
  end
  
  % include user-specified include paths
  if isfield(PAR, 'include_paths'),
    for i=1:length(PAR.include_paths),
      addpath(PAR.include_paths{i});
    end
  end
  
  for i=1:sum(pred_exms(p,:)),
    idx = find(pred_exms(p,:));
    idx = find(exm_id==exms(idx(i)));
    obs_seq = signal(:,idx);
    
    %%% Viterbi decoding
    pred_path = decode_viterbi(obs_seq, transition_scores, score_plifs, PAR);
    pred_label_seqs(idx) = pred_path.label_seq;
    pred_state_seqs(idx) = pred_path.state_seq;
  end
end

