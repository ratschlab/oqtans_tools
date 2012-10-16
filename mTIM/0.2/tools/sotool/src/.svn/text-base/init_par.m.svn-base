function PAR = init_par(PAR)

% PAR = set_default_par(PAR)
%
% Augments and checks a PAR (parameter) struct.
%
% PAR -- a struct to configure the HM-SVM (for specification see
%   setup_hmsvm_training.m)
% returns a PAR struct augmented with fields set to default values 
% (see comments below)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2011


hmsvm_home = which('init_par');
if hmsvm_home(end) == '/',
  hmsvm_home(end+1) = '/';
end

if ~isfield(PAR, 'include_paths'),
  PAR.include_paths = {};
end


% option to enable/disable some extra consistency checks
if ~isfield(PAR, 'extra_checks'),
  PAR.extra_checks = 0;
end

% option to control the amount of output
if ~isfield(PAR, 'verbose'),
  PAR.verbose = 1;
end

% option to enable/disable performance checks during training
if ~isfield(PAR, 'check_acc'),
  PAR.check_acc = 1;
end

% stopping criterion: constraint generation is terminated if no more
% margin violations are found or the relative change of the objective
% function is smaller than this parameter...
if ~isfield(PAR, 'min_rel_obj_change'),
  PAR.min_rel_obj_change = 10^-3;
end
% ... or if the maximum number of iterations is exceeded
if ~isfield(PAR, 'max_num_iter'),
  PAR.max_num_iter = 1000;
end

% margin constraints are only added if the example is predicted with an
% accuracy below this parameter
if ~isfield(PAR, 'max_accuracy'),
  PAR.max_accuracy = 0.99;
end
% and if the max margin violator incurs a loss at least as high as this
% parameter
if ~isfield(PAR, 'min_loss'),
  PAR.min_loss = 1;
end

% numerical tolerance to check consistent score calculation, constraint
% satisfaction and monotonictiy of the objective function
if ~isfield(PAR, 'epsilon'),
  PAR.epsilon = 10^-5;
end

% option to only solve partial intermediate training problems which do
% not contain constraints satisfied with a margin at least as large as
% the parameter value. Such constraints are however kept aside and
% checked in each iteration. Set to inf to always solve the full problem.
% Throwing away constraints is a HEURISTIC which speeds up training at
% the cost of losing the guarantee to converge to the correct solution!
if ~isfield(PAR, 'constraint_margin'),
  PAR.constraint_margin = inf;
end

% optimization software used to solve the (intermediate) training
% problem(s). Currently there are two possibilities: 'cplex' or 'mosek'
if ~isfield(PAR, 'optimizer'),
  PAR.optimizer = 'cplex'
end

% subsample examples for performance checks
if ~isfield(PAR, 'max_num_vald_exms'),
  PAR.max_num_vald_exms = 100;
end

% by default, do not submit any cluster jobs from within HM-SVM training
if ~isfield(PAR, 'submit_jobs'),
  PAR.submit_jobs = 0;
end

% only for plain hmm
if ~isfield(PAR, 'hmm_min_level'),
 PAR.hmm_min_level = 0;
end

% seed for random number generation
% be careful! very nasty bug if not wanted
%rand('seed', 11081979);

% mandatory fields of the parameter struct
assert(isfield(PAR, 'C_small'));
assert(isfield(PAR, 'C_smooth'));
assert(isfield(PAR, 'C_coupling'));
assert(isfield(PAR, 'out_dir'));

if ~exist(PAR.out_dir, 'dir'),
  mkdir(PAR.out_dir);
end

% eof
