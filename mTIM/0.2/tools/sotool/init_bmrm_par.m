function params = init_bmrm_par(PAR)

% BMRM PARAMETERS
% ---------------
% In order to evaluate the duality gap between the
% bundle method solution and the real value at the
% current solution a loop over all training examples
% and corresponding viterbi decoding at the current w is needed. 
% This is time-consuming but not neccessary after every iteration.
% This parameter sets the number of rounds where the gap
% is not re-evaluated: 
%
% (mod(iter, params.checkRealObjMod)==0) => re-evaluate
% default: 5
% min: 2 
params.checkRealObjMod=4;

% It is possible to either use the monotonicity constraints
% or not. The usage of monotonicity constraints slows down
% optimization significantly (at least factor 30 by rule of thumb).
% On the other hand, the monotonicity constraints restrict the
% solution space such that much less iterations are needed.
%
% default: 1 (use monotonicity constraints)
params.useMonotConstr = 0;

% The training finishes if the relative duality gap drops
% below this value:
%
% ((upperBound-lowerBound)/upperBound)<params.eps => training finished
% default: 1e-1;
params.eps = 1e-1;
%params.eps = 1e-4;

% Aggregation method aka memory-limited method makes
% the bundle methods applicable: instead of a all-time
% increasing bundle use only the params.K last subgradients.
%
% default: 1 (use aggregation)
% default: 118 (arbitrarily >2)
params.useAggregation = 1;
params.K = 78;

% Parameter for the USE_HEURISTIC_MONOT_SET_SELECTION:
% monotonicity shrinkage parameter 
% only use constraints where d-U*w < params.monotMargin
params.monotMargin = 0.01;

% linesearch bounds
params.useLs = 1;
params.lsEps = 1e-3;
params.lsSteps = 50;
params.lsTheta = 0.25;

% load parameters if possible
if (isfield(PAR,'bmrm'))
  fprintf('\nBMRM: Using new parameters:\n');
  params = PAR.bmrm;
  disp(params);
end

% Multi-Task learning parameters. 
% (default): disabled.
params.mtl_enable = 0; 

% parameter for convex combination for the additional
% 2-norm distance regularizer:
% mtl_b * w'Qw + (1-mtl_b) * ||w-mtl_wstar||^2
params.mtl_b = 0.1;
params.mtl_wstar = [];

% load parameters if possible
if (isfield(PAR,'mtl'))
  fprintf('\nMTL: Setting new parameters and w.\n');
  params.mtl_enable = PAR.mtl.mtl_enable;
  params.mtl_b = PAR.mtl.mtl_b;
  params.mtl_wstar = PAR.mtl.mtl_wstar;
  disp(PAR.mtl);
end

% HACK: currently it is neccessary to de-activate the
% line search if multi-task learning is enabled.
if (params.mtl_enable),
  warning('MTL is enabled therefore turn off linesearch.');
  params.useLs = 0;
end


% (Multiple) Losses
% loss function and corresponding subgradient script
params.loss_fct = @loss_hinge_native;
params.loss_sg_fct = @loss_hinge_sg;
params.loss = [];

if (isfield(PAR,'loss_fct')),
  fprintf('\nLOSS: Setting new loss function and parameters.\n');
  params.loss_fct = PAR.loss.fct;
  params.loss_sg_fct = PAR.loss.sg_fct;
  params = PAR.loss;
end


% init values
params.lsObj = inf;
params.lsW = -1;

params.a_tilt = [];
params.b_tilt = [];



