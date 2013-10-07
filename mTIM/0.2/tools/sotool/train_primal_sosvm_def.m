function [information] = train_primal_sosvm_def(PAR, model, trainset)
% Trains a structured output svm using the original yet modified 
% active set method. This version is able to handle both: margin
% rescaling and slack rescaling approaches. 
%
% PAR -- a struct to configure the SO-SVM
% returns a struct recording the training progress
%
% based on code by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

% init optimization problem
[PAR, model, Q, U,d] = init_op(PAR, model);

% init the bmrm values
params = init_bmrm_par(PAR);

% start iterative training
% iteration counter
iter = 1;

% record elapsed time
t_start = clock();

% slack variables for each training instance
nTrain = trainset.num_examples;
slacks = zeros(1,nTrain);

% dimensionality of the problem
DIM = size(Q,1);
fprintf('Problem dimensionality: %i\n',DIM);

% The whole function is build to get this
% optimal weight vector.
% Start with zero weight vector. Otherwise
% random initialization could be possible.
w = zeros(DIM,1);

% subgradients and corresponding offsets
information = {};

if (~params.useMonotConstr),
  fprintf('>NOT< using monotonicity constraints.\n');
  d = [];
  U = [];
else
  % can not handle monotonicity constraints by now
  d = [];
  U = [];
  warning('Not able to handle monotonicity constraints.');
end

normalizer = get_scaling_factor(trainset);
fprintf('Set normalizer to: %2.2f\n',normalizer);
warning('The active set method does not use the normalizer by now.');


% loss function and corresponding subgradient script
% loss function and corresponding subgradient script
loss_fct = params.loss_fct;
loss_sg_fct = params.loss_sg_fct;
lossParams = params.loss
params.loss_fct = loss_fct;
fprintf('Choose loss-fct(%s) and corresponding sg-fct(%s).\n',func2str(loss_fct),func2str(loss_sg_fct));
warning('The active set method does not allow for other loss functions than hinge loss.');

% init Multi-Task Learning 
mtl_sg_linear = zeros(DIM,1);
if (params.mtl_enable),
  mtl_b = params.mtl_b;
  mtl_wstar = params.mtl_wstar;
  % assert a convex combination therefore 0<=b<=1
  assert(mtl_b>=0.0 && mtl_b<=1.0);
  % and assert same dimensionality
  assert(length(mtl_wstar)==DIM);
  
  Q = Q*mtl_b + 2*PAR.C_small*(1-mtl_b)*eye(size(Q,1));
  mtl_sg_linear = -2*PAR.C_small*(1-mtl_b)*mtl_wstar;
  
  w = mtl_wstar;
  fprintf('Multi-Task learning enabled.\nSet new starting point as wstar.\n');
end

% some init values
realObjective = inf;
isConverged=0;
progressTooSmall=0;
trainingTimeStart = clock();

% init objective function values
obj = inf;

% constraints Ax <= b
A = [];
b = [];

% init optimization constants 
res = [w; slacks'];

% lower bound and upper bound of res
lb = [-inf(DIM,1); zeros(nTrain,1)];
ub = inf(length(res),1);

% linear objective constant

% mtl_sg_linear is initialized with zeros in case mtl is disabled
% otherwise it contains the linear part of the mtl regularizer
f = [mtl_sg_linear; ones(nTrain,1)];

% quadratic objective constant
Qs = zeros(length(res));
Qs(1:DIM,1:DIM) = Q;  


% Main loop
while (iter<=PAR.max_num_iter && ~isConverged && ~progressTooSmall),
  fprintf('\n\nActiveSet Iteration %i (%s):\n', iter, datestr(now,'yyyy-mm-dd_HH:MM'));
  new_constraints = zeros(1,PAR.num_train_exm);

  for i = 1:trainset.num_examples,
    % Get the current training example and the corresponding label.
    % Remember that the types of Xi,Yi can be arbitrary 
    % (vector,struct,cell,..).
    [Xi, Yi] = get_training_example(i, trainset);

    % calculate maximum violating example
    Xi.idx = i;
    Xi.iter = iter;
    [PAR, model, w_p, w_n, loss, trn_acc(i)] = argmax_sosvm(PAR, model, w, Xi, Yi);

    weight_delta = w_p - w_n;   
    assert(length(weight_delta) == PAR.num_param);
    if (norm(weight_delta)==0), assert(loss < PAR.epsilon); end

    % is slack rescaling enabled?
    slack_mul = 1;
    if (isfield(PAR,'slack_rescaling') && PAR.slack_rescaling),
        slack_mul = loss;
        fprintf('slack_mul: %f\n',slack_mul);
    end

    % add constraints for examples which have not been decoded correctly
    % and for which a max-margin violator has been found
    currentSlack = loss - slack_mul*w'*weight_delta';
    if ((currentSlack-slacks(i))>PAR.epsilon),
      % new max slack found
      new_constraints(i) = 1;
      v = zeros(1,nTrain);
      v(i) = 1;

      A = [A; -slack_mul*weight_delta, -v];
      b = [b; -loss];
    end
  end
  fprintf('\n');
  fprintf('Generated %i new constraints\n', sum(new_constraints));
  fprintf('Mean training accuracy (prior to solving): %2.1f%%\n', 100*mean(trn_acc));

  % solve intermediate optimization problem
  res = [w; slacks'];
  c_diff = b - A*res;
  part_idx = find(c_diff <= PAR.constraint_margin);
 warning('Infinit loops without progress possible due to heuristics.' );
%  part_idx = 1:length(b);
  fprintf('Solving problem with %2.1f%% of constraints\n\n',100*length(part_idx)/length(b));
  
  % solve intermediate optimization problem
  [res, fobj, exitflag] = quadprog(Qs, f, sparse(A(part_idx,:)), b(part_idx), [], [], lb, ub);
  obj = 0.5*res'*Qs*res + f'*res;
  fprintf('Objectives(exitflag) = %1.4f/%1.4f (%i)\n',obj,fobj,exitflag);
  
  w = res(1:DIM);
  slacks = res(DIM+1:end)';
                                                 
  % warning if minimizer is greater than real objective value
  % (that would be weird)
  diff = 0;
  if (diff < -PAR.epsilon),
    warning('Decrease in objective function %f by %f', obj, diff);
    keyboard;
  end

  % check prediction accuracy on training and holdout examples;
  accs = [];
  if (PAR.check_acc), accs=check_accuracy(PAR, model, trainset, w); end
  
  % collect bundle method debug information and store it
  information{iter}.newConstr = sum(new_constraints);
  information{iter}.trainAcc = trn_acc;
  information{iter}.accs = accs;
  information{iter}.objReal = obj;
  information{iter}.diff = diff;
  information{iter}.slacks = slacks;
  information{iter}.w = w;
  information{iter}.model = model;
  
	% save current solution
  if (iter<7 || mod(iter,5)==0),
    fname = sprintf('sosvm_iter%i',iter);
		fprintf('Saving intermediate result ("%s")...\n\n\n',fname);
		save('-v7',[PAR.out_dir fname], 'PAR', 'information','params');
  end
  
  if (sum(new_constraints)==0), isConverged=1; end;

  % check for too small progress
 % steps = 3;
 % if (iter>10 && iter>steps && abs(information{iter-steps}.objReal-obj)<params.eps),
 %   warning('Terminate due to too small progress.\n');
 %   progressTooSmall=1;
 % end
 % 
  iter = iter + 1;
end



% saving final result
fname = 'sosvm_final';
% name final result file depending wether the training has converged 
% or the number of iterations has exceeded
if (~isConverged), fname = sprintf('sosvm_unconverged_final_iter%i',iter); end
fprintf('Saving final result to "%s".\n\n\n', fname);
save('-v7',[PAR.out_dir fname], 'PAR', 'information','params','lossParams','isConverged','');

% print the overall training time
trainTime = etime(clock(),trainingTimeStart);
fprintf('\nComplete training took % iterations and %8.2fsec (~%ih)\n',iter,trainTime,round(trainTime/3600));

information{end}.PAR = PAR;
information{end}.params = params;
information{end}.isConverged = isConverged;
% eof
