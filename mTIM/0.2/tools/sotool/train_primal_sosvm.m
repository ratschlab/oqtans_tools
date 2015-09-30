function [information] = train_primal_sosvm(PAR, model, trainset)
% Trains a structured output svm using bundle methods.
%
% PAR -- a struct to configure the SO-SVM
% returns a struct recording the training progress
%
% based on code by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

USE_PRESOLVER = 0;
USE_HEURISTIC_TRAIN_EXM_SELECTION = 1;
USE_HEURISTIC_MONOT_SET_SELECTION = 1;
fprintf('USE_HEURISTIC_TRAIN_EXM_SELECTION: %i\n',USE_HEURISTIC_TRAIN_EXM_SELECTION);
fprintf('USE_HEURISTIC_MONOT_SET_SELECTION: %i\n',USE_HEURISTIC_MONOT_SET_SELECTION);

% init optimization problem
[PAR, model, Q,U,d] = init_op(PAR, model);

% init the bmrm values
params = init_bmrm_par(PAR);

% start iterative training
% iteration counter
iter = 1;

% record elapsed time
t_start = clock();


% Set of active constraints (indices)
dpsis = [];
dpsisi = [];
useCount = [];

% number of argmax-calls (in total)
num_argmax = 0;

% slack variables for each training instance
nTrain = trainset.num_examples;
losse = [];
slacks = zeros(1,nTrain);
slacksIdxs = zeros(1,nTrain);

% dimensionality of the problem
DIM = size(Q,1);
fprintf('Problem dimensionality: %i\n',DIM);

% The whole function is build to get this
% optimal weight vector.
% Start with zero weight vector. Otherwise
% random initialization could be possible.
w = zeros(DIM,1);
w_lb = w;

cneg = cell(1,nTrain);
closs = cell(1,nTrain);
cpos = cell(1,nTrain);

% loss function and corresponding subgradient script
loss_fct = params.loss_fct;
loss_sg_fct = params.loss_sg_fct;
lossParams = params.loss
params.loss_fct = loss_fct;
fprintf('Choose loss-fct(%s) and corresponding sg-fct(%s).\n',func2str(loss_fct),func2str(loss_sg_fct));

% find a good initialization
%SKIP PRE-SOLVER
if (USE_PRESOLVER),
    fprintf('Starting pre-solver...\n');
    f = zeros(DIM,1);
    fprintf(' Gathering feature-vectors.\n' );
    for i=1:nTrain,
      % Get the current training example and the corresponding label.
      % Remember that the types of Xi,Yi can be arbitrary 
      % (vector,struct,cell,..).
      [Xi, Yi] = get_training_example(i, trainset);

      % calculate maximum violating example
      Xi.idx = i;
      Xi.iter = 0;
      [fooPAR, foomodel, w_p, w_n, loss, footrn_acc(i)] = argmax_sosvm(PAR, model, w, Xi, Yi);
      num_argmax = num_argmax+1;

      weight_delta = w_p - w_n;   
      f = f + weight_delta';

      % calculate correct slack value for this max margin violator
      lossParams.num_examples = 1;
      [currentSlack foo_idxs fooParams] = loss_fct(lossParams, w, loss, weight_delta', 1); 
    end
    fprintf(' Solving qp.\n');
    % solve by increasing regularization until the solution lies in the
    % interior of the problem
    mul = 1;
    for i = 1:10,
      fprintf(' Pre-solve iteration %i\n',i);
      K2 = 10;
      w = quadprog(mul*Q,-f,[eye(DIM);-eye(DIM);U],[K2*ones(DIM,1);K2*ones(DIM,1);d]);
      if (sum(abs(w))<(0.7*DIM)), break;
      else mul = mul*10; end
    end
    w_init = w;
    w_lb = w;
    fprintf('Done (max/min=%1.2f,%1.2f).\n ',max(w),min(w));
    %keyboard
else
    fprintf('Not using pre-solver for hot-start.\n');
end

% subgradients and corresponding offsets
ai = [];
bi = [];
wi = []; % corresponding positions
information = {};

if (~params.useMonotConstr),
  fprintf('>NOT< using monotonicity constraints.\n');
  d = [];
  U = [];
  USE_HEURISTIC_MONOT_SET_SELECTION = 0;
end

% amount of monotonicity constraints
num_monot = length(d);

% heuristic:
% process all training examples where the corresponding
% value is not 0.
process_exm = 1:trainset.num_examples;

checkRealObj = 0;

normalizer = get_scaling_factor(trainset);
fprintf('Set normalizer to: %2.2f\n',normalizer);
fprintf('Dimensionality: %i\n',DIM);

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

% calculate inverse Q-matrix and normalize
Qinv = inv(Q/normalizer);

% some init values
realObjective = inf;
isConverged=0;
progressTooSmall=0;
isZeroSubgradient=0;
trainingTimeStart = clock();

% init objective function values
objBmrmTilt = -inf; 
objBmrm = -inf;
obj = inf;
objEst = 0;

% Main loop
PAR.save_now = 1;

w_store = [];



while (iter<=PAR.max_num_iter && ~isConverged && ~progressTooSmall),
  PAR.save_now_iter = iter;    
  PAR.iter = iter;
    
  fprintf('\n\nBMRM(#bundle=%i) Iteration %i (%s):\n', params.K, iter, datestr(now,'yymmdd_HHMM'));
  new_constraints = zeros(1,PAR.num_train_exm);

  % either perform bundle step or check the 
  % function value at FIXED w.
  if (checkRealObj), checkRealObj=0; end;
  if (sum(process_exm)==0 || iter==1 || mod(iter,params.checkRealObjMod)==0),
    % re-active all training examples
    process_exm = 1:trainset.num_examples;
    % set checkRealObj-flag
    checkRealObj = 1;
  end  

  % profile the methods used
  profile on;
  
  % some times for the optimization
  smoTime = 0; solverTime = 0; buildTime = 0;

  num_used_monot = 0;
  isZeroSubgradient = 0;

  % enter gap closing mode if estimated
  % objective and lower bound are too far away.
  gapIter = 0;
  isGapClosingMode = 0;
  isBigEstTiltGap = ((objEst-objBmrmTilt)/objEst)>params.eps;
  if (~checkRealObj && iter>3 && isBigEstTiltGap),
    fprintf('Close gap between estimated obj and lower bound.\n');
    isGapClosingMode = 1;
  end
  
  idx = 1;
  exms = process_exm(process_exm>0);
  trn_acc = -ones(1,trainset.num_examples);
  while ((idx<=length(exms) || isGapClosingMode) && gapIter<50),
    % in gap closing mode don't call argmax function
    % only do the bundle method step
    if (isGapClosingMode),
      gapIter = gapIter+1;
      objBmrmTilt = 0.5*w_lb'*Q*w_lb + (params.a_tilt'*w_lb + params.b_tilt)*normalizer; 
      objEst = objReg + sum(slacks) + w'*mtl_sg_linear;
      if (((objEst-objBmrmTilt)/objEst)<=params.eps),
        fprintf('Is close enough after %i iterations.\n',gapIter);
        isGapClosingMode=0;
      end
    else
      % choose example
      % for each training example search for the maximum violator 
      % and do a bundle step
      i = exms(idx);
      idx = idx+1;

      % Get the current training example and the corresponding label.
      % Remember that the types of Xi,Yi can be arbitrary 
      % (vector,struct,cell,..).
      [Xi, Yi] = get_training_example(i, trainset);

      % calculate maximum violating example
      Xi.idx = i;
      Xi.iter = iter;
      [PAR, model, w_p, w_n, loss, trn_acc(i)] = argmax_sosvm(PAR, model, w, Xi, Yi);
      num_argmax = num_argmax+1;

      % w_n is row-vector
      if iter==1,
         cpos{i} = w_p;
      end
      cneg{i} = [cneg{i}; w_n];
      closs{i} = [closs{i}; loss];
 

      weight_delta = w_p - w_n;   
      assert(length(weight_delta) == PAR.num_param);

      % calculate correct slack value for this max margin violator
      lossParams.num_examples = 1;
      [currentSlack foo_idxs fooParams] = loss_fct(lossParams, w, loss, weight_delta', 1); 

      % add constraints for examples which have not been decoded correctly
      % and for which a max-margin violator has been found
      if ((currentSlack-slacks(i))>0),
        % new max slack found
        new_constraints(i) = 1;

        % w_n is row-vector
        cneg{i} = [cneg{i}; w_n];
        closs{i} = [closs{i}; loss];


        % update
        losse = [losse, loss];
        dpsis = [dpsis, weight_delta'];
        dpsisi = [dpsisi, i];
        useCount = [useCount, iter];

        slacks(i) = currentSlack;
        slacksIdxs(i) = length(losse);
      else
        
        if (~checkRealObj && USE_HEURISTIC_TRAIN_EXM_SELECTION),
            process_exm(i) = 0;
        end

      end
    end

    % IF checkRealObj==1 then 
    % only check the function value at fixed w
    % ELSE
    % generate new subgradient/offset and calculate the new optimal w
    if (checkRealObj),
      staticRealSlacks(i) = currentSlack;
    else
      % ATTENTION: an exotic error occurs if
      % there is no sequence till now that violates the margin
      % (e.g. isempty(losse)==1).
      % Simple but not best solution is: to go directly (without
      % optimizing) to the next sequence.
      % TODO: check for the unlikly case that we start with an
      % init value of w that is already perfect..
      if (isempty(losse)),
        fprintf('Current example and all predecessors dont violate margin and therefore no optimization is possible/neccessary.\n');
        continue;
      end

      % calculate the subgradient and offset at the current position
      lossParams.num_examples = trainset.num_examples;
      [a,b] = loss_sg_fct(lossParams, w, losse, dpsis, dpsisi, slacks, slacksIdxs);
      % if multi-task learning is enabled ->
      % update the subgradient and re-calculate the bias
      if (params.mtl_enable),
        a = a + mtl_sg_linear;
        b = (sum(slacks) + w'*mtl_sg_linear) - w'*a;
      end

      % curious situation: no subgradient means
      % that one optimal solution is achieved (but only on the
      % approx empirical risk)
      if (norm(a)<eps),
        isZeroSubgradient=1;
        fprintf('Zero subgradient.\n');
        if (isGapClosingMode), 
          isGapClosingMode=0; 
        end;
      else
        % re-calculate offset of subgradient at w
        % Remp = a'w + b => Remp-a'w = b
        if (params.useAggregation && length(bi)>params.K),
          ai = [ai(2:end,:); a'./normalizer];
          bi = [bi(2:end), b/normalizer];        
          wi = [wi(2:end,:); w'];
        else
          ai(end+1,:) = a'./normalizer;
          bi(end+1) = b/normalizer;
          wi(end+1,:) = w';
        end

        % monotonicity working set selection
        part_idx=1:num_monot;
        if (USE_HEURISTIC_MONOT_SET_SELECTION), 
          d_diff = d-U*w;
          part_idx = find(d_diff < params.monotMargin);
        end;
        num_used_monot = num_used_monot + length(part_idx);
        if (params.useMonotConstr)
            [params, w] = argmin_bmrm(params, Qinv, ai, bi, U(part_idx,:), d(part_idx));
        else
            [params, w] = argmin_bmrm(params, Qinv, ai, bi, [],[]);
        end
        w_lb = w;

        smoTime = smoTime + params.smoTime;
        solverTime = solverTime + params.solverTime;
        buildTime = buildTime + params.buildTime;

        if (params.useLs),
          params.loss_fct = loss_fct;
          lossParams.num_examples = trainset.num_examples;
          params.lossParams = lossParams;
          if (params.useMonotConstr)
            [params, w] = ls_constrained_parabola(params, Q, losse,w,dpsis,dpsisi,U(part_idx,:));
          else
            [params, w] = ls_constrained_parabola(params, Q, losse,w,dpsis,dpsisi,[]);
          end
          lossParams = params.lossParams;
        end

        % calculate the subgradient and offset at the current position
        lossParams.num_examples = trainset.num_examples;
        [slacks,slacksIdxs,lossParams] = loss_fct(lossParams, w, losse, dpsis, dpsisi);
        useCount(slacksIdxs(slacksIdxs>0)) = useCount(slacksIdxs(slacksIdxs>0))+1;
      end
    end
  end
  profile off;
  
  stats = profile('info');
  [foo inds] = sort([stats.FunctionTable(:).TotalTime]);
  % show ranked profiling results
  for i=inds,
    method = stats.FunctionTable(i);
    name = method.FunctionName;
    time = method.TotalTime;
    %fprintf('%2.2fsec spend in %s\n',time,name);
  end
  
  fprintf('\n');
  fprintf('#monotonicity constraints: %i.\n#average of used monot-constraints: %i\n',...
		num_monot, round(num_used_monot/trainset.num_examples));

  fprintf('#argmax-calls so far: %i\n',num_argmax);
  fprintf('#training examples used: %i\n',length(process_exm(process_exm>0)));
	fprintf('Size of dpsis-matrix (#/bytes): %i/%4.1fMB\n',size(dpsis,2),(size(dpsis,1)*size(dpsis,2)*8)/(1024*1024));

  fprintf('Maximum value in dpsis is %1.2f.\n',max(max(abs(dpsis))));
  fprintf('Maximum value in delta_y_ybar is %1.2f.\n',max(abs(losse)));
  
  fprintf('Generated %i new constraints\n', sum(new_constraints));
  fprintf('Mean training accuracy (prior to solving): %2.1f%%\n', 100*mean(trn_acc(trn_acc>=0.0)));

  % calculate objective function values
  objReg = 0.5*w'*Q*w;
  if (~isempty(params.a_tilt)), 
    objBmrmTilt = 0.5*w_lb'*Q*w_lb + (params.a_tilt'*w_lb + params.b_tilt)*normalizer; 
    objBmrm = objReg + max(max(ai*w+bi'), (params.a_tilt'*w_lb + params.b_tilt))*normalizer;
  end;

  % the linear multi-task part is zero if not used
  objEst = objReg + sum(slacks);
  obj = objReg + sum(staticRealSlacks);

  if (params.mtl_enable),
    objEst = objEst + w'*mtl_sg_linear;
    obj = obj + w'*mtl_sg_linear;
  end
    
  % adapt the linesearch parameter
  diff = abs(obj - objBmrmTilt);
  relGap = abs(diff/obj);
 
  txtUpdate = '';
  if (checkRealObj), 
	txtUpdate = ' (*REAL VALUE UPDATED*) ';
    realObjective = obj;
  end
  obj = realObjective;

  fprintf('\nStats%s:\n',txtUpdate);
  fprintf(' -objectives(real,est,ls,lb) = %1.2f/%1.2f/%1.2f/%1.2f\n',obj,objEst,objBmrm,objBmrmTilt);
  fprintf(' -gap(abs[rel]) = %1.2f[%1.4f]\n', diff, relGap);
  fprintf(' -difference w/o linesearch: %1.2f\n',norm(w-w_lb));    
  fprintf(' -impact of quad. reg.: %1.4f\n',abs(objReg/objBmrmTilt));

  % warning if minimizer is greater than real objective value
  % (that would be weird)
  if (diff < -PAR.epsilon),
    warning('Decrease in objective function %f by %f', obj, diff);
    keyboard;
  end

  % check prediction accuracy on training and holdout examples;
  accs = [];
  dbg = [];
  if (PAR.check_acc && checkRealObj), 
      [accs, dbg] = check_accuracy(PAR, model, trainset, w); 
      fprintf('Approximate convex loss L(train)=%f\n',sum(slacks));
  end
  
  % collect bundle method debug information and store it
  information{iter}.trainTime = etime(clock(),trainingTimeStart);
  information{iter}.newConstr = sum(new_constraints);
  information{iter}.smoTime = smoTime;
  information{iter}.solverTime = solverTime;
  information{iter}.buildTime = buildTime;
  
  information{iter}.trainAcc_pe = trn_acc;
  information{iter}.trainAcc = mean(trn_acc);
  information{iter}.accs = accs;
  information{iter}.dbg = dbg;
  
  information{iter}.checkRealObj = checkRealObj;
  information{iter}.staticRealSlacks = staticRealSlacks;
  
  information{iter}.objReal = obj;
  information{iter}.objEst = objEst;
  information{iter}.objBmrm = objBmrm;
  information{iter}.objBmrmTilt = objBmrmTilt;
  
  information{iter}.diff = diff;
  information{iter}.diffRel = (diff/obj);
  
  information{iter}.w = w;
  information{iter}.w_lb = w_lb;
  information{iter}.model = model;
  
	% save current solution
  if (iter<7 || mod(iter,5)==0),
    fname = sprintf('sosvm_iter%i',iter);
		fprintf('Saving intermediate result ("%s")...\n\n\n',fname);
		save([PAR.out_dir fname], 'PAR', 'information','params');
  end

  % ATTENTION: an exotic error occurs if
  % there is no sequence till now that violates the margin
  % (e.g. isempty(losse)==1).
  % Optimization finished..
  if (isempty(losse)),
    warning('Exotic: we started with an already perfect solution vector w.\n');
    isConverged = 1;
  end

  % save and terminate training the gap between the real
  % objective and it's lower bound (bmrmTilt) is smaller
  % or equal to params.eps (like 1e-2).
  if (checkRealObj && abs(diff/obj)<params.eps && ~isinf(diff)), isConverged=1; end;

  % check for too small progress
  steps = params.checkRealObjMod*3+1;
  if (checkRealObj && iter>steps && abs(information{iter-steps}.objReal-obj)<params.eps),
    warning('Terminate due to too small progress.\n');
    progressTooSmall=1;
  end
  
  iter = iter + 1;
end

% test final training and holdout accuracies
fprintf('Checking final accuracies..\n');
[accs, dbg] = check_accuracy(PAR, model, trainset, w); 
information{end}.accs = accs;
information{end}.dbg = dbg;
fprintf('Approximate convex loss L(train)=%f\n',sum(slacks));
fprintf('Done!\n');

information{end}.PAR = PAR;
information{end}.params = params;
information{end}.isConverged = isConverged;

% saving final result
fname = 'sosvm_final';
% name final result file depending wether the training has converged 
% or the number of iterations has exceeded
if (~isConverged), 
    PAR.unconverged = 1;
    PAR.iteration = iter;
    warning('Not converged.');
end
fprintf('Saving final result to "%s".\n\n\n', fname);
save([PAR.out_dir fname], 'PAR', 'information','params','lossParams','isConverged','cpos','cneg','closs','Q');

% print the overall training time
trainTime = etime(clock(),trainingTimeStart);
fprintf('\nComplete training took % iterations and %8.2fsec (~%ih)\n',iter,trainTime,round(trainTime/3600));


% eof
