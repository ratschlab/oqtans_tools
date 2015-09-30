function [w, thetas, inds] = train_primal_svm(PAR, lambda, X, y)
% Trains a primal svm using bundle methods.
%
% PAR -- a struct to configure the SVM
% returns a struct recording the training progress
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

% init the bmrm values
params = init_bmrm_par(PAR);


% start iterative training
% iteration counter
iter = 1;

% record elapsed time
t_start = clock();


% slack variables for each training instance
nTrain = length(y);
slacks = ones(1,nTrain);

% dimensionality of the problem
DIM = size(X,1);
fprintf('Problem dimensionality: %i\n',DIM);

% The whole function is build to get this
% optimal weight vector.
% Start with zero weight vector. Otherwise
% random initialization could be possible.
w = zeros(DIM,1);
w_lb = w;

% subgradients and corresponding offsets
ai = [];
bi = [];
wi = []; % corresponding positions
information = {};

% no monotonicity constraints here
d = [];
U = [];


normalizer = 1.0;
fprintf('Set normalizer to: %2.2f\n',normalizer);

% calculate inverse Q-matrix and normalize
Q = lambda * eye(DIM);
Qinv = inv(Q/normalizer);

isConverged=0;
trainingTimeStart = clock();

% loss function and corresponding subgradient script
loss_fct = PAR.loss.fct;
loss_sg_fct = PAR.loss.sg_fct;
fprintf('Choose loss-fct(%s) and corresponding sg-fct(%s).\n',func2str(loss_fct),func2str(loss_sg_fct));

lossParams = PAR.loss


% Main loop
losses = ones(1,nTrain);
dpsis = repmat(y,DIM,1).*X;
dpsisi = 1:nTrain;

params.loss_fct = loss_fct;

fprintf('\n\nBMRM(#bundle=%i) Iteration %i (%s):\n', params.K, iter, datestr(now,'yyyy-mm-dd_HHhMM'));
lossParams.num_examples = nTrain;
[slacks slacksIdxs lossParams] = loss_fct(lossParams, w, losses, dpsis, dpsisi); 

while (iter<=PAR.max_num_iter && ~isConverged),
    
  % calculate correct slack value for this max margin violator
  [a,b] = loss_sg_fct(lossParams, w, losses, dpsis, dpsisi, slacks, slacksIdxs);

  if (norm(a)>1e-08),
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
  else
    % there is no valid subgradient because
    % all slacks are zeros -> break.
    fprintf('Separable data set. Could not determine valid subgradient.\n');
    isConverged = 1;
    continue;
  end
  
  [params, w] = argmin_bmrm(params, Qinv, ai, bi, [],[]);
  w_lb = w;

  params.lossParams = lossParams;
  [params, w] = ls_constrained_parabola(params, Q, losses,w,dpsis,dpsisi,[]);
  lossParams = params.lossParams;

  % calculate the subgradient and offset at the current position
  lossParams.num_examples = nTrain;
  [slacks,slacksIdxs,lossParams] = loss_fct(lossParams, w, losses, dpsis, dpsisi);
  
  % calculate objective function values
  %objBmrmTilt = 0.5*w_lb'*Q*w_lb + normalizer*max(max(ai*w_lb+bi'), params.a_tilt'*w_lb + params.b_tilt);
  objBmrmTilt = 0.5*w_lb'*Q*w_lb + normalizer*(params.a_tilt'*w_lb + params.b_tilt);
  obj = 0.5*w'*Q*w + sum(slacks);
    
  % adapt the linesearch parameter
  diff = obj - objBmrmTilt;
  relGap = diff/obj;

  fprintf('%3i: objectives(%1.2f/%1.2f) gap(%1.4f)\n',iter, obj,objBmrmTilt,relGap);

  % warning if minimizer is greater than real objective value
  % (that would be weird)
  if (diff < -PAR.epsilon),
    warning('Decrease in objective function %f by %f', obj, diff);
    keyboard;
  end
  
  % save and terminate training the gap between the real
  % objective and it's lower bound (bmrmTilt) is smaller
  % or equal to params.eps (like 1e-2).
  if (relGap<params.eps), isConverged=1; end;
  iter = iter + 1;
end



thetas = ones(nTrain,1);
inds = 1:nTrain;

% calculate the weightings for each data point
if (strcmpi(func2str(loss_fct),'loss_mll')),
  num_losses = length(lossParams.mll_losses);
  thetas = ones(nTrain,num_losses);
  
  pstar = lossParams.mll_p / (lossParams.mll_p-1);
  
  
  lossParams.num_examples = nTrain;
  val = 1/(pstar-1);
  for m=1:num_losses,
    [slacks,foo1,foo2] = lossParams.mll_losses{m}(lossParams, w, losses, dpsis, dpsisi);
    thetas(:,m) = slacks.^val;
  end
  
  if (isinf(pstar)), z = sum(thetas,2);
  else z = sum(thetas.^pstar,2).^(1/pstar); end
  z = repmat(z,1,num_losses);

  thetas = thetas./z;

  % check p norm
  
  min(min(thetas))
  max(max(thetas))
  any(abs(sum(thetas.^pstar,2).^(1/pstar)-1)>eps)
  
  [vals, inds] = sort(y.*(w'*X));
end


% print the overall training time
trainTime = etime(clock(),trainingTimeStart);
fprintf('\nComplete training took % iterations and %8.2fsec (~%ih)\n',iter,trainTime,round(trainTime/3600));

% eof
