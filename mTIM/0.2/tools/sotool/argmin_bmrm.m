function [params, w, p] = argmin_bmrm(params, Qinv, as, bs, C,d)
% Single optimization step with CRBM as proposed by TMT Do 2010 for either
% vanilla crbm and crbm with aggregation. This is a special/adapted version 
% that also takes linear inequalities into account.
%
% IN
%   params : structure
%            params.useAggregation (0/1) 
%            params.K (K>2 if useAggregation==1) : size of aggregation
%                     model (=number of planes used-1 plus aggregation
%                     plane)
%            params.a_tilt : aggr plane gradient vector
%            params.b_tilt : aggr plane bias 
%   Qinv   : invers Q (regularization) matrix 
%   as     : cutting plane gradient vectors [as(num,dim)]
%   bs     : cutting plane biases [bs(num)]
%   C      : plif regularization inequation Cw <= d 
%   d      : ~ (for mTim: expected to be zero)
%
% OUT
%   params : updated parameters
%            params.a_tilt : updated aggregated gradient
%            params.b_tilt : updated aggregated bias
%   w      : new weight vector 
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, 2010

% choose solver
SOLVER_MOSEK  = 0;
SOLVER_SMO    = 0;
SOLVER_PRLOQO = 1;

% methods also delivers solver times
params.smoTime = 0;
params.solverTime = 0;

% indicate errors or inconsistencies 
params.isError = 0;


% number of cutting planes so far
n_elems = length(bs);

% solve optimization problem in dual
if (params.useAggregation && n_elems >= (params.K-1))
  % 1. choose last (agg-1) sg's
  params.inds = (n_elems-params.K+2):n_elems;
  % 2. chooses highes xis
  %[foo, params.inds] = sort(-xis);
  params.inds = params.inds(1:params.K-1);
  as = [as(params.inds,:); params.a_tilt'];
  bs = [bs(params.inds), params.b_tilt];
end

% some indices
asEnd = length(bs);

tb = clock();
H = [as; C]*Qinv'*[as; C]'; % objectives quadratic term
f = [-bs'; zeros(length(d),1)]; % objectives linear term
B = [ones(1,length(bs)), zeros(1,length(d))]; % equality linear term
lb = zeros(length(f),1);
params.buildTime = etime(clock(),tb);

% -- Mosek
%[x,fval,exitflag,output,lambda]=quadprog(H,f,A,b,B,c,l,u,x0,options)
% minimize     0.5*x'*H*x+f'*x    
% subject to         A*x <= b 
%                    B*x  = c
%                 l <= x <= u 
%
top = clock();
if (SOLVER_MOSEK),
  [alphas,fval,exitflag,output,lambda] = quadprog(H, f, [], [], B, 1, lb, []);
  params.solverTime = etime(clock(),top); 

  % check if the optimizer delivered an correct solution
  % ATTENTION!!  sometimes mosek6 returns a zero vector
  if (sqrt(sum(alphas.*alphas)) <= length(f)*eps),
    warning('OPTIMIZER DELIVERED AN  INVALID SOLUTION!');
    keyboard;

    % and try again until some valid solution was found.
    % if there is a solution, then mosek is guilty otherwise
    % it is likely that there is some bug in our formulation
    maxIters = 100;
    iters = 0;
    foundNewSolution = 0;
    fprintf('Calling again...\n');
    while (iters<maxIters && ~foundNewSolution),
      alphasNew = quadprog(H, f, [], [], B, 1, lb, []);

      iters = iters+1;
      if (sqrt(sum(alphas.*alphas))>length(f)*eps), 
        foundNewSolution=1
        fprintf('Round %i/%i delivered a valid solution..\n',iters,maxIters);
      else
        fprintf('Invalid solution in round %i/%i..\n',iters,maxIters);
      end
    end

    keyboard;
    alphas = alphasNew;
  end
end


% -- PrLoqo
%    minimize   c'*x + 1/2 x'*H*x
%    subject to A*x = b
%               l <= x <= u
prloqo_time = clock();
if (SOLVER_PRLOQO),
  UBinf = 1e6; 
  ub = UBinf*ones(length(f),1);
  % shogun pr-loqo
  %[alphas, y] = sg('pr_loqo',f', H, A, b', lb', ub');
  %alphas = alphas';
  % matlab pr-loqo
  %warning off;
  alphas = zeros(length(f),1);
  inds = find(diag(abs(full(H))) >= 1e-16);
  
  if (isempty(inds)),
    fprintf('ERROR! inds is empty.');
    save('-v7','prloqo_debug.mat','inds','alphas','f','H','B','lb','ub','C','as','bs','Qinv','params');
    [alphas, foo] = pr_loqo2(f, H, B, 1, lb, ub);
    params.isError = 1;
  else
    [alphas(inds), foo] = pr_loqo2(f(inds), (H(inds,inds)), B(inds), 1, lb(inds), ub(inds));
  end
  params.solverTime = etime(clock(),prloqo_time);
end

%{
% SMO
% EXPERIMENTAL
% triple smo test
A = as*Qinv'*as';
%B = as*Qinv'*C';
%Cc = C*Qinv'*C';
%[alphas2,betas2] = triple_smo_grad_wss(A,B,Cc,bs,alphas, -Qinv*[as;C]', C,d);
smo_time = clock();
%[alphas2,betas2] = triple_smo_grad_wss(A,B,Cc,bs,alphas, [],[],[]);
alphas2 = zeros(length(f),1);
inds = find(diag(abs(full(A))) >= 1e-18);
[alphas2(inds)] = simple_smo(A(inds,inds),-f(inds),alphas(inds));
d_ftsmo = sum((alphas2-alphas).^2);
if (any(isnan(alphas2)) || any(isinf(alphas2)) || d_ftsmo>1e-3)
  fprintf('Difference (alphas): %f\n ',d_ftsmo);
  keyboard
end
%keyboard
params.smoTime = etime(clock(),smo_time);
%}


% new aggregated hyperplanes parameter
params.b_tilt = bs*alphas(1:asEnd);
params.a_tilt = [as; C]'*alphas;
 
% optimizer of above qp 
w = -Qinv*params.a_tilt;
