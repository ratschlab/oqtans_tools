function [a,b] = loss_logistic_sg(params, w, delta_y_ybar, delta_psis, ...
  delta_psis_idxs, losses, losses_idxs)
% Calculates the subgradient 'a' and corresponding offset 'b'
% at the current solution 'w' for the logistic regression loss of given 
% delta_psis and delta_y_y_bar-metric. 
%
% Attention!
% In contrast to the corresponding loss function (here: loss_logistic.m)
% we expect that delta_* only contains active constraints (therefore,
% #delta_psis<=t).
%
% Subgradient a and offset b:
% Remp = sum(losses)
% Remp = a'*w + b
% b = Remp - a'*w
% a = dRemp/dw
% 
% n : #all generated max-margin violators
% d : #dimensions
% t : #training examples (or less in the first iterations)
%
% dPsis_i(ybar) = psi(x_i,y_i) - psi(x_i,ybar)
%
%
% IN
%   params          : parameter structure (not needed here)
%   w               : [d x 1] current solution vector
%   delta_y_ybar    : [1 x n] 'loss'-metric for all generated maximum
%                     margin violators of all training examples
%   delta_psis      : [d x n] corresponds to psi(x,y)-psi(x,ybar)
%   delta_psis_idxs : [1 x n] corresponding delta_psis and delta_y_ybar
%                     belong to this training example
%   losses          : [1 x t] empirical risk term value at w
%   losses_idxs     : [1 x t] empirical risk term value at w
%
% OUT
%   a               : [d x 1] subgradient at w
%   b               : (scalar) corresponding offset
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

cutoff = 5;
alpha = 1;
if (isfield(params,'logistic_alpha')), alpha = params.logistic_alpha; end
if (isfield(params,'logistic_cutoff')), cutoff = params.logistic_cutoff; end

% treshold
t = log(1+exp(cutoff));
% linear slope
s = exp(cutoff)/(1+exp(cutoff));

% calculate subgradient and offset
idxs = losses_idxs(losses_idxs>0);
losses = losses(losses_idxs>0);

% find linear and logistic losses
inds = find(losses<cutoff); % (1..t)
inds2 = idxs(inds); % logistic indices (1..n)
inds1 = setxor(idxs,inds2); % linear indices (1..n)

% calculate empirical error
Remp = sum(losses);

% Derivative is (for the logistic part):
%
% t = alpha*(delta_y_ybar - w'*delta_psis)
% -alpha*delta_psis * exp(t)/(1+exp(t))
%
% .. and for the linear part: -s*alpha*delta_psis;
dim = length(w);

mul = exp(losses(inds));
mul = (mul-1)./mul;
mul = repmat(mul,dim,1);

a = zeros(dim,1);

% linear losses
if (sum(inds1>0)==1), a=-alpha*s*delta_psis(:,inds1); end
if (sum(inds1>0)> 1), a=-alpha*s*sum(delta_psis(:,inds1),2); end

% logistic losses (extend a)
if (sum(inds2>0)==1), a=a - alpha*delta_psis(:,inds2).*mul; end
if (sum(inds2>0)> 1), a=a - alpha*sum(delta_psis(:,inds2).*mul,2); end

% subgradient bias
b = Remp - w'*a;  % Remp = a'w + b => Remp-a'w = b