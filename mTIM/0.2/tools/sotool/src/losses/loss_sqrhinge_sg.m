function [a,b] = loss_sqrhinge_sg(params, w, delta_y_ybar, delta_psis, ...
  delta_psis_idxs, losses, losses_idxs)
% Calculates the subgradient 'a' and corresponding offset 'b'
% at the current solution 'w' for the hinge loss of given delta_psis 
% and delta_y_y_bar-metric. 
%
% Attention!
% In contrast to the corresponding loss function (here: loss_linear.m)
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

alpha = 0.1;
if (isfield(params,'sqrhinge_alpha')), alpha = params.sqrhinge_alpha; end;

% IMPORTANT: bound the norm of the subgradient otherwise
% the values may become too big
NORM_BOUND = 1e5;

% calculate subgradient and offset
idxs = losses_idxs(losses>0);
losses = losses(losses>0);

Remp = sum(losses);

dim = length(w);
mul = repmat(sqrt(losses),dim,1);

% ATTENTION! zero values can cause problems
% if not handled correctly during optimization.
a = zeros(length(w),1);
if (sum(idxs>0)==1), a=-2*alpha*delta_psis(:,idxs).*mul; end
if (sum(idxs>0)> 1), a=-2*alpha*sum(delta_psis(:,idxs).*mul,2); end

% bound subgradient a if neccessary
norm_a = norm(a);
if (norm(a)>NORM_BOUND),
  %a = a * (NORM_BOUND/norm_a);
end

b = Remp - w'*a;  % Remp = a'w + b => Remp-a'w = b