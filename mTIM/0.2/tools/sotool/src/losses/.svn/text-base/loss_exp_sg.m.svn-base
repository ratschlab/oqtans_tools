function [a,b] = loss_exp_sg(params, w, delta_y_ybar, delta_psis, ...
  delta_psis_idxs, losses, losses_idxs)
% Calculates the subgradient 'a' and corresponding offset 'b'
% at the current solution 'w' for the exp loss of given delta_psis 
% and delta_y_y_bar-metric. 
%
% Attention!
% In contrast to the corresponding loss function (here: loss_exp.m)
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
if (isfield(params,'exp_alpha')), alpha = params.exp_alpha; end;

% calculate subgradient and offset
idxs = losses_idxs(losses_idxs>0);
losses = losses(losses_idxs>0);
Remp = sum(losses);

% ATTENTION! zero values can cause problems
% if not handled correctly during optimization.
dim = length(w);
mul = repmat(losses,dim,1);

a = zeros(dim,1);

if (sum(idxs>0)==1), a=-alpha*delta_psis(:,idxs).*mul; end
if (sum(idxs>0)> 1), a=-alpha*sum(delta_psis(:,idxs).*mul,2); end
b = Remp - w'*a;  % Remp = a'w + b => Remp-a'w = b