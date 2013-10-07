function [loss,idxs,params] = loss_logistic(params, w, delta_y_ybar, delta_psis, delta_psis_idxs)
% Calculates the logistic loss of given delta_psis and delta_y_y_bar-metric. 
% 
% t = alpha*(d-w'dpsi);
% l(t<cutoff) = log(1+exp(t))   // logistic
% l(t>=cutoff) = s*t            // linear (with same derivative in cutoff)
%
% n : #all generated max-margin violators
% d : #dimensions
% t : #training examples
%
% dPsis_i(ybar) = psi(x_i,y_i) - psi(x_i,ybar)
%
% IN
%   params          : parameter structure 
%                     params.num_examples
%                     params.logistic_alpha (default=1)
%                     params.logistic_cutoff (default=5)
%                     (num training examples == t)
%   w               : [d x 1] current solution vector
%   delta_y_ybar    : [1 x n] 'loss'-metric for all generated maximum
%                     margin violators of all training examples
%   delta_psis      : [d x n] corresponds to psi(x,y)-psi(x,ybar)
%   delta_psis_idxs : [1 x n] corresponding delta_psis and delta_y_ybar
%                     belong to this training example
%
% OUT
%   loss            : [1 x t], i=1,2,..,t
%                     loss_i=max_ybar(delta(y_i,ybar) - w'dPsis_i(ybar))
%   idxs            : [1 x t], points to the delta_* that
%                     that belongs to the corresponding loss
%   params          : (updated) parameter structure
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

t = params.num_examples;

cutoff = 5;
alpha = 1;
if (isfield(params,'logistic_alpha')), alpha = params.logistic_alpha; end
if (isfield(params,'logistic_cutoff')), cutoff = params.logistic_cutoff; end

% PART I: calculate losses
x = alpha*(delta_y_ybar - w'*delta_psis);
s = exp(cutoff)/(1+exp(cutoff));

inds1 = find(x>cutoff);
inds2 = setxor(1:length(delta_psis_idxs),inds1);

losses(inds2) = log(1+exp(x(inds2)));
losses(inds1) = (x(inds1)-cutoff)*s + log(1+exp(cutoff));

idxs = zeros(1,t);
loss = zeros(1,t);

for i=1:t,
  inds = find(delta_psis_idxs==i);
  if (~isempty(inds)),
    [loss(i) ind] = max(losses(inds));
    idxs(i) = inds(ind);
  end
end
