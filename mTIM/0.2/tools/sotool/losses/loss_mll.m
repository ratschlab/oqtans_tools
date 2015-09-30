function [loss,idxs,params] = loss_mll(params, w, delta_y_ybar, delta_psis, delta_psis_idxs)
% Multiple Loss Learning
% 
% n : #all generated max-margin violators
% d : #dimensions
% t : #training examples
%
% dPsis_i(ybar) = psi(x_i,y_i) - psi(x_i,ybar)
%
% IN
%   params          : parameter structure 
%                     params.num_examples  = num training examples = t
%                     params.mll_losses: cell array with loss_fct handles
%                     params.mll_p : p-norm
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
loss = zeros(1,t);
idxs = 1:t;

params.mll_ret_loss = {};
params.mll_ret_idxs = {};

p = params.mll_p;

losses_fct = params.mll_losses;

num_losses = length(params.mll_losses);
for i=1:num_losses,
  
  [l,inds,params] = losses_fct{i}(params, w, delta_y_ybar, delta_psis, delta_psis_idxs);
  params.mll_ret_loss{i} = l;
  params.mll_ret_idxs{i} = inds;
  
  loss = loss+l.^p;
end
loss = loss.^(1/p);