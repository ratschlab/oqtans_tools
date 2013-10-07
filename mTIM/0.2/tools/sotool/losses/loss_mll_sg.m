function [a,b] = loss_mll_sg(params, w, delta_y_ybar, delta_psis, ...
  delta_psis_idxs, losses, losses_idxs)
% Calculates the subgradient 'a' and corresponding offset 'b'
% at the current solution 'w' for the p-norm mutliple loss.
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

t = params.num_examples;

% special loss function call
[losses,losses_idxs,params] = loss_mll(params, w, delta_y_ybar, delta_psis, delta_psis_idxs);

% subgradient is calculated by:
%
% a = [(sum_m (mlloss_m^p))^(1/p-1)] .* [sum_m { (mlloss_m)^(p-1) mll_losses_sg_m(w,i)]
%


mlloss = params.mll_ret_loss;
mlidxs = params.mll_ret_idxs;

sgs = params.mll_losses_sg;
p = params.mll_p;

num_losses = length(params.mll_losses);

dim = length(w);
av = zeros(dim,t);

part1 = zeros(1,t);

for m=1:num_losses,
  % calculate individual subgradient and offset (offset is not needed)
  it = length(mlloss{m});
  am = zeros(dim,it);
  for i=1:it,
    params.num_examples = 1;
    [am(:,i) foo] = sgs{m}(params, w, delta_y_ybar, delta_psis, delta_psis_idxs, mlloss{m}(i), mlidxs{m}(i));
  end
    
  av(:,1:it) = av(:,1:it) + repmat((mlloss{m}.^(p-1)),dim,1).*am;

  
  part1 = part1 + mlloss{m}.^p;

end


part1(part1>0.0) = part1(part1>0.0).^(1/p-1);


a = sum(repmat(part1,dim,1) .* av,2);

if (any(isnan(a))),
keyboard;
end
  
  
if (norm(a)>1e3),
%  a = a./(norm(a)*1000);
end


% return subgradient offset
b = sum(losses) - w'*a;  % Remp = a'w + b => Remp-a'w = b
