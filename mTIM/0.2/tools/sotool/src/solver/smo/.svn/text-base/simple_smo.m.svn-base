function [alphas] = simple_smo(A,b,sol)
% Fast Triple SMO solves the following optimization problem:
%
% min  0.5 alpha' A alpha - alpha'b
% s.t. sum_i alpha_i = 1
%      alpha_i >= 0   forall i
%
% KKT (complementary slackness): 
%     alpha_i*(a_i*w + b_i - xi) = 0
% 
% A = as*Qinv'*as'
% w = -Qinv*as'*alphas;
%
% Primal problem formulation
% min  0.5 w'Qw + Xi
% s.t. a_i' w + bi <= Xi      \forall i
%      atilt' w + btilt <= Xi
%
% The algorithm uses maximum gradients for the working set selection:
%
% DO
%   choose alpha_i and alpha_j
%   optimize(alpha_i,alpha_j)
% UNTIL some stopping criterion
%
%
% WARNING: 
% - alpha-version 
%

eps = 1e-6;

num_alphas = size(A,1);

% init with valid values
alphas = ones(num_alphas,1)/num_alphas; 

iter = 0;
obj = inf;
dobj = inf;

if (num_alphas==1),
	eta = 0.01;
	diff = 1;
  iter = 0;
else
  noChanges = 0;

  %while (iter<=1000 && dobj>eps),
  while (iter<=1000 && ~noChanges),
		noChanges = 1;

		for k=1:num_alphas-1,
		  for l=k+1:num_alphas,
        %ainds = 1:num_alphas;
        %ainds = ainds(ainds~=k);
        %rpi = randperm(length(ainds));
        %ainds = ainds(rpi);
      	%l = ainds(1);
        [alphas,hasChanged] = simple_smo_update(k,l, alphas, A,b);            
        if (hasChanged), noChanges=0; end;
			end
    end

    distance = sum((alphas-sol).^2);
    fprintf('%3i: dist(%2.4f)\n',iter,distance);
    
		%oldObj = obj;
		%obj = 0.5*alphas'*A*alphas - alphas'*b';
		%dobj = oldObj-obj;
		% increase loop counter 
		iter = iter+1;
	end
end
if (iter>=1000), fprintf('tSmo: too many iterations.\n '); end

