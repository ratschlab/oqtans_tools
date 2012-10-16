function [alphas,betas] = triple_smo_grad_wss(A,B,C,b,sol, w_part,U,d)
% Fast Triple SMO solves the following optimization problem:
%
% min  0.5 alpha' A alpha + alpha' B beta + 0.5 beta' C beta - alpha'b
% s.t. sum_i alpha_i = 1
%      alpha_i >= 0   forall i
%      beta_l  >= 0   forall l
%
%
% The algorithm uses maximum gradients for the working set selection:
%
% DO
%   choose alpha_i and alpha_j
%   optimize(alpha_i,alpha_j)
%   calculate new betas
% UNTIL some stopping criterion
%
%
% WARNING: 
% - alpha-version 
%

eps = 1e-6;

num_alphas = size(A,1);
num_betas = size(C,1);

% init with valid values
alphas = ones(num_alphas,1)/num_alphas; 
betas = zeros(num_betas,1);

inds = find(diag(abs(full(C))) >= 1e-16);

%fprintf('The new Fast TripleSMO.\n');

iter = 0;
obj = inf;
dobj = inf;

betas(betas<0) = 0;

if (num_alphas==1),
	eta = 0.01;
	diff = 1;
  iter = 0;
	while (diff>eps),
		diff = 0;
		for i=1:length(inds),
			D = C(inds,inds)*betas(inds) + B(:,inds)';
			m = inds(i);
			foo = betas(m);
			betas(m) = betas(m) - eta*D(i);
			if (betas(m)<0), betas(m)=0; end
			diff = diff + abs(betas(m)-foo);
		end
		iter = iter+1;
	end
	fprintf('%i full gradient descent steps.\n',iter);

else
  noChanges = 0;
	while (iter<=1000 && dobj>eps),

		noChanges = 1;

		for k=1:num_alphas-1,
		  for l=k+1:num_alphas,
			%ainds = 1:num_alphas;
			%ainds = ainds(ainds~=k);
			%rpi = randperm(length(ainds));
			%ainds = ainds(rpi);
			%l = ainds(1);
				[alphas,hasChanged] = fast_tsmo_update(k,l, alphas,betas, A,B,C,b);            
  	    if (hasChanged), noChanges=0; end;
			end
		end

		eta = 0.01;
		diff = 1;
		while (diff>eps),
			diff = 0;
		 	for i=1:length(inds),
				D = C(inds,inds)*betas(inds) + B(:,inds)'*alphas;
				m = inds(i);
				foo = betas(m);
				betas(m) = betas(m) - eta*D(i);
				if (betas(m)<0), betas(m)=0; end
				diff = diff + abs(betas(m)-foo);
			end
		end

		oldObj = obj;
		obj = 0.5*alphas'*A*alphas + alphas'*B*betas + 0.5*betas'*C*betas - alphas'*b';
		dobj = oldObj-obj;
		% increase loop counter 
		iter = iter+1;
	end
end

if (iter>=1000), fprintf('tSmo: too many iterations.\n '); end

%fprintf('FINAL: iter=%i  diff2opt=%2.2f\n',iter,sqrt(sum(([alphas;betas]-sol).^2)));

