function [alphas,betas] = triple_smo(A,B,C,b,sol)
% Triple smo solves the following optimization problem:
%
% min  0.5 alpha' A alpha + alpha' B beta + 0.5 beta' C beta - alpha'b
% s.t. sum_i alpha_i = 1
%      alpha_i >= 0   forall i
%      beta_l  >= 0   forall l
%
%
% The algorithm works as follows:
%
% DO
%	choose a beta_m (heuristic)
%	choose alpha_i and alpha_j (heuristic)
%	optimize(beta_m,alpha_i,alpha_j)
% UNTIL some stopping criterion
%
% WARNING: 
% 	- until now, there is no heuristic for choosing the variables
%	  which makes the solving really slow
%	- ugly code
% 	- alpha-version (just to proof that it works)
%
warning('Triple SMO is only a proof of concept. Do not use it right now.');

n = size(A,1);
m = size(C,1);

% init with valid values
alphas = ones(n,1)/n; 
betas = zeros(m,1);

iters = 0;
numChanges = n*n*m;
while (iters<=100 && numChanges>0),
  numChanges = 0;
  incidents = 0;
  case1 = 0;
  case2 = 0;
  for i=1:m
    for k=1:n
      for l=k+1:n
        % choose alpha_k,alpha_l,beta_m
        check=0;
        if (k~=l),
          if (abs(C(i,i))>1e-8),
            [alphas,betas,check,hasChanged] = triple_smo_update1(k,l,i, alphas,betas, A,B,C,b);
            case1=case1+1;
          else 
            [alphas,betas,check,hasChanged] = triple_smo_update2(k,l,i, alphas,betas, A,B,C,b);            
            case2=case2+1;
          end
          %fprintf('i=%i,k=%i,l=%i   ak=%2.4f  al=%2.4f  bm=%2.4f   (%i)\n',i,k,l,alphas(k),alphas(l),betas(i),check);
          if (hasChanged), numChanges=numChanges+1; end
        end
      end
    end
  end

  fprintf('iters=%i  changes=%i/%i  diff2opt=%2.2f  case1/2=%i/%i\n',...
    iters,numChanges,n*n*m, sqrt(sum(([alphas;betas]-sol).^2)),case1,case2 );
  
 iters = iters+1;
end
fprintf('FINAL: iters=%i  diff2opt=%2.2f\n',iters,sqrt(sum(([alphas;betas]-sol).^2)));


