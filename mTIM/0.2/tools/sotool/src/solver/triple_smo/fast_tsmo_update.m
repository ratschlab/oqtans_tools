function [alphas,hasChanged] = fast_tsmo_update(k,l, alphas,betas, A,B,C,b)

d = alphas(k) + alphas(l);

htiltk = B(k,:)*betas;
htiltl = B(l,:)*betas;

vtiltk = alphas'*A(:,k);
vtiltl = alphas'*A(:,l);

% calculate new optimal alpha_k
rho = A(k,k)-2*A(k,l)+A(l,l);
ak = alphas(k) + (b(k)-b(l)-htiltk+htiltl-vtiltk+vtiltl)/rho;

% check boundaries
if (ak<0), ak = 0.0; end
if (ak>d), ak = d; end;
al = d - ak;

% set sign if some value has changed within eps 
eps = 1e-4;
hasChanged = (abs(alphas(k)-ak)>eps || abs(alphas(l)-al)>eps);

% set new values
alphas(k) = ak;
alphas(l) = al;
