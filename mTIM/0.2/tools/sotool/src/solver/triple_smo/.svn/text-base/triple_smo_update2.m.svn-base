function [alphas,betas,check,hasChanged] = triple_smo_update2(k,l,m, alphas,betas, A,B,C,b)
check = 0;
d = alphas(k) + alphas(l);

htiltk = B(k,:)*betas;
htiltl = B(l,:)*betas;

vtiltk = alphas'*A(:,k);
vtiltl = alphas'*A(:,l);


% calculate new optimal alpha_k
Bklm = B(k,m)-B(l,m);
rho = A(k,k)-2*A(k,l)+A(l,l);

ak = alphas(k) ...
  + betas(m)*Bklm/rho ...
  + (b(k)-b(l)-htiltk+htiltl-vtiltk+vtiltl)/rho;

% check boundaries
if (ak<0), ak = 0.0; end
if (ak>d), ak = d; end;
al = d - ak;

% set new beta_m
bm = 0.0;

% set sign if some value has changed within eps 
eps = 1e-2;
hasChanged = (abs(alphas(k)-ak)>eps || abs(alphas(l)-al)>eps || abs(betas(m)-bm)>eps);

% set new values
alphas(k) = ak;
alphas(l) = al;
betas(m) = bm;
