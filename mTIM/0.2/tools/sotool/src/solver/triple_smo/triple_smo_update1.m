function [alphas,betas,check,hasChanged] = triple_smo_update1(k,l,m, alphas,betas, A,B,C,b)

d = alphas(k) + alphas(l);

ntilt = alphas'*B(:,m);
wtilt = betas'*C(:,m);

htiltk = B(k,:)*betas;
htiltl = B(l,:)*betas;

vtiltk = alphas'*A(:,k);
vtiltl = alphas'*A(:,l);


% calculate new optimal alpha_k
Bklm = B(k,m)-B(l,m);
Akl = A(k,k)-2*A(k,l)+A(l,l);
rho = Bklm*Bklm - Akl*C(m,m);

ak = alphas(k) ...
  - (Bklm*(ntilt+wtilt))/rho ...
  + C(m,m)*(htiltk-htiltl+vtiltk-vtiltl-b(k)+b(l))/rho;

% check boundaries
if (ak<0), ak = 0.0; end
if (ak>d), ak = d; end;
al = d - ak;

% calculate new beta_m
bm = betas(m) - (ntilt+wtilt- (alphas(k)-ak)*B(k,m) - (alphas(l)-al)*B(l,m))/C(m,m);

check = 0;
if (bm<0.0), bm = 0.0; end

% set sign if some value has changed within eps 
eps=1e-2;
hasChanged = (abs(alphas(k)-ak)>eps || abs(alphas(l)-al)>eps || abs(betas(m)-bm)>eps);

% set new values
alphas(k) = ak;
alphas(l) = al;
betas(m) = bm;


