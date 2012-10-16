function [alphas,hasChanged] = simple_smo_update(k,l, alphas, A,b)

d = alphas(k) + alphas(l);

vtiltk = alphas'*A(:,k);
vtiltl = alphas'*A(:,l);

% calculate new optimal alpha_k
rho = A(k,k)-2*A(k,l)+A(l,l);
if (abs(rho)<1e-12),
  hasChanged = 3;
  return;
end

ak = alphas(k) + (b(k)-b(l)-vtiltk+vtiltl)/rho;
if isnan(ak),
  keyboard,
end

% check boundaries
if (ak<0), ak = 0.0; end
if (ak>d), ak = d; end;
al = d - ak;

% set sign if some value has changed within eps 
eps = 1e-4;
hasChanged = (abs(alphas(k)-ak)>eps);

% set new values
alphas(k) = ak;
alphas(l) = al;
