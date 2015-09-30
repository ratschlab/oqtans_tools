function [params, w] = ls_constrained_parabola(params, Q,losse,w,dpsis,dpsisi,U)
% Parabola linesearch with constraints w'U <= 0.
% Nico Goernitz'2010

params.lsIters = 0;
params.lsNstar = 0.0;

lossParams = params.lossParams;
loss_fct = params.loss_fct;

% if this is the first function call and the last valid w is
% therefore unavailable store the current w and exit.
if (params.lsW == -1)
    params.lsW = w;
    return;
end

% the 3 parabel points
dir = (params.lsW-w);
na = 0.0; a = dir*na + w; 
nb = 0.5; b = dir*nb + w;
nc = 1.0; c = dir*nc + w;

% calculate function values
% constraint U*w<=0 is fullfilled (convex set and 2 valid solutions)
[slacks,slacksIdxs,lossParams] = loss_fct(lossParams,a,losse,dpsis,dpsisi);
obj_a = 0.5*a'*Q*a + sum(slacks);
[slacks,slacksIdxs,lossParams] = loss_fct(lossParams,b,losse,dpsis,dpsisi);
obj_b = 0.5*b'*Q*b + sum(slacks);
[slacks,slacksIdxs,lossParams] = loss_fct(lossParams,c,losse,dpsis,dpsisi);
obj_c = 0.5*c'*Q*c + sum(slacks);

% obj_c should be the lowest value
if (obj_c > obj_a)
  foo=obj_c; obj_c=obj_a; obj_a=foo;
  foo=c; c=a; a=foo;
  foo=nc; nc=na; na=foo;
end
if (obj_c > obj_b)
  foo=obj_c; obj_c=obj_b; obj_b=foo;
  foo=nc; nc=nb; nb=foo;
end

params.lsObj = obj_c;
params.lsW = dir*nc + w;
params.lsNstar = nc;


% iterate
hasConverged = 0;
iter = 0;
while (iter<params.lsSteps && ~hasConverged)
  % abbriviation
  naa = na*na;
  nbb = nb*nb;
  ncc = nc*nc;
  
  % calculate parameters of the parabola according 
  % to the formula: f_x = p1 x^2 + p2 x + p3
  if (abs(na)<eps)
    % if 'na' is zero -> '1*p3 = obj_a'
    m = (nbb*nc - nb*ncc);
    n = (nbb - ncc);
    o = obj_c*nbb - obj_b*ncc;    
    
    p3 = obj_a;
    p2 = (o - n*p3) / m;
    p1 = (obj_b - p3 - nb*p2) / nbb;
  elseif (abs(nb)<eps)
    % if 'nb' is zero -> '1*p3 = obj_b'
    m = (naa*nc - na*ncc);
    n = (naa - ncc);
    o = obj_c*naa - obj_a*ncc;    
    
    p3 = obj_b;
    p2 = (o - n*p3) / m;
    p1 = (obj_a - p3 - na*p2) / naa;
  elseif (abs(nc)<eps)
    % if 'nc' is zero -> '1*p3 = obj_c'
    m = (naa*nb - na*nbb);
    n = (naa - nbb);
    o = obj_b*naa - obj_a*nbb;    
    
    p3 = obj_c;
    p2 = (o - n*p3) / m;
    p1 = (obj_a - p3 - na*p2) / naa;
  else
    % do full calculation: na,nb,nc ~= 0
    m = nb*naa - na*nbb;
    n = naa - nbb;
    o = obj_b*naa - obj_a*nbb;
    
    p = nc*naa - na*ncc;
    q = naa - ncc;
    r = obj_c*naa - obj_a*ncc;
    
    p3 = (r*m - o*p)/(q*m - n*p);
    p2 = (o - n*p3) / m;
    p1 = (obj_a - na*p2 - p3) / naa;
  end
  
  % new optimal step size
  nstar = -p2/(2*p1);

  % new position
  wstar = nstar*dir + w;
   
  % next iteration
  iter = iter + 1;

  % check if the constraints are fullfilled
  % at the new possition
  if (~isempty(U) && any(U*wstar>0))
    % constraints not fullfilled -> 
    % store old values and return
    break;
  else
    % constraints are fullfilled ->
    % calculate new objective function
    % and store new position, start next iteration
    [slacks,slacksIdxs,lossParams] = loss_fct(lossParams,wstar,losse,dpsis,dpsisi);
    objstar = 0.5*wstar'*Q*wstar + sum(slacks);
  
    % check if converged:
    % difference in objective is now small enough 
    % to break up
    if (abs(params.lsObj-objstar)<params.lsEps)
      hasConverged = 1;
    end

    % is the new objective value at wstar smaller
    % than obj_c (last lowest value)?
    if (objstar<obj_c) 
      params.lsObj = objstar;
      params.lsW = wstar;
      params.lsNstar = nstar;

      % change values = skip highest value, 
      % store new lowest value at c
      if (obj_a>obj_b)
        na=nb; a=b; obj_a=obj_b;
        nb=nc; b=c; obj_b=obj_c;
        nc=nstar; c=wstar; obj_c=objstar;
      else
        nb=na; b=a; obj_b=obj_a;
        na=nc; a=c; obj_a=obj_c;
        nc=nstar; c=wstar; obj_c=objstar;
      end
    else
      % objstar is not smaller than obj_c
      if (obj_a>obj_b)
        na=nstar; a=wstar; obj_a=objstar;
      else
        nb=nstar; b=wstar; obj_b=objstar;
      end
    end
  end  
end

params.lsIters = iter;
% the new solution is somewhere in between the 
% optimal solution obtained by the bundle method and
% the linesearch solution.
w = (1-params.lsTheta)*params.lsW + params.lsTheta*w;
