function [x, y] = pr_loqo2(c, H, A, b, l, u)
% pr_loqo2 -- optimizer for quadratic programs
%
% -- Synopsis
%    [x, y] = pr_loqo2(c, H, A, b, l, u)
%
% -- Input arguments
%    c        Linear coefficients
%    H        Quadratic coefficients 
%    A        Equality constraints
%    b        A * x = b
%    l        Lower bounds
%    u        Upper bounds
%
% -- Output arguments
%    x        Variables
%    y        Output
%
% -- Description
%    minimize   c' * x + 1/2 x' * H * x
%    subject to A*x = b
%               l <= x <= u
%
% -- Machine Learning for Computer Security Toolbox
%    Berlin Institute of Technology (TU Berlin)
%    $Id: pr_loqo2.m 10532 2010-03-06 19:32:01Z KonradRieck $
%
% originally written by Alex Smola

error('')
% the fudge factors
margin = 0.05;
bound  = 100;
sigfig_max = 8;
counter_max = 50;

% gather some starting data
[m, n] = size(A);
% this is done in order not to mess up H but only to tamper with H_x
H_x    = H;
H_diag = diag(H);

b_plus_1 = 1;
c_plus_1 = norm(c) + 1;
one_x = -ones(n,1);
one_y = -ones(m,1);

% starting point
for i = 1:n H_x(i,i) = H_diag(i) + 1; end;
H_y = eye(m);
c_x = c;
c_y = 0;

% and solve the system [-H_x A'; A H_y] [x, y] = [c_x; c_y]
R = chol(H_x);
H_Ac = R \ ([A; c_x'] / R)';
H_A = H_Ac(:,1:m);
H_c = H_Ac(:,(m+1):(m+1));
A_H_A = A * H_A;
A_H_c = A * H_c;
H_y_tmp = (A_H_A + H_y);
y = H_y_tmp \ (c_y + A_H_c);
x = H_A * y - H_c;

g = max(abs(x - l), bound);
z = max(abs(x), bound);
t = max(abs(u - x), bound);
s = max(abs(x), bound);

mu = (z' * g + s' * t)/(2 * n);

% set some default values
sigfig = 0;
counter = 0;
alfa = 1;

while ((sigfig < sigfig_max) * (counter < counter_max)),

    %update the iteration counter
    counter = counter + 1;

    %central path (predictor)
    H_dot_x = H * x;

    rho = - A * x + b;%%%
    nu = l - x + g;
    tau = u - x - t;
    sigma = c - A' * y - z + s + H_dot_x;

    gamma_z = - z;
    gamma_s = - s;

    % instrumentation
    x_dot_H_dot_x = x' * H_dot_x;

    primal_infeasibility = norm([tau; nu]) / b_plus_1;
    dual_infeasibility = norm([sigma]) / c_plus_1;

    primal_obj = c' * x + 0.5 * x_dot_H_dot_x;
    dual_obj = - 0.5 * x_dot_H_dot_x + l' * z - u' * s + b'*y; %%%

    old_sigfig = sigfig;
    sigfig = max(-log10(abs(primal_obj - dual_obj)/(abs(primal_obj) + 1)), 0);

    report = sprintf('counter %i p_i %e d_ii %e sigfig %f, alpha %f, p_o %f d_o %f mu %e', counter, primal_infeasibility, dual_infeasibility, sigfig, alfa, primal_obj, dual_obj, mu);
    %disp(report);

    % some more intermediate variables (the hat section)
    hat_nu = nu + g .* gamma_z ./ z;
    hat_tau = tau - t .* gamma_s ./ s;

    % the diagonal terms
    d = z ./ g + s ./ t;

    % initialization before the big cholesky
    for i = 1:n H_x(i,i) = H_diag(i) + d(i); end;
    H_y = 0;
    c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
    c_y = rho;

    % and solve the system [-H_x A'; A H_y] [delta_x, delta_y] = [c_x; c_y]
    R = chol(H_x);
    H_Ac = R \ ([A; c_x'] / R)';
    H_A = H_Ac(:,1:m);
    H_c = H_Ac(:,(m+1):(m+1));
    A_H_A = A * H_A;
    A_H_c = A * H_c;
    H_y_tmp = (A_H_A + H_y);
    delta_y = H_y_tmp \ (c_y + A_H_c);
    delta_x = H_A * delta_y - H_c;

    %backsubstitution
    delta_s = s .* (delta_x - hat_tau) ./ t;
    delta_z = z .* (hat_nu - delta_x) ./ g;

    delta_g = g .* (gamma_z - delta_z) ./ z;
    delta_t = t .* (gamma_s - delta_s) ./ s;

    %central path (corrector)
    gamma_z = mu ./ g - z - delta_z .* delta_g ./ g;
    gamma_s = mu ./ t - s - delta_s .* delta_t ./ t;

    % some more intermediate variables (the hat section)
    hat_nu = nu + g .* gamma_z ./ z;
    hat_tau = tau - t .* gamma_s ./ s;

    % the diagonal terms
    %d = z ./ g + s ./ t;

    % initialization before the big cholesky
    %for i = 1:n H_x(i,i) = H_diag(i) + d(i);
    %H_y = 0;
    c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
    c_y = rho;

    % and solve the system [-H_x A'; A H_y] [delta_x, delta_y] = [c_x; c_y]
    % R = chol(H_x);
    H_Ac = R \ ([A; c_x'] / R)';
    H_A = H_Ac(:,1:m);
    H_c = H_Ac(:,(m+1):(m+1));
    A_H_A = A * H_A;
    A_H_c = A * H_c;
    H_y_tmp = (A_H_A + H_y);
    delta_y = H_y_tmp \ (c_y + A_H_c);
    delta_x = H_A * delta_y - H_c;

    %backsubstitution

    delta_s = s .* (delta_x - hat_tau) ./ t;
    delta_z = z .* (hat_nu - delta_x) ./ g;

    delta_g = g .* (gamma_z - delta_z) ./ z;
    delta_t = t .* (gamma_s - delta_s) ./ s;

    %compute the updates
    alfa = - 0.95 / min([delta_g ./ g; delta_t ./ t;
        delta_z ./ z; delta_s ./ s; -1]);

    mu = (z' * g + s' * t)/(2 * n);
    mu = mu * ((alfa - 1) / (alfa + 10))^2;

    x = x + delta_x * alfa;
    g = g + delta_g * alfa;
    t = t + delta_t * alfa;
    y = y + delta_y * alfa;
    z = z + delta_z * alfa;
    s = s + delta_s * alfa;

end

%disp(report);


