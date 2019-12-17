
L = 3.0;
N = 300;

max_itr = 20;
err_tol = 1e-14;

x0 = 1/4; % 5.0; % 1/4;

a0 = zeros(N + 1, 1);
a0(1) = x0;

% Define new functions f and Df to be F and
% DF as functions of only the variable a
f = @(a) F(a, x0, L);
Df = @(a) DF(a, L);

[a, num_itr, converged] = Newton(f, Df, a0, max_itr, err_tol);

nu = 1.0;

[r_min, r_max, Z2, Z1 , Z0, Y0] = RadiiPolynomial(a, L, x0, nu);

clr = [1 0 0]; % Red
PlotSolutionTaylor(a, x0, L, nu, clr)
