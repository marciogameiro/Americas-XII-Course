function J = Df(x, pars)
% Computes the Jacobian matrix of the Lorenz vector field

% Parameters
sigma = pars(1);
rho = pars(2);
beta = pars(3);

J = [-sigma, sigma, 0;
     rho - x(3), -1, -x(1);
     x(2), x(1), -beta];

end
