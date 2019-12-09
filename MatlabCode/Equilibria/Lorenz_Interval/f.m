function y = f(x, pars)
% Computes the Lorenz vector field

% Parameters
sigma = pars(1);
rho = pars(2);
beta = pars(3);

y = [sigma * (x(2) - x(1));
     rho * x(1) - x(2) - x(1) * x(3);
     -beta * x(3) + x(1) * x(2)];

end
