function [I, verified] = RigorousVerificationInterval(x_bar, pars)
% This function uses the radii polynomials approach to verify
% rigorously that there is solution to f(x)=0, near the
% approximate solution x_bar provided, where f is the function
% Lorenz vector field
%
% The function returns verified = true if the method was
% successful in verifying rigorously the existence of a solution.
% In this case I returns the interval of radii for which the radii
% polynomials are negative. If the rigorous verificatin failed, the
% function returns verified = false and I = [-1, -1];

% Compute A, the numerical inverse of Df
A = inv(Df(x_bar, pars));

% Compute the Y0 bound
Y0 = norm(A * f(x_bar, pars), Inf);

% Compute the Z0 bound
Z0 = norm(eye(3) - A * Df(x_bar, pars), Inf);

% Compute the Z2 bound
Z2 = 2 * max([abs(A(1,2)) + abs(A(1,3)), ...
              abs(A(2,2)) + abs(A(2,3)), ...
              abs(A(3,2)) + abs(A(3,3))]);

% Define the radii polynomial
p = [Z2, -(1-Z0), Y0];

% Compute and sort the roots
r = sort(roots(p));

% Need to decrease the interval
% by a small value delta_r
delta_r = 1e-15;

% Compute the interval I
I = [r(1) + delta_r, r(2) - delta_r];

% Now we need to check that the radii polynomials
% are negative in I.

% To get a rigorous proof we use intrval arithmetic
% We need to recompute Y0, Z0 and Z2

% Redefine x_bar as interval
x_bar = intval(x_bar);

% Redefine pars as interval
% Careful, we may get the wrong
% parameter values if we define
% pars = intval(pars);
sigma = intval(10);
rho = intval(28);
beta = intval(8)/intval(3);
pars = [sigma rho beta];

% Redefine A as interval
% No need to recompute A
A = intval(A);

% Redefine Y0 as interval
Y0 = norm(A * f(x_bar, pars), Inf);

% Redefine Z0 as interval
Z0 = norm(eye(3) - A * Df(x_bar, pars), Inf);

% Redefine Z2 as interval
Z2 = 2 * max([abs(A(1,2)) + abs(A(1,3)), ...
              abs(A(2,2)) + abs(A(2,3)), ...
              abs(A(3,2)) + abs(A(3,3))]);

% Evaluate the raddi polynomials at the
% end points using interval arithmetic
r1 = intval(I(1));
r2 = intval(I(2));
p1 = Z2 * r1^2 -(1 - Z0) * r1 + Y0;
p2 = Z2 * r2^2 -(1 - Z0) * r2 + Y0;

% Check if they are negative
% if p1 >= 0  % p1 not negative
if sup(p1) >= 0  % p1 not negative
  I = [-1, -1];
  verified = false;
  return;
end

% if p2 >= 0  % p1 not negative
if sup(p2) >= 0  % p1 not negative
  I = [-1, -1];
  verified = false;
  return;
end

% Solution verified if we get here
verified = true;

end
