function [I, verified] = RigorousVerification(x_bar, lambda)
% This function uses the radii polynomials approach to verify
% rigorously that there is solution to f(x)=0, near the
% approximate solution x_bar provided, where f is the function
% defined on the file f.m
%
% The function returns verified = true if the method was
% successful in verifying rigorously the existence of a solution.
% In this case it returns the interval I of radii for which the radii
% polynomials are negative. If the rigorous verificatin failed, the
% function returns verified = false and I = [-1, -1];


% Compute the matrix A
% A = inv(Df(x_bar,lambda));
A = (1 / (16 * x_bar(1) * x_bar(2) - 1)) * [2 * x_bar(2) -1; -1 8 * x_bar(1)];

% Compute the Y0 bound using the sup norm
Y0 = norm(A * f(x_bar,lambda), Inf);

% Compute the Z0 bound using the sup norm
% Z0 = norm(eye(2) - A*Df(x_bar,lambda), Inf);
Z0 = 0; % Since we use the exact inverse

% Compute the Z2 bound using the sup norm
Z2 = max(16*abs(x_bar(2))+2, 8+16*abs(x_bar(1)))/(abs(16*x_bar(1)*x_bar(2)-1));
  
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

% Evaluate the raddi polynomials at the end points
r1 = I(1); r2 = I(2);
p1 = Z2 * r1^2 -(1 - Z0) * r1 + Y0;
p2 = Z2 * r2^2 -(1 - Z0) * r2 + Y0;

% Check if they are negative
if p1 >= 0  % p1 not negative
  I = [-1, -1];
  verified = false;
  return;
end

if p2 >= 0  % p1 not negative
  I = [-1, -1];
  verified = false;
  return;
end

% Solution verified if we get here
verified = true;

end
