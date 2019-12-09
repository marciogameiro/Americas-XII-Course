
% Calls the function RigorousVerification to rigorously verify the
% existence of a solution in a ball about the approximated equilibrium
% solution for the Lorenz vector field

% Define the parameter values for Lorenz
sigma = 10; rho = 28; beta = 8/3;

% Parameter for vector field
pars = [sigma rho beta];

% Set the value of an equilibria for Lorenz
x_bar = [8.4853; 8.4853; 27];
% x_bar = [8.4; 8.4; 26.5];
% x_bar = [sqrt(beta * (rho - 1)); sqrt(beta * (rho - 1)); rho - 1];

% Try to rigorously verify the solution x_bar
[I, verified] = RigorousVerification(x_bar, pars);

% Display a message and finish if verificatioin failed
if verified == false
  disp('Verification Falied!');
  return;
end

% If we got here, verificatioin was successful
disp('Verification was successful!');

disp(' ')

% Print information about the solution
disp('The verified numerical approximation is:');
disp(['x = (' num2str(x_bar(1), 16) ', ' num2str(x_bar(2), 16) ', ' num2str(x_bar(3), 16) ')'])

disp(' ')

disp('The existence interval of verified radii is:');
disp(['I = (' num2str(I(1), 16) ', ' num2str(I(2), 16) ')'])
