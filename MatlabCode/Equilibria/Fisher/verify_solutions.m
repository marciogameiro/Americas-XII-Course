% Compute a numerical solutions and call the function
% RigorousVerification to rigorously verify the existence
% of a true solution close to the approximated solution

% Set the parameter values
lambda = 50; N = 100;

% Initial guess for Newton
x0 = zeros(N + 1, 1);
x0(1) = 1.0;
x0(2) = 0.3;

max_itr = 100;    % Max number of Newton iterations
err_tol = 1e-10;  % Error tolerance for Newton

% Compute x_bar using Newton's method
x_bar = Newton(@f, @Df, x0, lambda, max_itr, err_tol);

% Try to rigorously verify the solution x_bar
[I, verified] = RigorousVerification(x_bar, lambda);

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
format long
x_bar

disp(' ')

disp('The existence interval of verified radii is:');
disp(['I = (' num2str(I(1), 16) ', ' num2str(I(2), 16) ')'])

plot(x_bar, '.')
