
% This script uses Newton's method implemented in the file Newton.m
% to compute four numerical zeros of the function in the file f.m,
% for lambda = 3, and saves the solutions to the file solns.mat

% Initial guesses for Newton
x1_0 = [-1; 1];
x2_0 = [1; 1];
x3_0 = [1; -1];
x4_0 = [-1; -1];

lambda = 3;         % Value of parameter to compute solution
max_itr = 20;       % Maximum number of Newton iterations
err_tol = 1.0e-10;  % Error tolerance for Newton
% err_tol = 1.0e-14;  % Error tolerance for Newton

% Run Newton's method using x1_0 as initial guess to compute solution x1
[x1_bar, num_itr1, converged] = Newton(@f, @Df, x1_0, lambda, max_itr, err_tol);

% Stop if Newton failed
if converged == false
  disp('Newton did not convege.');
  return;
end

% Run Newton's method using x2_0 as initial guess to compute solution x2
[x2_bar, num_itr2, converged] = Newton(@f, @Df, x2_0, lambda, max_itr, err_tol);

% Stop if Newton failed
if converged == false
  disp('Newton did not convege.');
  return;
end

% Run Newton's method using x3_0 as initial guess to compute solution x3
[x3_bar, num_itr3, converged] = Newton(@f, @Df, x3_0, lambda, max_itr, err_tol);

% Stop if Newton failed
if converged == false
  disp('Newton did not convege.');
  return;
end

% Run Newton's method using x4_0 as initial guess to compute solution x4
[x4_bar, num_itr4, converged] = Newton(@f, @Df, x4_0, lambda, max_itr, err_tol);

% Stop if Newton failed
if converged == false
  disp('Newton did not convege.');
  return;
end

% Print the result of computations
disp(' ')
disp(['For the initial guess x1_0 Newton converged in ' num2str(num_itr1) ' iterations to x1 = (' num2str(x1_bar') ')'])
disp(['For the initial guess x2_0 Newton converged in ' num2str(num_itr2) ' iterations to x2 = (' num2str(x2_bar') ')'])
disp(['For the initial guess x3_0 Newton converged in ' num2str(num_itr3) ' iterations to x3 = (' num2str(x3_bar') ')'])
disp(['For the initial guess x4_0 Newton converged in ' num2str(num_itr4) ' iterations to x4 = (' num2str(x4_bar') ')'])

% Save solution to solns.mat
save solns.mat lambda x1_bar x2_bar x3_bar x4_bar
