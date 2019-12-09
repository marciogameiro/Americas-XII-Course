function [x1, num_itr, converged] = Newton(f, Df, x0, lambda, max_itr, err_tol)
  x1 = x0; % Initialize x1 as the initial guess x0

  f_x1 = f(x1, lambda);   % Compute f at the initial guess
  norm_f_x1 = norm(f_x1); % Norm of f at the initial guess

  num_itr = 0;            % Intialize iteration counter

  % Set converged to true (and return x1 as the solution)
  % if norm of f_x1 is less than tolerance
  if norm_f_x1 < err_tol
    converged = true;
    return
  end

  %%%%%%%%%%%%%%% Start Newton method %%%%%%%%%%%%%%%
  display(['Norm of f(x0) = ' num2str(norm_f_x1)]);

  % Start Newton iterations
  while num_itr < max_itr
    x0 = x1;                 % Set x0 as x1
    f_x0 = f_x1;             % Set f_x0 as f_x1

    Df_x0 = Df(x0, lambda);  % Jacobian matrix at x0
    u = Df_x0 \ f_x0;        % Newton iteration (solve linear sistem Df_x0 * u = f_x0)
    x1 = x0 - u;             % Newton iteration (compute x1)

    f_x1 = f(x1, lambda);    % Compute f at x1
    norm_f_x1 = norm(f_x1);  % Norm of f at x1

    display(['Norm of f(x) = ' num2str(norm_f_x1)]);

    num_itr = num_itr + 1;   % Increment iteration counter

    % Check for convergence
    if norm_f_x1 < err_tol
      converged = true;
      return
    end
  end

  % Newton did not converged if we got here
  % Set converged to false (and return)
  converged = false;
end
