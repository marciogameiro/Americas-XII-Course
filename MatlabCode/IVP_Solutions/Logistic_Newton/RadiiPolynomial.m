function [r_min, r_max, Z2, Z1 , Z0, Y0] = RadiiPolynomial(a_bar, L, x0, nu)
  % Compute the radii polynomial and its roots

  % Get the system size N
  N = length(a_bar) - 1;

  % Compute the operator A_N
  DF_a_bar = DF(a_bar, L);
  A_N = inv(DF_a_bar);

  %%%%%%%%%% Compute Y0 %%%%%%%%%%

  % Pad a_bar with N+1 zeros for the correct
  % computation of the Cauchy product
  a_bar_N = [a_bar; zeros(N+1, 1)];

  F_a_N = F(a_bar_N, x0, L); % Compute F(a_bar_N)
  F_a_bar = F_a_N(1:N+1);    % Get the first N entries

  % Compute A * F(a_bar_N)
  AF_a_N = [A_N * F_a_bar; F_a_N(N+2:2*N+2) ./ (N+1:2*N+1)'];

  % Y0 is the weighted ell_1 norm of AF_a_N
  Y0 = WeightedL1VectorNorm(AF_a_N, nu);

  %%%%%%%%%% Compute Z0 %%%%%%%%%%

  % Compute the matrix B
  B = eye(N+1) - A_N * DF_a_bar;

  % Z0 is the weighted ell_1 norm of B
  Z0 = WeightedL1MatrixNorm(B, nu);

  %%%%%%%%%% Compute Z1 %%%%%%%%%%

  % Compute vector g
  g = -2 * a_bar;
  g(1) = g(1) + 1;

  % Compute the weighted ell_1 norm of g
  norm_g = WeightedL1VectorNorm(g, nu);

  % Z1 is a multiple of the norm of g
  Z1 = (L * nu / (N + 1)) * norm_g;

  %%%%%%%%%% Compute Z2 %%%%%%%%%%

  % Compute the weighted ell_1 norm of A_N
  norm_A_N = WeightedL1MatrixNorm(A_N, nu);

  % Compute upper bound for norm of A
  norm_A = max(norm_A_N, 1 / (N + 1));

  % Z2 is a multiple of the norm of A
  Z2 = 2 * L * nu * norm_A;

  %%%%%%%%%% Roots of radii poly %%%%%%%%%%

  % Get polynomial coefficients
  p_coeff = [Z2, Z1 + Z0 - 1, Y0];

  % Compute the roots
  p_roots = roots(p_coeff);

  if isreal(p_roots) % if rooots are real
    p_roots = sort(p_roots); % sort the roots
    r_min = p_roots(1); % get min r
    r_max = p_roots(2); % get max r

    if r_min >= 0 % if roots are positive
      disp('The proof was sucessful!')
    else
      disp('Failure! Radii polynomial has negative roots!')
    end
  else
    disp('Failure! Radii polynomial has complex roots!')
    r_min = -1;
    r_max = -1;
  end
end
