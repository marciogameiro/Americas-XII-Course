function norm_A = WeightedL1MatrixNorm(A, nu)
  % Compute the weighted ell_1 matrix norm
  % Assumes that A is a square matrix

  % Number of rows of A
  n = size(A, 1);

  % Compute the weights nu^k
  omega = nu .^ (0:n-1);

  % Compute the weighted norm of A
  norm_A = max((omega * abs(A)) ./ omega);
end
