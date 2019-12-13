function norm_a = WeightedL1VectorNorm(a, nu)
  % Compute the weighted ell_1 vector norm

  % Get the vector length
  n = length(a);

  % Compute the weights nu^k
  omega = nu .^ (0:n-1);

  % Compute the weighted norm of a
  norm_a = sum(abs(a) .* omega');
end
