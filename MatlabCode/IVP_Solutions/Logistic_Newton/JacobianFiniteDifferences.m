function J = JacobianFiniteDifferences(f, x, h)
  % Compute the Jocobian matrix using finite differences

  m = size(x, 1);

  J = zeros(m);
  I = eye(m);

  for j = 1:m
    x_h = x + h * I(:, j);
    J(:, j) = (f(x_h) - f(x)) / h;
  end
end
