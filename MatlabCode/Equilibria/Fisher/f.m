function f_x = f(x, lambda)
  N = length(x) - 1;

  x2 = zeros(size(x));  % The quadratic term
  f_x = zeros(size(x)); % The computed f(x)

  for k = 0:N
    % Compute the quadratic term
    for k1 = -N:N
      k2 = k - k1;
      if abs(k2) <= N
        x2(k + 1) = x2(k + 1) + x(abs(k1) + 1) * x(abs(k2) + 1);
      end
    end

    f_x(k + 1) = (lambda - k^2) * x(k + 1) - lambda * x2(k + 1);
  end
end
