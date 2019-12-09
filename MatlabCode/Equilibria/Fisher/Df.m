function Df_x = Df(x, lambda)
  N = length(x) - 1;

  Df_x = zeros(N + 1, N + 1); % Jacobian matrix Df(x)

  for k = 0:N
    for l = 0:N
      % Compute the derivative of the non-linear term
      if (k + l <= N) && (l > 0)
        Df_x(k + 1, l + 1) = - 2.0 * lambda * (x(abs(k - l) + 1) + x(k + l + 1));
      else
        Df_x(k + 1, l + 1) = - 2.0 * lambda * x(abs(k - l) + 1);
      end
    end

    % Add derivative of linear term to diagonal
    Df_x(k + 1, k + 1) = Df_x(k + 1, k + 1) + lambda - k^2;
  end
end
