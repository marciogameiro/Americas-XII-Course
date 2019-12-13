function F_a = F(a, x0, L)
  % The function F(a)

  % Get the system size N
  N = length(a) - 1;

  % Initilaize F_a as zero
  F_a = zeros(N + 1, 1);

  % Compute the Cauchy product
  a2 = zeros(N + 1, 1);

  for k = 0:N
    a2(k + 1) = sum(a(1:k+1) .* a(k+1:-1:1));
  end

  % % Compute using a double loop
  % for k = 0:N
  %   for j = 0:k
  %     a2(k+1) = a2(k+1) + a(j+1) * a(k-j+1);
  %   end
  % end

  % Set first entry of F
  k = 0;
  F_a(k + 1) = a(k + 1) - x0;

  % Set remaning entries of F
  for k = 1:N
    F_a(k+1) = k * a(k+1) + L * (a(k) - a2(k));
  end
end
