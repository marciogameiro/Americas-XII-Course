function DF_a = DF(a, L)
  % Compute the derivative DF(a)

  % Get the system size N
  N = length(a) - 1;

  % Initilaize DF_a as zero
  DF_a = zeros(N + 1, N + 1);

  % Set non-linear terms
  for k1 = 1:N  % rows 1:N
    for k2 = 0:k1-1  % columns 0:k1-1
      DF_a(k1 + 1, k2 + 1) = -2.0 * L * a(k1 - k2);
    end
  end

  % Add main diagonal part and lower diagonal part
  DF_a = DF_a + diag([1, 1:N]) + L * diag(ones(1, N), -1);
end
