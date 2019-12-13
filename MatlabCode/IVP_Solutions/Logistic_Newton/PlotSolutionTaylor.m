function PlotSolutionTaylor(a, L, nu, clr)
  % Plot solution from Taylor coefficients

  % Get the system size N
  N = length(a) - 1;

  % Get the values of t
  t = (-nu:0.001:nu)';

  % Number of t values
  m = length(t);

  % Compute the solution u
  u = zeros(m, 1);

  for k = 1:m
    u(k) = sum(a .* t(k) .^ (0:N)');
  end

  % Rescale t by L
  t = L * t;

  % plot solutoion
  plot(t, u, 'Color', clr, 'Linewidth', 3)
  xlabel('$$t$$', 'Interpreter', 'Latex', 'FontSize', 25)
  ylabel('$$x(t)$$', 'Interpreter', 'Latex', 'FontSize', 25)
end
