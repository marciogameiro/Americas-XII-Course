function PlotEnclosures(X, Iks)

% Number of solutions
num_solns = size(X, 2);

figure; hold on;

% Plot an enclosure for each solution
for k = 1:num_solns
  x = X(:, k);     % k-th solution
  I = Iks(k, :);  % Interval of the k-th solutions

  R = I(2);  % Largest value of r

  % Plot the point corresponding to the solution
  plot(x(1), x(2), '.', 'MarkerSize', 10, 'Color', [0 0 0])

  % Compute the four corners of the square
  p1 = x + R * [-1; -1];
  p2 = x + R * [1; -1];
  p3 = x + R * [1; 1];
  p4 = x + R * [-1; 1];

  % Plot the sides of the square
  plot([p1(1) p2(1)], [p1(2) p2(2)], 'LineWidth', 2);  % Line connecting p1 and p2
  plot([p2(1) p3(1)], [p2(2) p3(2)], 'LineWidth', 2);  % Line connecting p2 and p3
  plot([p3(1) p4(1)], [p3(2) p4(2)], 'LineWidth', 2);  % Line connecting p3 and p4
  plot([p4(1) p1(1)], [p4(2) p1(2)], 'LineWidth', 2);  % Line connecting p4 and p1
end

% Add labels to axis
xlabel('x', 'FontSize', 18)
ylabel('y', 'FontSize', 18)

% Set the axis size
axis([-2.5 2.5 -2.5 2.5])

% Save plot to a file
print('-depsc2', 'enclosures_2d_ex.eps')
