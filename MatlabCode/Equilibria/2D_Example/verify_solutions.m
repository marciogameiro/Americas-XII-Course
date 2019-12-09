
% Read the approximate solutions from the file solns.mat and call
% the function RigorousVerification to rigorously verify the solutions

% Read the paramter value lambda and the approximate
% solutions x1, x2, x3, and x4 from the file solns.mat
load solns.mat

% Try to rigorously verify the solution x1_bar
[I1, verified] = RigorousVerification(x1_bar, lambda);

if verified == true
  disp('Verification was successful for x1_bar!');
  disp(['x1_bar = (' num2str(x1_bar(1), 16) ', ' num2str(x1_bar(2), 16) ')']);
  disp(['I1 = (' num2str(I1(1), 16) ', ' num2str(I1(2), 16) ')']);
else
  disp('Verification Falied for x1_bar!');
end

% Try to rigorously verify the solution x2_bar
[I2, verified] = RigorousVerification(x2_bar, lambda);

if verified == true
  disp('Verification was successful for x2_bar!');
  disp(['x2_bar = (' num2str(x2_bar(1), 16) ', ' num2str(x2_bar(2), 16) ')']);
  disp(['I2 = (' num2str(I2(1), 16) ', ' num2str(I2(2), 16) ')']);
else
  disp('Verification Falied for x2_bar!');
end

% Try to rigorously verify the solution x3_bar
[I3, verified] = RigorousVerification(x3_bar, lambda);

if verified == true
  disp('Verification was successful for x3_bar!');
  disp(['x3_bar = (' num2str(x3_bar(1), 16) ', ' num2str(x3_bar(2), 16) ')']);
  disp(['I3 = (' num2str(I3(1), 16) ', ' num2str(I3(2), 16) ')']);
else
  disp('Verification Falied for x3_bar!');
end

% Try to rigorously verify the solution x4_bar
[I4, verified] = RigorousVerification(x4_bar, lambda);

if verified == true
  disp('Verification was successful for x4_bar!');
  disp(['x4_bar = (' num2str(x4_bar(1), 16) ', ' num2str(x4_bar(2), 16) ')']);
  disp(['I4 = (' num2str(I4(1), 16) ', ' num2str(I4(2), 16) ')']);
else
  disp('Verification Falied for x4_bar!');
end

% Next plot the enclosures for each solution
X = [x1_bar, x2_bar, x3_bar, x4_bar];  % All solutions
Iks = [I1; I2; I3; I4];                % The intervals for all solutions

% Plot the enclosures
PlotEnclosures(X, Iks);
