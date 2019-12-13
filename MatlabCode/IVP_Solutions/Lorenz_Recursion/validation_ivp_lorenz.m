function [I] = validation_ivp_lorenz(x0,L,nu,N,par)

sigma = par(1);
rho = par(2);
beta = par(3);

a1 = zeros(N+1,1);
a2 = zeros(N+1,1);
a3 = zeros(N+1,1);

a1(1)=x0(1);
a2(1)=x0(2);
a3(1)=x0(3);

for k=2:N+1
    a1(k)=L/(k-1)*sigma*(a2(k-1)-a1(k-1));
    a2(k)=L/(k-1)*(rho*a1(k-1)-a2(k-1)-sum(a1(1:k-1).*a3(k-1:-1:1)));
    a3(k)=L/(k-1)*(-beta*a3(k-1)+sum(a1(1:k-1).*a2(k-1:-1:1)));
end

norm_a1 = sum(abs(a1).*nu.^(0:N)');
norm_a2 = sum(abs(a2).*nu.^(0:N)');
norm_a3 = sum(abs(a3).*nu.^(0:N)');

a1 = [a1;zeros(N,1)];
a2 = [a2;zeros(N,1)];
a3 = [a3;zeros(N,1)];

a1a3 = quadratic_cauchy_product(a1,a3);
a1a2 = quadratic_cauchy_product(a1,a2);

k = (N:2*N)';
Y0_1 = L/(N+1)*abs(sigma*(a2(N+1)-a1(N+1)))*nu^(N+1);
Y0_2 = L*nu*sum(1./(k+1).*abs(rho*a1(k+1)-a2(k+1)-a1a3(k+1)).*nu.^k);
Y0_3 = L*nu*sum(1./(k+1).*abs(-beta*a3(k+1)+a1a2(k+1)).*nu.^k);
Y0 = max([Y0_1 Y0_2 Y0_3]);

z0_1 = 2*L*abs(sigma);
z0_2 = L*(abs(rho)+1+norm_a1+norm_a3);
z0_3 = L*(abs(beta)+norm_a1+norm_a2);
Z1 = nu/(N+1)*max([z0_1 z0_2 z0_3]);
Z2 = 2*L*nu/(N+1);

% fprintf('\n')
% display(['Z2 = ',num2str(Z2),', Z1 = ',num2str(Z1),', Y0 = ',num2str(Y0)])
% fprintf('\n')

p = [Z2 -(1-Z1) Y0];
disc = (1-Z1)^2-4*Y0*Z2;

I = sort(roots(p));

if disc>0 && Z1<1
    plot_a_taylor_lorenz(a1,a2,a3,x0,L,nu)
    display(I)
else
    disp('Failure!')
end

end

