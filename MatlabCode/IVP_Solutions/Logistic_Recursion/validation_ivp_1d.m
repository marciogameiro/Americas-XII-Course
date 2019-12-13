function [I] = validation_ivp_1d(x0,L,nu,N)

a = zeros(N+1,1);

a(1)=x0;

for n=2:N+1
    a(n)=-L/(n-1)*(a(n-1)-sum(a(1:n-1).*a(n-1:-1:1)));
end

norm_a = sum(abs(a).*nu.^(0:N)');

a = [a;zeros(N,1)];

a2 = quadratic_cauchy_product(a,a);

n = (N:2*N)';
Y0 = L*nu*sum(1./(n+1).*abs(-a(n+1)+a2(n+1)).*nu.^n);

Z1 = (1+2*norm_a)*L*nu/(N+1);
Z2 = 2*L*nu/(N+1);

fprintf('\n')
display(['Z2 = ',num2str(Z2),', Z1 = ',num2str(Z1),', Y0 = ',num2str(Y0)])
fprintf('\n')

p = [Z2 -(1-Z1) Y0];
disc = (1-Z1)^2-4*Y0*Z2;

I = sort(roots(p));

if disc>0
    plot_a_taylor(a,x0,L,nu);
end

end

