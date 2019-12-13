function [I] = validation_ivp_bootstrap_1d(x0,L,nu,N)

a = zeros(N+1,1);

a(1)=x0;

for n=2:N+1
    a(n)=-L/(n-1)*(a(n-1)-sum(a(1:n-1).*a(n-1:-1:1)));
end

norm_a = sum(abs(a).*nu.^(0:N)');

a  = [a;zeros(2*N,1)];
a2 = quadratic_cauchy_product(a,a);
a3 = quadratic_cauchy_product(a2,a);

n = (N-1:3*N)';
Y0 = (L*nu)^2*sum(1./((n+2).*(n+1)).*abs(2*a3(n+1)-3*a2(n+1)+a(n+1)).*nu.^n);

Z1 = (L*nu)^2/(N*(N+1))*(6*norm_a^2+6*norm_a+1);
Z2 = 6*(L*nu)^2/(N*(N+1))*(1+2*norm_a);
Z3 = 6*(L*nu)^2/(N*(N+1));

fprintf('\n')
display(['Z3 = ',num2str(Z3),', Z2 = ',num2str(Z2),', Z1 = ',num2str(Z1),', Y0 = ',num2str(Y0)])
fprintf('\n')

p = [Z3 Z2 -(1-Z1) Y0];

I = sort(roots(p));

if norm(imag(I))>0 || Z1>1
    disp('failure')
    I=[-1 -1];
else
    I(1)=[];
    plot_a_taylor(a,x0,L,nu);
end

end

