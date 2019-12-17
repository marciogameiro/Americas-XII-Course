function [DF]=finite_diff_f_comp_lorenz(x,gamma0_phase,par,k0)

h=1e-6;
m=length(x);
E=eye(m);
DF=zeros(m);
for j=1:m
    xh=x+h*E(:,j);
    DF(:,j)=(f_comp_lorenz(xh,gamma0_phase,par,k0)-f_comp_lorenz(x,gamma0_phase,par,k0))/h;
end
end