function [f] = f_comp_lorenz(x,gamma0_phase,par,k0)

% input: x = [omega;a1;a2;a3], where for i=1,2,3, 
%        ai=(ai_{-N},...,ai_{-1},ai_0,ai_1,...,ai_{N}) \in C^{2N+1}
%        k0 : 2*k0+1 is the number of Fourier modes we keep to evaluate the
%             phase condition 

% output  f = [f0;f1;f2;f3], where 
%        f0 = phase condition, and for i=1,2,3
%        fi = (fi_{-N},...,fi_{-1},fi_0,fi_1,...,fi_{N}) \in C^{2N+1}
% f : R x C^(6N+3) --> R x C^(6N+3).

sigma=par(1);
rho=par(2);
beta=par(3);

gamma0_phase_bar=gamma0_phase(:,1);
gamma0_phase_dot=gamma0_phase(:,2);

m=(length(x)+2)/6; %%% N = m-1
omega=x(1);

a1=x(2:2*m);       %%% Complex Fourier coefficients of u1
a2=x(2*m+1:4*m-1); %%% Complex Fourier coefficients of u2
a3=x(4*m:6*m-2);   %%% Complex Fourier coefficients of u3

conv12=quadratic_sum_complex(a1,a2);
conv13=quadratic_sum_complex(a1,a3);

%%% Phase condition
gamma_0=[sum(a1(m-k0:m+k0));sum(a2(m-k0:m+k0));sum(a3(m-k0:m+k0))];
phase=gamma0_phase_dot'*(gamma0_phase_bar-gamma_0);

%%% Linear part
ikomega=1i*(-m+1:m-1)'*omega;

%%% Evaluating f
f=[phase;-ikomega.*a1+sigma*(a2-a1);-ikomega.*a2+(rho*a1-a2-conv13);-ikomega.*a3+(-beta*a3+conv12)];

end
