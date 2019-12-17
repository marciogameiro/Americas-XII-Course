function [x]=padding_lorenz_comp(x,N)

%%% inputs: x = (omega,a1,a2,a3), 

m=(length(x)+2)/6; %%% N = m-1
omega=x(1);

a1=x(2:2*m);       %%% Complex Fourier coefficients of u1
a2=x(2*m+1:4*m-1); %%% Complex Fourier coefficients of u2
a3=x(4*m:6*m-2);

if N>0   
    a1=[zeros(N,1);a1;zeros(N,1)];
    a2=[zeros(N,1);a2;zeros(N,1)];
    a3=[zeros(N,1);a3;zeros(N,1)];   
else % if N<=0
    a1=a1(m-N:m+N);
    a2=a2(m-N:m+N);
    a3=a3(m-N:m+N);
end

x=[omega;a1;a2;a3];


end