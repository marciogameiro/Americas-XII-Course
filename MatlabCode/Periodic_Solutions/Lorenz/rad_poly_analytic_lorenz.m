function [I,success,Ranalytic]=rad_poly_analytic_lorenz(x,gamma0_phase,par,k0,nu)

sigma=par(1);
rho=par(2);
beta=par(3);

display(['nu = ',num2str(nu)])

m=(length(x)+2)/6; %%% N = m-1
omega=x(1);

a1=x(2:2*m);       %%% Complex Fourier coefficients of u1
a2=x(2*m+1:4*m-1); %%% Complex Fourier coefficients of u2
a3=x(4*m:6*m-2);   %%% Complex Fourier coefficients of u3

nu_power = nu.^abs(-m+1:m-1);
a1_nu_norm=sum(abs(a1).*nu_power');
a2_nu_norm=sum(abs(a2).*nu_power');
a3_nu_norm=sum(abs(a3).*nu_power');
%display(['The nu-norm of a1 is = ',num2str(a1_nu_norm)])
%display(['The nu-norm of a2 is = ',num2str(a2_nu_norm)])
%display(['The nu-norm of a3 is = ',num2str(a3_nu_norm)])

Df = Df_comp_lorenz(x,gamma0_phase,par,k0);
Am=inv(Df);
%display(['norm of Am is = ',num2str(norm(Am))])

%%%%%%%%%
%%% Y %%%
%%%%%%%%%

%%% The use of m1 should be used to compute the f and the norms

a1_extended=[zeros(m-1,1);a1;zeros(m-1,1)];
a2_extended=[zeros(m-1,1);a2;zeros(m-1,1)];
a3_extended=[zeros(m-1,1);a3;zeros(m-1,1)];

x_extended = [omega;a1_extended;a2_extended;a3_extended];
f = f_comp_lorenz(x_extended,gamma0_phase,par,k0); 
fN = [f(1);f(m+1:3*m-1);f(5*m-2:7*m-4);f(9*m-5:11*m-7)];

yN=Am*fN;

yN1=yN(2:2*m);       
yN2=yN(2*m+1:4*m-1); 
yN3=yN(4*m:6*m-2);

Y0 = abs(yN(1));

yN1_nu_norm=sum(abs(yN1).*nu_power');
yN2_nu_norm=sum(abs(yN2).*nu_power');
yN3_nu_norm=sum(abs(yN3).*nu_power');

nu_power_extended = [nu.^abs(-2*m+2:-m)';nu.^abs(m:2*m-2)'];
omega_k_extended = omega*[abs(-2*m+2:-m)';abs(m:2*m-2)'];
y_ext1_nu_norm=sum(abs([f(2:m);f(3*m:4*m-2)]).*nu_power_extended./omega_k_extended);
y_ext2_nu_norm=sum(abs([f(4*m-1:5*m-3);f(7*m-3:8*m-5)]).*nu_power_extended./omega_k_extended);
y_ext3_nu_norm=sum(abs([f(8*m-4:9*m-6);f(11*m-6:12*m-8)]).*nu_power_extended./omega_k_extended);

Y1 = yN1_nu_norm + y_ext1_nu_norm;
Y2 = yN2_nu_norm + y_ext2_nu_norm;
Y3 = yN3_nu_norm + y_ext3_nu_norm;

%display(['The nu-norm of Y0 is = ',num2str(Y0)])
%display(['The nu-norm of Y1 is = ',num2str(Y1)])
%display(['The nu-norm of Y2 is = ',num2str(Y2)])
%display(['The nu-norm of Y3 is = ',num2str(Y3)])

Y=[Y0;Y1;Y2;Y3];

%%%%%%%%%%%%%%%%%%%%%%%
%%% Z0, Z1, Z2 & Z3 %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Z0 %%%

B=abs(eye(6*m-2)-Am*Df);

Z0=zeros(4,1);

B_omega_0 = B(1,1);
B_a1_0 = max(B(1,2:2*m)./nu_power);
B_a2_0 = max(B(1,2*m+1:4*m-1)./nu_power);
B_a3_0 = max(B(1,4*m:6*m-2)./nu_power);

Z0(1) = B_omega_0 + B_a1_0 + B_a2_0 + B_a3_0;

B_omega_1 = sum(B(2:2*m,1).*nu_power');
B_omega_2 = sum(B(2*m+1:4*m-1,1).*nu_power');
B_omega_3 = sum(B(4*m:6*m-2,1).*nu_power'); 

B11 = B(2:2*m,2:2*m);
B21 = B(2:2*m,2*m+1:4*m-1);
B31 = B(2:2*m,4*m:6*m-2);
B12 = B(2*m+1:4*m-1,2:2*m);
B22 = B(2*m+1:4*m-1,2*m+1:4*m-1);
B32 = B(2*m+1:4*m-1,4*m:6*m-2);
B13 = B(4*m:6*m-2,2:2*m);
B23 = B(4*m:6*m-2,2*m+1:4*m-1);
B33 = B(4*m:6*m-2,4*m:6*m-2);

B_a1_1 = max(sum(B11.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a2_1 = max(sum(B21.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a3_1 = max(sum(B31.*repmat(nu_power,2*m-1,1)')./nu_power);

B_a1_2 = max(sum(B12.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a2_2 = max(sum(B22.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a3_2 = max(sum(B32.*repmat(nu_power,2*m-1,1)')./nu_power);

B_a1_3 = max(sum(B13.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a2_3 = max(sum(B23.*repmat(nu_power,2*m-1,1)')./nu_power);
B_a3_3 = max(sum(B33.*repmat(nu_power,2*m-1,1)')./nu_power);

Z0(2) = B_omega_1 + B_a1_1 + B_a2_1 + B_a3_1;
Z0(3) = B_omega_2 + B_a1_2 + B_a2_2 + B_a3_2;
Z0(4) = B_omega_3 + B_a1_3 + B_a2_3 + B_a3_3;

%display(['The nu-norm of Z0_1 is = ',num2str(Z0(1))])
%display(['The nu-norm of Z0_2 is = ',num2str(Z0(2))])
%display(['The nu-norm of Z0_3 is = ',num2str(Z0(3))])
%display(['The nu-norm of Z0_4 is = ',num2str(Z0(4))])

%%% Z1 %%%

Z1=zeros(4,1);

A20=abs(Am(1,2*m+1:4*m-1));
A30=abs(Am(1,4*m:6*m-2));

A11 = abs(Am(2:2*m,2:2*m));
A21 = abs(Am(2:2*m,2*m+1:4*m-1));
A31 = abs(Am(2:2*m,4*m:6*m-2));
A12 = abs(Am(2*m+1:4*m-1,2:2*m));
A22 = abs(Am(2*m+1:4*m-1,2*m+1:4*m-1));
A32 = abs(Am(2*m+1:4*m-1,4*m:6*m-2));
A13 = abs(Am(4*m:6*m-2,2:2*m));
A23 = abs(Am(4*m:6*m-2,2*m+1:4*m-1));
A33 = abs(Am(4*m:6*m-2,4*m:6*m-2));

ell_hat_a1=zeros(2*m-1,1);
ell_hat_a2=zeros(2*m-1,1);
ell_hat_a3=zeros(2*m-1,1);

for k=-m+1:-1
    ell_hat_a1(k+m)=max(nu.^-(m-1+abs((k:-1))').*abs(a1(k-(k:-1)+(m-1)+m)));
    ell_hat_a2(k+m)=max(nu.^-(m-1+abs((k:-1))').*abs(a2(k-(k:-1)+(m-1)+m)));
    ell_hat_a3(k+m)=max(nu.^-(m-1+abs((k:-1))').*abs(a3(k-(k:-1)+(m-1)+m)));
end

for k=1:m-1
    ell_hat_a1(k+m)=max(nu.^-(m-1+(1:k)').*abs(a1(k-(1:k)-(m-1)+m)));
    ell_hat_a2(k+m)=max(nu.^-(m-1+(1:k)').*abs(a2(k-(1:k)-(m-1)+m)));
    ell_hat_a3(k+m)=max(nu.^-(m-1+(1:k)').*abs(a3(k-(1:k)-(m-1)+m)));
end

beta3_a2_0=sum((A20*ell_hat_a3).*nu_power');
beta2_a3_0=sum((A30*ell_hat_a2).*nu_power');
beta1_a3_0=sum((A30*ell_hat_a1).*nu_power');
beta1_a2_0=sum((A20*ell_hat_a1).*nu_power');

beta3_a2_2=sum((A22*ell_hat_a3).*nu_power');
beta2_a3_2=sum((A32*ell_hat_a2).*nu_power');
beta3_a2_1=sum((A21*ell_hat_a3).*nu_power');
beta2_a3_1=sum((A31*ell_hat_a2).*nu_power');
beta2_a3_3=sum((A33*ell_hat_a2).*nu_power');
beta3_a2_3=sum((A23*ell_hat_a3).*nu_power');
beta1_a3_1=sum((A31*ell_hat_a1).*nu_power');
beta1_a3_2=sum((A32*ell_hat_a1).*nu_power');
beta1_a3_3=sum((A33*ell_hat_a1).*nu_power');
beta1_a2_1=sum((A23*ell_hat_a1).*nu_power');
beta1_a2_2=sum((A22*ell_hat_a1).*nu_power');
beta1_a2_3=sum((A23*ell_hat_a1).*nu_power');

Z1(1)=beta3_a2_0+beta2_a3_0+beta1_a3_0+beta1_a2_0;
Z1(2)=2*abs(sigma)/(omega*m)+beta3_a2_1+beta2_a3_1+beta1_a3_1+beta1_a2_1;
Z1(3)=(abs(rho)+1+a3_nu_norm+a1_nu_norm)/(omega*m)+beta3_a2_2+beta2_a3_2+beta1_a3_2+beta1_a2_2;
Z1(4)=(abs(beta)+a1_nu_norm+a2_nu_norm)/(omega*m)+beta3_a2_3+beta2_a3_3+beta1_a3_3+beta1_a2_3;

%display(['The nu-norm of Z1_1 is = ',num2str(Z1(1))])
%display(['The nu-norm of Z1_2 is = ',num2str(Z1(2))])
%display(['The nu-norm of Z1_3 is = ',num2str(Z1(3))])
%display(['The nu-norm of Z1_4 is = ',num2str(Z1(4))])

%%% Z2 %%%

Z2=zeros(4,1);

hat_B_a1_0 = max(abs(Am(1,2:2*m))./nu_power.*abs(-m+1:m-1));
hat_B_a2_0 = max(abs(Am(1,2*m+1:4*m-1))./nu_power.*abs(-m+1:m-1));
hat_B_a3_0 = max(abs(Am(1,4*m:6*m-2))./nu_power.*abs(-m+1:m-1));

hat_B_a1_1 = max(sum(A11.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a2_1 = max(sum(A21.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a3_1 = max(sum(A31.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));

hat_B_a1_2 = max(sum(A12.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a2_2 = max(sum(A22.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a3_2 = max(sum(A32.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));

hat_B_a1_3 = max(sum(A13.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a2_3 = max(sum(A23.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));
hat_B_a3_3 = max(sum(A33.*repmat(nu_power,2*m-1,1)')./nu_power.*abs(-m+1:m-1));

Z2(1)=2*(hat_B_a1_0+hat_B_a2_0+hat_B_a3_0+B_a2_0+B_a3_0);
Z2(2)=2*(max(hat_B_a1_1,1/omega)+hat_B_a2_1+hat_B_a3_1+B_a3_1+B_a2_1);
Z2(3)=2*(max(hat_B_a2_2,1/omega)+hat_B_a1_2+hat_B_a3_2+max(B_a2_2,1/(omega*m))+B_a3_2);
Z2(4)=2*(max(hat_B_a3_3,1/omega)+hat_B_a1_3+hat_B_a2_3+max(B_a3_3,1/(omega*m))+B_a2_3);

%display(['The nu-norm of Z2_1 is = ',num2str(Z2(1))])
%display(['The nu-norm of Z2_2 is = ',num2str(Z2(2))])
%display(['The nu-norm of Z2_3 is = ',num2str(Z2(3))])
%display(['The nu-norm of Z2_4 is = ',num2str(Z2(4))])

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Radii polynomial %%%
%%%%%%%%%%%%%%%%%%%%%%%%

p=[Z2 Z0+Z1-1 Y];

I=zeros(4,2);
for k=1:4
    R=roots(p(k,:));
    R=sort(R);
    Int=R';
    I(k,:)=Int;
end

I=[max(I(:,1)) min(I(:,2))];
Ranalytic=0;
if norm(imag(I))>0
    display('FAILURE !')
    I=[-1 1];
    success=0;
elseif I(1)<0
    display('FAILURE !')
    I=[-1 1];
    success=0;
else    
    Ranalytic=log(nu)/omega;
    display(['SUCCESS ! The radius of analyticity is = ',num2str(Ranalytic)])
    display(['I = ',num2str(I)])
    success=1;
end

end