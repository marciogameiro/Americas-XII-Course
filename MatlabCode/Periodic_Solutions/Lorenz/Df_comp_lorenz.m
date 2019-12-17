function [Df] = Df_comp_lorenz(x,gamma0_phase,par,k0)

% input: x = [omega,a1,a2,a3], where for i=1,2,3, 
%  ai=(ai_{-N},...,ai_{-1},ai_0,ai_1,...,ai_{N}) \in C^{2N+1}
%
% output Df = [ df_0/d_omega df_0/da_1 df_0/da_2 df_0/da_3 ;
%               df_1/d_omega df_1/da_1 df_1/da_2 df_1/da_3 ;
%               df_2/d_omega df_2/da_1 df_2/da_2 df_2/da_3 ;
%               df_3/d_omega df_3/da_1 df_3/da_2 df_3/da_3  ]
%
%  where for i=1,2,3,
%       fi=(fi_{-N},...,fi_{-1},fi_0,fi_1,...,fi_{N}) \in C^{2N+1}
    
sigma=par(1);
rho=par(2);
beta=par(3);

gamma0_phase_dot=gamma0_phase(:,2);

m=(length(x)+2)/6; %%% N = m-1

Df = zeros(6*m-2);

omega=x(1);        %%% Frequency variable
a1=x(2:2*m);       %%% Complex Fourier coefficients of u1
a2=x(2*m+1:4*m-1); %%% Complex Fourier coefficients of u2
a3=x(4*m:6*m-2);   %%% Complex Fourier coefficients of u3

ikomega=1i*(-m+1:m-1)'*omega;
ik=1i*(-m+1:m-1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Derivatives of the quadratic terms in df_i/da_j %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

da1a3_da1=zeros(2*m-1);
da1a3_da3=zeros(2*m-1);
da1a2_da1=zeros(2*m-1);
da1a2_da2=zeros(2*m-1);

a1_extended=[zeros(1,m-1) reshape(a1,1,2*m-1) zeros(1,m-1)];
a2_extended=[zeros(1,m-1) reshape(a2,1,2*m-1) zeros(1,m-1)];
a3_extended=[zeros(1,m-1) reshape(a3,1,2*m-1) zeros(1,m-1)];

for k=-m+1:m-1
    da1a3_da1(k+m,:)=a3_extended(k+2*m-1-(-m+1:m-1));
    da1a3_da3(k+m,:)=a1_extended(k+2*m-1-(-m+1:m-1));
    da1a2_da1(k+m,:)=a2_extended(k+2*m-1-(-m+1:m-1));
    da1a2_da2(k+m,:)=a1_extended(k+2*m-1-(-m+1:m-1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_0/da_i, i=1,2,3 [First row] %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ones_k0=ones(1,2*k0+1);

%%% df_0/da_1 %%%
Df(1,m-k0+1:m+k0+1)=-gamma0_phase_dot(1)*ones_k0;

%%% df_0/da_2 %%%
Df(1,3*m-k0:3*m+k0)=-gamma0_phase_dot(2)*ones_k0;

%%% df_0/da_3 %%%
Df(1,5*m-k0-1:5*m+k0-1)=-gamma0_phase_dot(3)*ones_k0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_i/d_omega, i=1,2,3 [First column] %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Df(2:2*m,1)=-ik.*a1;
Df(2*m+1:4*m-1,1)=-ik.*a2;
Df(4*m:6*m-2,1)=-ik.*a3;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_1/da_1 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(2:2*m,2:2*m)=diag(-ikomega-sigma);

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_1/da_2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(2:2*m,2*m+1:4*m-1)=diag(sigma*ones(1,2*m-1));

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_2/da_1 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(2*m+1:4*m-1,2:2*m)=diag(rho*ones(1,2*m-1))-da1a3_da1;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_2/da_2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(2*m+1:4*m-1,2*m+1:4*m-1)=diag(-ikomega-1);

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_2/da_3 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(2*m+1:4*m-1,4*m:6*m-2)=-da1a3_da3;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_3/da_1 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(4*m:6*m-2,2:2*m)=da1a2_da1;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_3/da_2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(4*m:6*m-2,2*m+1:4*m-1)=da1a2_da2;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% df_3/da_3 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Df(4*m:6*m-2,4*m:6*m-2)=diag(-ikomega-beta);

end
