
%clear
close all

z0=27; % Fixed plane at height z=27

%%%% Data from the paper of D. Viswanath %%%%
% AB
%x0=-13.763610682134; y0=-19.578751942452; T=1.5586522107162;

%AAB
%x0=-12.595115397689; y0=-16.970525307084; T=2.3059072639398;

%ABB
%x0=-14.426408025035; y0=-21.111230056994; T=2.3059072639398; 

%AAAB
%x0=-11.998523280062; y0=-15.684254096883; T=3.0235837034339; 

%AABB
%x0=-12.915137970311; y0=-17.673100172646; T=3.0842767758221;

%ABBB
%x0=-14.839800671717; y0=-22.086807160089; T=3.0235837034339;

%AAAAB
%x0=-11.586951663722; y0=-14.814615288122; T=3.7256417715558;

%AAABB
%x0=-12.231406596273; y0=-16.182730598894; T=3.8202541634368;

%AABAB
x0=-12.698941349915; y0=-17.197497247713; T=3.8695391125646;

%AABBB
%x0=-13.056930146345; y0=-17.987214049281; T=3.8202541634368;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0=[x0;y0;z0];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

[t,y]=ode45(@lorenz_classic,[0 T],u0,options);   

%plot3(y(:,1),y(:,2),y(:,3));

% First approach based on computing the integrals %

%N=length(t)-1;
m=90; % # of Fourier coefficients we keep
omega=2*pi/T;
i=sqrt(-1);

a1=zeros(2*m-1,1);
a2=zeros(2*m-1,1);
a3=zeros(2*m-1,1);

for k=-m+1:-1
    v=exp(-i*k*omega*t(2:end))-exp(-i*k*omega*t(1:end-1));
    a1(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,1).*v);
    a2(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,2).*v);
    a3(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,3).*v);
end

a1(m)=omega/(2*pi)*sum(y(1:end-1,1).*(t(2:end)-t(1:end-1)));
a2(m)=omega/(2*pi)*sum(y(1:end-1,2).*(t(2:end)-t(1:end-1)));
a3(m)=omega/(2*pi)*sum(y(1:end-1,3).*(t(2:end)-t(1:end-1)));

for k=1:m-1
    v=exp(-i*k*omega*t(2:end))-exp(-i*k*omega*t(1:end-1));
    a1(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,1).*v);
    a2(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,2).*v);
    a3(k+m)=-1/(2*pi*k*i)*sum(y(1:end-1,3).*v);
end

% Second approach based on FFT %

% w1=ifft(y(:,1));
% w2=ifft(y(:,2));
% w3=ifft(y(:,3));
% 
% %%% The ordering of the frequencies is given by 
% %%%  [0 1 2 ... N -N+1 -N+2 ... -1], if length(w1) = 2N
% %%%  [0 1 2 ... N-1 -N+1 -N+2 ... -1], if length(w1) = 2N-1
% 
% m=10; % 2m-1 is the number of Fourier coefficients we keep
% 
% a1=[w1(end-m+2:end);w1(1:m)];
% a2=[w2(end-m+2:end);w2(1:m)];
% a3=[w3(end-m+2:end);w3(1:m)];
% omega=2*pi/T;

gamma0_phase_bar=u0;
gamma0_phase_dot=lorenz_classic(1,u0)/norm(lorenz_classic(1,u0));

gamma0_phase(:,1)=gamma0_phase_bar;
gamma0_phase(:,2)=gamma0_phase_dot;

clear gamma0_phase_bar gamma0_phase_dot y t

sigma=10; rho=28; beta=8/3;

par=[sigma;rho;beta];

k0=10;

clear w1 w2 w3 m options T x0 y0 z0 u0 beta sigma rho

x=[omega;a1;a2;a3];

clear N a1 a2 a3 i k omega v

[x]=newton_f_comp_lorenz(x,gamma0_phase,par,k0);

while norm(x(end-4:end))>1e-14
    [x]=padding_lorenz_comp(x,10);
    [x]=newton_f_comp_lorenz(x,gamma0_phase,par,k0);
end

x(1)=real(x(1));
    
nu=1.01;

 extra_dim=10;
 x=padding_lorenz_comp(x,extra_dim);

[I,success,Ranalytic]=rad_poly_analytic_lorenz(x,gamma0_phase,par,k0,nu);

if success==1
    plot_periodic_complex(x)
end

%save AABBB x gamma0_phase par k0 I nu
