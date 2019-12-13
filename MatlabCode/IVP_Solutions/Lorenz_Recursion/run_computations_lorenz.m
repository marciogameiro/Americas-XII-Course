
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y]=ode45(@lorenz_classic,[0 50],[1;1;1],options);

figure; hold on
plot3(y(:,1),y(:,2),y(:,3),'color','green')
set(gca,'FontSize',10)
axis tight
xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 10)
zlabel('$$x_3$$', 'Interpreter', 'latex', 'FontSize', 10)

view(40,20);

sigma = 10; rho = 28; beta = 8/3;

L = 0.2; nu = 0.9; N=200;

par = [sigma; rho; beta];

M = length(t);
m = ceil(M/20);
index = 1:m:M;
index = index(2:end-1);
n = length(index);

for k=1:n
   x0=[y(index(k),1);y(index(k),2);y(index(k),3)];
   [I] = validation_ivp_lorenz(x0,L,nu,N,par);
end
