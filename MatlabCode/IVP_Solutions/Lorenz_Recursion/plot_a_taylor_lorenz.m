function [] = plot_a_taylor_lorenz(a1,a2,a3,x0,L,nu)

%close all

N=length(a1)-1;

t=(-nu:.001:nu);

m=length(t);
u1=zeros(1,m);
u2=zeros(1,m);
u3=zeros(1,m);

for k=1:m
    u1(k)=sum(a1.*t(k).^(0:N)');
    u2(k)=sum(a2.*t(k).^(0:N)');
    u3(k)=sum(a3.*t(k).^(0:N)');
end

% plot(t*L,u1,'Linewidth',3,'color','blue')
% hold on
% plot(t*L,u2,'Linewidth',3,'color','green')
% plot(t*L,u3,'Linewidth',3,'color','red')
% set(gca,'FontSize',15)
% 
% plot(0,x0(1),'*','Linewidth',3,'color',[0 0 0])
% plot(0,x0(2),'*','Linewidth',3,'color',[0 0 0])
% plot(0,x0(3),'*','Linewidth',3,'color',[0 0 0])
% axis tight
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
%hold on
plot3(u1,u2,u3,'Linewidth',2,'color','blue')
plot3(x0(1),x0(2),x0(3),'*','Linewidth',1,'color',[0 0 0])
%options = odeset('RelTol',1e-6,'AbsTol',1e-6);
%[t,y]=ode45(@lorenz_classic,[0 50],[1;1;1],options);  
%plot3(y(:,1),y(:,2),y(:,3),'color','green')
% set(gca,'FontSize',15)
% axis tight
% xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 20)
% ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 20)
% zlabel('$$x_3$$', 'Interpreter', 'latex', 'FontSize', 20)
% 
% view(40,20);

end

    