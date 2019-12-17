function []=plot_periodic_complex(x)

close all

L=real(x(1));
period=2*pi/abs(L);

m=(length(x)+2)/6;

c1=x(2:2*m);       %%% Complex Fourier coefficients of u1
c2=x(2*m+1:4*m-1); %%% Complex Fourier coefficients of u2
c3=x(4*m:6*m-2);   %%% Complex Fourier coefficients of u3

a1=real(c1(m:end)); b1=imag(c1(m:end));
a2=real(c2(m:end)); b2=imag(c2(m:end));
a3=real(c3(m:end)); b3=imag(c3(m:end));

t=(0:.0001:period);
u1=periodic(L,a1,b1,t);
u2=periodic(L,a2,b2,t);
u3=periodic(L,a3,b3,t);

plot3(u1,u2, u3,'LineWidth',2)
view(20, 20)

figure
plot(t,u1,'Color',[0 0 1],'LineWidth',2)
box off
hold on
plot(t,u2,'Color',[1 0 0],'LineWidth',2)
plot(t,u3,'Color',[0 1 0],'LineWidth',2)

end

function [periodic]=periodic(L,a,b,t)

m=length(a);

sum=0;
for k=1:m-1
    sum=sum+a(k+1,1)*cos(k*L*t)-b(k+1,1)*sin(k*L*t);
end
periodic=a(1,1)+2*sum;

end