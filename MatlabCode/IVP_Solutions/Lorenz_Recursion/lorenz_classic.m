function du = lorenz_classic(t,u)

sigma=10; rho=28; beta=8/3;

du=[
    sigma*(u(2)-u(1));
    rho*u(1)-u(2)-u(1)*u(3);
    -beta*u(3)+u(1)*u(2)
    ];