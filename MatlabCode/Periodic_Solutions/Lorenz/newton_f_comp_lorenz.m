function [x]=newton_f_comp_lorenz(x,gamma0_phase,par,k0)

tol=1e-12;

[f] = f_comp_lorenz(x,gamma0_phase,par,k0);

display(['At the beggining ||f|| = ',num2str(norm(f))])
Df_inv = inv(Df_comp_lorenz(x,gamma0_phase,par,k0));

k=0;
while (k<=40) && (norm(f)> tol),
    x = x - Df_inv * f_comp_lorenz(x,gamma0_phase,par,k0);
    Df_inv = inv(Df_comp_lorenz(x,gamma0_phase,par,k0));
    f = f_comp_lorenz(x,gamma0_phase,par,k0);
    display(['||f|| = ',num2str(norm(f))])
    k=k+1;
end

%plot_periodic_complex(x)

return
