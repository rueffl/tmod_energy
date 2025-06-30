function [e_loss,as_loss] = objective_loss(x, k_tr, Omega, vr, gamma, len)
N = round(x(1));
epsilon_kappa = x(2);
omega = x(3);
epsilon_rho = x(4);

phase_rho = 2*pi*rand(1,N); 
phase_kappa = 2*pi*rand(1,N); 

rs = []; ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))/2, 1, epsilon_rho*exp(1i*phase_rho(j))/2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))/2, 1, epsilon_kappa*exp(1i*phase_kappa(j))/2];
    rs = [rs; rs_j];
    ks = [ks; ks_j];
end

[~, ~, ~, ~, e_loss, as_loss] = get_Energy_regimes(k_tr, omega, Omega, rs, ks, vr, gamma, len, N);

end