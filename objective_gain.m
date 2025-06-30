function [e_gain,as_gain] = objective_gain(x, k_tr, Omega, vr, gamma, len)
% x = [N, epsilon_kappa, omega, epsilon_rho]
N = round(x(1)); % N must be integer
epsilon_kappa = x(2);
omega = x(3);
epsilon_rho = x(4);

% Create random phases (or define them externally)
phase_rho = 2*pi*rand(1,N); 
phase_kappa = 2*pi*rand(1,N); 

% Define rs and ks
rs = []; ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))/2, 1, epsilon_rho*exp(1i*phase_rho(j))/2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))/2, 1, epsilon_kappa*exp(1i*phase_kappa(j))/2];
    rs = [rs; rs_j];
    ks = [ks; ks_j];
end

% Call main function
[e_gain, as_gain, ~, ~, ~, ~] = get_Energy_regimes(k_tr, omega, Omega, rs, ks, vr, gamma, len, N);

% Negative because we want to maximize
e_gain = -e_gain;
end
