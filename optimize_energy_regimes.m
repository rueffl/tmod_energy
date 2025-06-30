% Main script for optimizing e_gain and e_loss
% N is iterated externally.
% epsilon_kappa, epsilon_rho, and omega (real & imag) are optimized by fmincon.
% Uses results from N-1 as initial guess for N ("warm start").
format long

% --- CONFIGURATION ---
% Set range of N values to test
N_values_to_test = 2:2:100; % Example: Test N from 3 to 7. Adjust as needed.

% Initial guesses for the first N optimization
% These will be updated for subsequent N values (warm start)
epsilon_kappa_init_guess = 0.5; % Initial guess for epsilon_kappa
epsilon_rho_init_guess = 0.5; % Initial guess for epsilon_rho
omega_init_guess = 0.01 + 0i; % Initial guess for omega
gamma_init_guess = 0.05; % Initial guess for gamma

% Optimization options for fmincon
options = optimoptions('fmincon', ...
                       'Display', 'iter', ...
                       'Algorithm', 'sqp', ...
                       'MaxIterations', 100, ...
                       'StepTolerance', 1e-8, ...
                       'OptimalityTolerance', 1e-8, ...
                       'ConstraintTolerance', 1e-10); % Added constraint tolerance

% Define bounds for the optimization variables: [epsilon_kappa, epsilon_rho, real(omega), imag(omega)]
params_lower_bounds = [0,   0,   0,  0, 0];  % Lower bounds for eps_k, eps_rho, real(om), imag(om)
params_upper_bounds = [0.99,   0.99,   0.08, 0, 50]; % Upper bounds for eps_k, eps_rho, real(om), imag(om)
                          % Ensure epsilon_rho <= 1 (or other physically meaningful upper bound)

% --- OPTIMIZATION FOR MAXIMUM e_gain ---
disp('Optimizing for maximum e_gain...');
disp('Iterating through N values with warm starts...');

all_results_gain = struct('N', num2cell(N_values_to_test), ...
                          'epsilon_kappa', NaN, ...
                          'epsilon_rho', NaN, ...
                          'omega', NaN+1i*NaN, ...
                          'gamma', NaN, ...
                          'fval', Inf, ... % fval stores -e_gain
                          'as_val', []);
best_fval_gain_overall = Inf;
best_result_gain_idx = 0;

% Initial guess for the continuous parameters for the very first N
current_initial_guess_params_gain = [epsilon_kappa_init_guess, epsilon_rho_init_guess, real(omega_init_guess), imag(omega_init_guess), gamma_init_guess];

for i = 1:length(N_values_to_test)
    N_current = N_values_to_test(i);
    fprintf('--- Optimizing e_gain for N = %d ---\n', N_current);
    fprintf('Initial guess for [eps_k, eps_rho, re(om), im(om), gamma]: [%.2f, %.2f, %.2f, %.2f, %.2f]\n', ...
            current_initial_guess_params_gain(1), current_initial_guess_params_gain(2), ...
            current_initial_guess_params_gain(3), current_initial_guess_params_gain(4), current_initial_guess_params_gain(5));

    [opt_params_subset_gain, fval_current_gain, as_current_gain] = ...
        optimize_params_for_gain(N_current, current_initial_guess_params_gain, ...
                                 params_lower_bounds, params_upper_bounds, options);

    all_results_gain(i).epsilon_kappa = opt_params_subset_gain(1);
    all_results_gain(i).epsilon_rho   = opt_params_subset_gain(2);
    all_results_gain(i).omega         = opt_params_subset_gain(3) + 1i*opt_params_subset_gain(4);
    all_results_gain(i).gamma         = opt_params_subset_gain(5);
    all_results_gain(i).fval          = fval_current_gain; 
    all_results_gain(i).as_val        = as_current_gain;

    if fval_current_gain < best_fval_gain_overall
        best_fval_gain_overall = fval_current_gain;
        best_result_gain_idx = i;
    end

    % Update initial guess for the next N (warm start)
    current_initial_guess_params_gain = opt_params_subset_gain;
end

if best_result_gain_idx > 0
    optimized_N_gain = all_results_gain(best_result_gain_idx).N;
    optimized_epsilon_kappa_gain = all_results_gain(best_result_gain_idx).epsilon_kappa;
    optimized_epsilon_rho_gain = all_results_gain(best_result_gain_idx).epsilon_rho;
    optimized_omega_gain_complex = all_results_gain(best_result_gain_idx).omega;
    optimized_gamma_gain = all_results_gain(best_result_gain_idx).gamma;
    final_fval_gain = all_results_gain(best_result_gain_idx).fval;
    final_as_gain = all_results_gain(best_result_gain_idx).as_val;

    disp(sprintf('\nOverall Optimized Parameters for e_gain:'));
    disp(['Optimized N: ', num2str(optimized_N_gain)]);
    disp(['Optimized epsilon_kappa: ', num2str(optimized_epsilon_kappa_gain)]);
    disp(['Optimized epsilon_rho: ', num2str(optimized_epsilon_rho_gain)]);
    disp(['Optimized omega: ', num2str(optimized_omega_gain_complex)]);
    disp(['Optimized gamma: ', num2str(optimized_gamma_gain)]);
    disp(['Maximized e_gain value: ', num2str(-final_fval_gain)]); % -fval is e_gain
    disp(['Corresponding incident wave field (as_gain): ', num2str(final_as_gain)]);
else
    disp('No successful optimization for e_gain found.');
end

% --- OPTIMIZATION FOR MINIMUM e_loss ---
disp(sprintf('\nOptimizing for minimum e_loss...'));
disp('Iterating through N values with warm starts...');

all_results_loss = struct('N', num2cell(N_values_to_test), ...
                          'epsilon_kappa', NaN, ...
                          'epsilon_rho', NaN, ...
                          'omega', NaN+1i*NaN, ...
                          'gamma', NaN, ...
                          'fval', Inf, ... % fval stores e_loss
                          'as_val', []);
best_fval_loss_overall = Inf;
best_result_loss_idx = 0;

% Initial guess for the continuous parameters for the very first N
current_initial_guess_params_loss = [epsilon_kappa_init_guess, epsilon_rho_init_guess, real(omega_init_guess), imag(omega_init_guess), gamma_init_guess];

for i = 1:length(N_values_to_test)
    N_current = N_values_to_test(i);
    fprintf('--- Optimizing e_loss for N = %d ---\n', N_current);
    fprintf('Initial guess for [eps_k, eps_rho, re(om), im(om), gamma]: [%.2f, %.2f, %.2f, %.2f, %.2f]\n', ...
            current_initial_guess_params_loss(1), current_initial_guess_params_loss(2), ...
            current_initial_guess_params_loss(3), current_initial_guess_params_loss(4), current_initial_guess_params_loss(5));

    [opt_params_subset_loss, fval_current_loss, as_current_loss] = ...
        optimize_params_for_loss(N_current, current_initial_guess_params_loss, ...
                                 params_lower_bounds, params_upper_bounds, options);

    all_results_loss(i).epsilon_kappa = opt_params_subset_loss(1);
    all_results_loss(i).epsilon_rho   = opt_params_subset_loss(2);
    all_results_loss(i).omega         = opt_params_subset_loss(3) + 1i*opt_params_subset_loss(4);
    all_results_loss(i).gamma         = opt_params_subset_loss(5);
    all_results_loss(i).fval          = fval_current_loss; % This is e_loss
    all_results_loss(i).as_val        = as_current_loss;

    if fval_current_loss < best_fval_loss_overall
        best_fval_loss_overall = fval_current_loss;
        best_result_loss_idx = i;
    end
    
    % Update initial guess for the next N (warm start)
    current_initial_guess_params_loss = opt_params_subset_loss;
end

if best_result_loss_idx > 0
    optimized_N_loss = all_results_loss(best_result_loss_idx).N;
    optimized_epsilon_kappa_loss = all_results_loss(best_result_loss_idx).epsilon_kappa;
    optimized_epsilon_rho_loss = all_results_loss(best_result_loss_idx).epsilon_rho;
    optimized_omega_loss_complex = all_results_loss(best_result_loss_idx).omega;
    optimized_gamma_loss = all_results_loss(best_result_loss_idx).gamma;
    final_fval_loss = all_results_loss(best_result_loss_idx).fval;
    final_as_loss = all_results_loss(best_result_loss_idx).as_val;

    disp(sprintf('\nOverall Optimized Parameters for e_loss:'));
    disp(['Optimized N: ', num2str(optimized_N_loss)]);
    disp(['Optimized epsilon_kappa: ', num2str(optimized_epsilon_kappa_loss)]);
    disp(['Optimized epsilon_rho: ', num2str(optimized_epsilon_rho_loss)]);
    disp(['Optimized omega: ', num2str(optimized_omega_loss_complex)]);
    disp(['Optimized gamma: ', num2str(optimized_gamma_loss)]);
    disp(['Minimized e_loss value: ', num2str(final_fval_loss)]);
    disp(['Corresponding incident wave field (as_loss): ', num2str(final_as_loss)]);
else
    disp('No successful optimization for e_loss found.');
end

%% Create Plots

fval_gain = -1.*[all_results_gain.fval];
fval_loss = [all_results_loss.fval];
omega_gain = [all_results_gain.omega];
omega_loss = [all_results_loss.omega];
gamma_gain = [all_results_gain.gamma];
gamma_loss = [all_results_loss.gamma];

figure()
hold on
% subplot(1,2,1)
% plot(N_values_to_test,fval_gain,'*-g')
% grid on
% xlabel('$N$','FontSize',18,'Interpreter','latex')
% ylabel('$E_{\mathrm{gain}}$','FontSize',18,'Interpreter','latex')
% subplot(1,2,2)
plot(N_values_to_test,fval_loss,'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$E_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')

figure()
hold on
% subplot(1,2,1)
% plot(N_values_to_test,real(omega_gain),'*-g')
% grid on
% xlabel('$N$','FontSize',18,'Interpreter','latex')
% ylabel('$\omega_{\mathrm{gain}}$','FontSize',18,'Interpreter','latex')
% subplot(1,2,2)
plot(N_values_to_test,real(omega_loss),'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$\omega_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')

figure()
hold on
% subplot(1,2,1)
% plot(N_values_to_test,gamma_gain,'*-g')
% grid on
% xlabel('$N$','FontSize',18,'Interpreter','latex')
% ylabel('$\gamma_{\mathrm{gain}}$','FontSize',18,'Interpreter','latex')
% subplot(1,2,2)
plot(N_values_to_test,gamma_loss,'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$\gamma_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')

figure()
subplot(1,3,1)
hold on
plot(N_values_to_test,fval_loss,'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$E_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')
ylim([0.79,0.81])
subplot(1,3,2)
hold on
plot(N_values_to_test,real(omega_loss),'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$\omega_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')
subplot(1,3,3)
hold on
plot(N_values_to_test,gamma_loss,'*-r')
grid on
xlabel('$N$','FontSize',18,'Interpreter','latex')
ylabel('$\gamma_{\mathrm{loss}}$','FontSize',18,'Interpreter','latex')

%% Functions


% --- Function to optimize [eps_k, eps_rho, omega] for e_gain for FIXED N using FMINCON ---
function [optimized_params_subset, fval, as_val] = optimize_params_for_gain(N_fixed, initial_guess_subset, lb, ub, options_in)
    % initial_guess_subset = [epsilon_kappa, epsilon_rho, real(omega), imag(omega)]
    objective_fn = @(params_subset) energy_objective_gain_fixed_N_ek_er(params_subset, N_fixed);

    [optimized_params_subset, fval, exitflag, output] = fmincon(objective_fn, initial_guess_subset, ...
                                               [], [], [], [], lb, ub, [], options_in);
    if exitflag <= 0
        fprintf('Warning: fmincon for e_gain (N=%d) might not have converged. Exitflag: %d. Output: %s\n', N_fixed, exitflag, output.message);
    end

    epsilon_kappa_opt = optimized_params_subset(1);
    epsilon_rho_opt   = optimized_params_subset(2);
    omega_opt_complex = optimized_params_subset(3) + 1i * optimized_params_subset(4);
    gamma_opt         = optimized_params_subset(5);

    % Recalculate as_val with optimized parameters
    k_tr = 4; len = 0.01; xi = 0.2; Omega = xi * len; 
    phase_rho = ones(1, N_fixed).*pi/2; phase_kappa = ones(1, N_fixed).*pi/2;
    vr = 1; gamma = gamma_opt;
    rs = []; ks = [];
    for j = 1:N_fixed
        rs_j = [epsilon_rho_opt*exp(-1i*phase_rho(j))./2,1,epsilon_rho_opt*exp(1i*phase_rho(j))./2]; % Use optimized epsilon_rho
        ks_j = [epsilon_kappa_opt*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa_opt*exp(1i*phase_kappa(j))./2]; % Use optimized epsilon_kappa
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    [~, as_val, ~, ~, ~, ~] = get_Energy_regimes(k_tr, omega_opt_complex, Omega, rs, ks, vr, gamma, len, N_fixed);
end

% --- Function to optimize [eps_k, eps_rho, omega] for e_loss for FIXED N using FMINCON ---
function [optimized_params_subset, fval, as_val] = optimize_params_for_loss(N_fixed, initial_guess_subset, lb, ub, options_in)
    objective_fn = @(params_subset) energy_objective_loss_fixed_N_ek_er(params_subset, N_fixed);

    [optimized_params_subset, fval, exitflag, output] = fmincon(objective_fn, initial_guess_subset, ...
                                               [], [], [], [], lb, ub, [], options_in);
    if exitflag <= 0
        fprintf('Warning: fmincon for e_loss (N=%d) might not have converged. Exitflag: %d. Output: %s\n',N_fixed, exitflag, output.message);
    end
    
    epsilon_kappa_opt = optimized_params_subset(1);
    epsilon_rho_opt   = optimized_params_subset(2);
    omega_opt_complex = optimized_params_subset(3) + 1i * optimized_params_subset(4);
    gamma_opt         = optimized_params_subset(5);

    k_tr = 4; len = 0.01; xi = 0.2; Omega = xi * len; 
    phase_rho = ones(1, N_fixed).*pi/2; phase_kappa = ones(1, N_fixed).*pi/2;
    vr = 1; gamma = gamma_opt;
    rs = []; ks = [];
    for j = 1:N_fixed
        rs_j = [epsilon_rho_opt*exp(-1i*phase_rho(j))./2,1,epsilon_rho_opt*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa_opt*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa_opt*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    [~, ~, ~, ~, ~, as_val] = get_Energy_regimes(k_tr, omega_opt_complex, Omega, rs, ks, vr, gamma, len, N_fixed); % Get as_loss
end

% --- Objective function for maximizing e_gain (fixed N) ---
% params_subset = [epsilon_kappa, epsilon_rho, real(omega), imag(omega)]
function objective_val = energy_objective_gain_fixed_N_ek_er(params_subset, N_fixed)
    epsilon_kappa = params_subset(1);
    epsilon_rho   = params_subset(2);
    omega_complex = params_subset(3) + 1i * params_subset(4);
    gamma         = params_subset(5);

    large_penalty_value = 1e20; 
    % Bounds are handled by fmincon. This penalty is for NaN/Inf from get_Energy_regimes.

    k_tr = 4; len = 0.01; xi = 0.2; Omega = xi * len; 
    phase_rho_val = ones(1, N_fixed).*pi/2; % Renamed to avoid conflict with phase_rho in outer scope
    phase_kappa_val = ones(1, N_fixed).*pi/2; % Renamed
    vr = 1; 
    rs = []; ks = [];
    for j = 1:N_fixed
        % Use the epsilon_rho from params_subset
        rs_j = [epsilon_rho*exp(-1i*phase_rho_val(j))./2,1,epsilon_rho*exp(1i*phase_rho_val(j))./2];
        % Use the epsilon_kappa from params_subset
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa_val(j))./2,1,epsilon_kappa*exp(1i*phase_kappa_val(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    [e_gain, ~, ~, ~, ~, ~] = get_Energy_regimes(k_tr, omega_complex, Omega, rs, ks, vr, gamma, len, N_fixed);

    if isnan(e_gain)
        objective_val = large_penalty_value; return;
    end
    if isinf(e_gain)
        objective_val = -e_gain; return; % If e_gain=+Inf, obj=-Inf (good). If e_gain=-Inf, obj=+Inf (bad).
    end
    objective_val = -e_gain;
end

% --- Objective function for minimizing e_loss (fixed N) ---
% params_subset = [epsilon_kappa, epsilon_rho, real(omega), imag(omega)]
function objective_val = energy_objective_loss_fixed_N_ek_er(params_subset, N_fixed)
    epsilon_kappa = params_subset(1);
    epsilon_rho   = params_subset(2);
    omega_complex = params_subset(3) + 1i * params_subset(4);
    gamma         = params_subset(5);

    large_penalty_value = 1e20;

    k_tr = 4; len = 0.01; xi = 0.2; Omega = xi * len;
    phase_rho_val = ones(1, N_fixed).*pi/2; % Renamed
    phase_kappa_val = ones(1, N_fixed).*pi/2; % Renamed
    vr = 1; 
    rs = []; ks = [];
    for j = 1:N_fixed
        rs_j = [epsilon_rho*exp(-1i*phase_rho_val(j))./2,1,epsilon_rho*exp(1i*phase_rho_val(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa_val(j))./2,1,epsilon_kappa*exp(1i*phase_kappa_val(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    [~, ~, ~, ~, e_loss, ~] = get_Energy_regimes(k_tr, omega_complex, Omega, rs, ks, vr, gamma, len, N_fixed);

    if isnan(e_loss)
        objective_val = large_penalty_value; return;
    end
    if isinf(e_loss)
        objective_val = e_loss; return; % If e_loss=+Inf, obj=+Inf (bad). If e_loss=-Inf, obj=-Inf (good).
    end
    objective_val = e_loss;
end

% NOTE: The function get_Energy_regimes(k_tr, omega, Omega, rs, ks, vr, gamma, len, N)
% is assumed to be defined elsewhere and accessible on MATLAB's path.
% Its expected return signature is:
% [e_gain, as_gain_val, e_max_diss_val, e_min_diss_val, e_loss, as_loss_val]
