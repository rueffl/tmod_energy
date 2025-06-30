% Geometric setup
N = 2; % number of resonators
len = 0.1; li = ones(1,N).*len; % length of the resonator
spacing = 1*len; lij = ones(1,N-1).*spacing; % spacing between the resonators
xm = [0]; % left boundary points of the resonators
for i = 2:N
    xm = [xm,xm(end)+len+spacing];
end
xp = xm + li; % right boundary points of the resonators
zi = (xm+xp)./2; % centre points of the resonators
k_tr = 3; % truncation parameter

% Settings for the material parameters
delta = 1*10^(-3); % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Settings for modulation
Omega = 1; % modulation frequency
eks = linspace(0,0.99,50); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = ones(1,N).*pi/2; % modulation phases of rho

C = make_capacitance_finite(N,lij); % capacitance matrix
ik = 1;
omegas = zeros(length(eks),2*N);
for epsilon_kappa = eks

    % Compute the resonant frequencies 
    omegas_eks = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0);
%     omegas_eks = get_capacitance_approx_spec_hot(epsilon_kappa, phase_kappa, epsilon_rho, phase_rho, Omega, delta, vr, li, C, v0);
    [sorted_o, index] = sort(imag(omegas_eks)); omegas(ik,:) = omegas_eks(index);

    ik = ik+1;

end

% Create Plot
expt_idx = 17; % index of eks where an exceptional point occurs
cmap = parula(2*N+2);
figure()
hold on
subplot(1,2,1)
hold on
for i = 1:2*N
    if i == 1
        plot(eks,-real(omegas(:,i)),'.','Color','k','MarkerSize',8)
    end
    plot(eks,real(omegas(:,i)),'.','Color','k','MarkerSize',8)
end
plot([eks(expt_idx),eks(expt_idx)],[min(min(real(omegas)))-0.05,max(max(real(omegas)))+0.05],'--','LineWidth',1)
grid on
ylim([min(min(real(omegas)))-0.05,max(max(real(omegas)))+0.05])
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('Re$\left(\omega_0\right)$','Interpreter','latex')
subplot(1,2,2)
hold on
for i = 1:2*N
    plot(eks,imag(omegas(:,i)),'.','Color','k','MarkerSize',8)
end
plot([eks(expt_idx),eks(expt_idx)],[min(min(imag(omegas)))-0.05,max(max(imag(omegas)))+0.05],'--','LineWidth',1)
grid on
ylim([min(min(imag(omegas)))-0.05,max(max(imag(omegas)))+0.05])
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('Im$\left(\omega_0\right)$','Interpreter','latex')


%% Calculate E_{gain} and E_{loss} at the herein computed frequencies

ik = 1;
es_gain = zeros(length(eks),2*N); es_loss = zeros(length(eks),2*N);
for epsilon_kappa = eks
    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    for iN = 1:2*N
        if real(omegas(ik,iN))>0
            omega = real(omegas(ik,iN));
            % Compute the different regimes
            [es_gain(ik,iN), as_gain, es_cons, as_cons, es_loss(ik,iN), as_loss] = get_Energy_regimes_nonasympt(k_tr,omega,Omega,rs,ks,vr,delta,N,xm,xp);
        end
    end
    ik = ik+1;
end

% Create Plot
cmap = parula(2*N+2);
figure()
hold on
subplot(1,2,1)
hold on
for i = 1:2*N
    plot(eks(2:end),es_gain(2:end,i),'.','Color','k','MarkerSize',8)
end
plot([eks(expt_idx),eks(expt_idx)],[min(min(es_gain))-0.05,max(max(es_gain))+1],'--','LineWidth',1)
grid on
ylim([min(min(es_gain))-0.05,max(max(es_gain))+1])
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$E_{\mathrm{gain}}$','Interpreter','latex')
subplot(1,2,2)
hold on
for i = 1:2*N
    plot(eks(2:end),es_loss(2:end,i),'.','Color','k','MarkerSize',8)
end
plot([eks(expt_idx),eks(expt_idx)],[min(min(es_loss))-0.05,max(max(es_loss))+1],'--','LineWidth',1)
grid on
ylim([min(min(es_loss))-0.05,max(max(es_loss))+1])
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$E_{\mathrm{loss}}$','Interpreter','latex')

