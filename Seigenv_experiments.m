%% Iterate over epsilon_kappa

N = 1; % number of resonators
len = 0.1; li = ones(1,N).*len; % length of the resonator
U = 4000; % length of the domain
spacing = 2*len; lij = ones(1,N-1).*spacing; % spacing between the resonators
xm = [0]; % left boundary points of the resonators
for i = 2:N
    xm = [xm,xm(end)+len+spacing];
end
xp = xm + li; % right boundary points of the resonators
zi = (xm+xp)./2; % centre points of the resonators
L = spacing + len; % length of unit cell
k_tr = 4; % truncation parameter

% Settings for the material parameters
gamma = 0.05; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mu = 0.9; omega = mu*len; % operating frequency
kr = omega/vr; % wave number inside the resonator
k = omega/v0; % wave number outside of the resonator

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
eks = linspace(0,0.99,10); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = zeros(1,N); % modulation phases of rho

es_loss = zeros(1,length(eks)); es_gain = zeros(1,length(eks));
ik = 1;
Es = zeros(length(eks),2*k_tr+1);
cmap = parula(length(eks));
figure()
hold on
for epsilon_kappa = eks

    % fix the time dependent material parameters
    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    % Compute the total energy with the eigenvectors of S as the incident field
    Stilde_tot = eye(2*(2*k_tr+1));
    for i = 1:N
        C = get_Ci(k_tr,i,omega,Omega,rs,ks,vr);
        [F,lambdas_square] = eig(C); % take the eigenvalues and eigenvectors of C
        lambdas = diag(sqrt(lambdas_square));
        Si = get_tdepStilde(xm(i), xp(i), delta, k, F, lambdas);
        Stilde_tot = Si*Stilde_tot;
    end
    % Construct the scattering matrix
    S11 = Stilde_tot(1:(2*k_tr+1), 1:(2*k_tr+1))-Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S12 = Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot(2,2));
    S21 = -inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S22 = inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end));
    S_tot = [S11, S12; S21, S22];

    [V,D] = eig(S_tot); eigenvals = diag(D); eigenvals = eigenvals(1:(2*k_tr+1));
    [sorted_z, index] = sort(abs(eigenvals));
    eigenvals = eigenvals(index); V = V(:,index);
    G = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    for iv = 1:2*k_tr+1
        Vi_normalised = V(1:(2*k_tr+1),iv)./(norm(O*V(1:(2*k_tr+1),iv)));
        Es(ik,iv) = get_E(G,Vi_normalised,O);
    end
    es_gain(ik) = max(Es(ik,:)); es_loss(ik) = min(Es(ik,:));
    plot(abs(eigenvals), Es(ik,:), 'x-', 'Color', cmap(ik, :), 'LineWidth', 2, 'DisplayName', strcat('$\varepsilon_{\kappa}=$ ', num2str(epsilon_kappa)))
    ik = ik+1;

end

legend('Interpreter','latex')
ylabel('$E$','Interpreter','latex')

% Create Plot
figure()
hold on
plot(eks,es_gain,'-*g')
plot(eks,es_loss,'-*r')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')


%% Iterate over epsilon_rho

N = 1; % number of resonators
len = 0.1; li = ones(1,N).*len; % length of the resonator
U = 4000; % length of the domain
spacing = 2*len; lij = ones(1,N-1).*spacing; % spacing between the resonators
xm = [0]; % left boundary points of the resonators
for i = 2:N
    xm = [xm,xm(end)+len+spacing];
end
xp = xm + li; % right boundary points of the resonators
zi = (xm+xp)./2; % centre points of the resonators
L = spacing + len; % length of unit cell
k_tr = 4; % truncation parameter

% Settings for the material parameters
gamma = 0.05; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mu = 0.9; omega = mu*len; % operating frequency
kr = omega/vr; % wave number inside the resonator
k = omega/v0; % wave number outside of the resonator

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
epsilon_kappa = 0; % modulation amplitude of kappa
ers = linspace(0,0.99,10); % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = zeros(1,N); % modulation phases of rho

es_loss = zeros(1,length(ers)); es_gain = zeros(1,length(ers));
ir = 1;
Es = zeros(length(ers),2*k_tr+1);
cmap = parula(length(ers));
figure()
hold on
for epsilon_rho = ers

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    
    % Compute the total energy with the eigenvectors of S as the incident field
    Stilde_tot = eye(2*(2*k_tr+1));
    for i = 1:N
        C = get_Ci(k_tr,i,omega,Omega,rs,ks,vr);
        [F,lambdas_square] = eig(C); % take the eigenvalues and eigenvectors of C
        lambdas = diag(sqrt(lambdas_square));
        Si = get_tdepStilde(xm(i), xp(i), delta, k, F, lambdas);
        Stilde_tot = Si*Stilde_tot;
    end
    % Construct the scattering matrix
    S11 = Stilde_tot(1:(2*k_tr+1), 1:(2*k_tr+1))-Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S12 = Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot(2,2));
    S21 = -inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S22 = inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end));
    S_tot = [S11, S12; S21, S22];

    [V,D] = eig(S_tot); eigenvals = diag(D); eigenvals = eigenvals(1:(2*k_tr+1));
    [sorted_z, index] = sort(abs(eigenvals));
    eigenvals = eigenvals(index); V = V(:,index);
    G = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    for iv = 1:2*k_tr+1
        Vi_normalised = V(1:(2*k_tr+1),iv)./(norm(O*V(1:(2*k_tr+1),iv)));
        Es(ik,iv) = get_E(G,Vi_normalised,O);
    end
    es_gain(ir) = max(Es(ir,:)); es_loss(ir) = min(Es(ir,:));
    plot(abs(eigenvals), Es(ir,:), 'x-', 'Color', cmap(ir, :), 'LineWidth', 2, 'DisplayName', strcat('$\varepsilon_{\rho}=$ ', num2str(epsilon_rho)))
    ir = ir+1;

end

legend('Interpreter','latex')
ylabel('$E$','Interpreter','latex')

% Create Plot
figure()
hold on
plot(ers,es_loss,'-*r')
plot(ers,es_gain,'-*g')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$\varepsilon_{\rho}$','Interpreter','latex')

%% Iterate over N

Ns = 1:2:150; % number of resonators
len = 0.1;

es_loss = zeros(1,length(Ns)); es_gain = zeros(1,length(Ns));
iN = 1;
Es = zeros(length(Ns),2*k_tr+1);
cmap = parula(length(Ns));
figure()
hold on
for N = Ns

    % Settings for the material parameters
    gamma = 0.05; delta = gamma*len^2; % small contrast parameter
    vr = 1; % wave speed inside the resonators
    v0 = 1; % wave speed outside the resonators

    % Settings for modulation
    epsilon_kappa = 0.8; % modulation amplitude of kappa
    epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
    phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
    phase_rho = zeros(1,N); % modulation phases of rho
    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    len = 0.1; li = ones(1,N).*len; % length of the resonator
    U = 4000; % length of the domain
    spacing = 2*len; lij = ones(1,N-1).*spacing; % spacing between the resonators
    xm = [0]; % left boundary points of the resonators
    for i = 2:N
        xm = [xm,xm(end)+len+spacing];
    end
    xp = xm + li; % right boundary points of the resonators
    zi = (xm+xp)./2; % centre points of the resonators
    L = spacing + len; % length of unit cell
    k_tr = 4; % truncation parameter

    % Compute the resonant frequency as a reference
%     if N > 1
%         C = make_capacitance_finite(N,lij);
%         w_out = real(get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr));
%         w_out = w_out(w_out >= 0);
%     else
%         w_out = real(get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0));
%     end
%     omega = max(w_out);
    omega = 0.0002;
    
    % Settings for the material parameters
    mu = omega/len; % operating frequency
    kr = omega/vr; % wave number inside the resonator
    k = omega/v0; % wave number outside of the resonator
    xi = 0.2; Omega = xi*len; % modulation frequency
    O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);

    % Compute the total energy with the eigenvectors of S as the incident field
    Stilde_tot = eye(2*(2*k_tr+1));
    for i = 1:N
        C = get_Ci(k_tr,i,omega,Omega,rs,ks,vr);
        [F,lambdas_square] = eig(C); % take the eigenvalues and eigenvectors of C
        lambdas = diag(sqrt(lambdas_square));
        Si = get_tdepStilde(xm(i), xp(i), delta, k, F, lambdas);
        Stilde_tot = Si*Stilde_tot;
    end
    % Construct the scattering matrix
    S11 = Stilde_tot(1:(2*k_tr+1), 1:(2*k_tr+1))-Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S12 = Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot(2,2));
    S21 = -inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S22 = inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end));
    S_tot = [S11, S12; S21, S22];

    [V,D] = eig(S_tot); eigenvals = diag(D); eigenvals = eigenvals(1:(2*k_tr+1));
    [sorted_z, index] = sort(abs(eigenvals));
    eigenvals = eigenvals(index); V = V(:,index);
    G = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    for iv = 1:2*k_tr+1
        Vi_normalised = V(1:(2*k_tr+1),iv)./(norm(O*V(1:(2*k_tr+1),iv)));
        Es(ik,iv) = get_E(G,Vi_normalised,O);
    end
    es_gain(iN) = max(Es(iN,:)); es_loss(iN) = min(Es(iN,:));
    plot(abs(eigenvals), Es(iN,:), 'x-', 'Color', cmap(iN, :), 'LineWidth', 2, 'DisplayName', strcat('$N=$ ', num2str(N)))
    iN = iN+1;
    
end

legend('Interpreter','latex')
ylabel('$E$','Interpreter','latex')

% Create Plot
figure()
hold on
plot(Ns,es_gain,'-*g')
plot(Ns,es_loss,'-*r')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$N$','Interpreter','latex')


%% Iterate over omega

N = 2; % number of resonators
len = 0.1; li = ones(1,N).*len; % length of the resonator
U = 4000; % length of the domain
spacing = 2*len; lij = ones(1,N-1).*spacing; % spacing between the resonators
xm = [0]; % left boundary points of the resonators
for i = 2:N
    xm = [xm,xm(end)+len+spacing];
end
xp = xm + li; % right boundary points of the resonators
zi = (xm+xp)./2; % centre points of the resonators
L = spacing + len; % length of unit cell
k_tr = 4; % truncation parameter

% Settings for the material parameters
gamma = 0.05; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mus = linspace(0.001,1,20); % scaling coefficient of omega

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
epsilon_kappa = 0.5; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = zeros(1,N); % modulation phases of rho
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

es_loss = zeros(1,length(mus)); es_gain = zeros(1,length(mus));
omegas = mus.*len;
im = 1;
Es = zeros(length(mus),2*k_tr+1);
cmap = parula(length(mus));
figure()
hold on
for mu = mus
    
    omega = mu*len; % operating frequency
    kr = omega/vr; % wave number inside the resonator
    k = omega/v0; % wave number outside of the resonator
    O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
    % Compute the total energy with the eigenvectors of S as the incident field
    Stilde_tot = eye(2*(2*k_tr+1));
    for i = 1:N
        C = get_Ci(k_tr,i,omega,Omega,rs,ks,vr);
        [F,lambdas_square] = eig(C); % take the eigenvalues and eigenvectors of C
        lambdas = diag(sqrt(lambdas_square));
        Si = get_tdepStilde(xm(i), xp(i), delta, k, F, lambdas);
        Stilde_tot = Si*Stilde_tot;
    end
    % Construct the scattering matrix
    S11 = Stilde_tot(1:(2*k_tr+1), 1:(2*k_tr+1))-Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S12 = Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot(2,2));
    S21 = -inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
    S22 = inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end));
    S_tot = [S11, S12; S21, S22];

    [V,D] = eig(S_tot); eigenvals = diag(D); eigenvals = eigenvals(1:(2*k_tr+1));
    [sorted_z, index] = sort(abs(eigenvals));
    eigenvals = eigenvals(index); V = V(:,index);
    G = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    for iv = 1:2*k_tr+1
        Vi_normalised = V(1:(2*k_tr+1),iv)./(norm(O*V(1:(2*k_tr+1),iv)));
        Es(ik,iv) = get_E(G,Vi_normalised,O);
    end
    es_gain(im) = max(Es(im,:)); es_loss(im) = min(Es(im,:));
    plot(abs(eigenvals), Es(im,:), 'x-', 'Color', cmap(im, :), 'LineWidth', 2, 'DisplayName', strcat('$\omega=$ ', num2str(omega)))
    im = im+1;
    
end

legend('Interpreter','latex')
ylabel('$E$','Interpreter','latex')

% Compute the resonant frequency as a reference
C = make_capacitance_finite(N,lij);
w_out = real(get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr));

% Create Plot
figure()
hold on
plot(omegas,es_gain,'-*g')
plot(omegas,es_loss,'-*r')
% plot(max(w_out).*ones(2,1),[min(es_loss)-0.01,max(es_gain)+0.01],'--k')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','$\omega_0$','Interpreter','latex')
xlabel('$\omega$','Interpreter','latex')
grid on