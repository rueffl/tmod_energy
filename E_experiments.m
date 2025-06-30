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
eks = linspace(0,0.95,10); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = zeros(1,N); % modulation phases of rho

es_loss = zeros(1,length(eks)); es_gain = zeros(1,length(eks));
ik = 1;
for epsilon_kappa = eks

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    % Compute the different regimes
    [es_gain(ik), as_gain, e_cons, as_cons, es_loss(ik), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    ik = ik+1;

end

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
phase_rho = ones(1,N).*pi/2; % modulation phases of rho

es_loss = zeros(1,length(ers)); es_gain = zeros(1,length(ers));
ir = 1;
for epsilon_rho = ers

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    % Compute the different regimes
    [es_gain(ir), as_gain, e_cons, as_cons, es_loss(ir), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    ir = ir+1;

end

% Create Plot
figure()
hold on
plot(ers,es_gain,'-*g')
plot(ers,es_loss,'-*r')
grid on
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$\varepsilon_{\rho}$','Interpreter','latex')


%% Iterate over epsilon_rho and epsilon_kappa

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
eks = linspace(0,0.99,100); % modulation amplitude of kappa
ers = linspace(0,0.99,100); % modulation amplitude of rho
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = ones(1,N).*pi/2; % modulation phases of rho

es_loss = zeros(length(eks),length(ers)); es_gain = zeros(length(eks),length(ers));
ir = 1;
for epsilon_rho = ers
    ik = 1;
    for epsilon_kappa = eks

        rs = []; % Fourier coefficients of 1/rho
        ks = []; % Fourier coefficients of 1/kappa
        for j = 1:N
            rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
            ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
            ks = [ks; ks_j];
            rs = [rs; rs_j];
        end
        % Compute the different regimes
        [es_gain(ik,ir), as_gain, e_cons, as_cons, es_loss(ik,ir), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
        ik = ik+1;

    end
    ir = ir+1;

end

% Create Plot
figure()
hold on
subplot(1,2,1)
surf(eks,ers,es_gain)
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$\varepsilon_{\rho}$','Interpreter','latex')
colorbar
shading interp
xlim([0,0.99])
ylim([0,0.99])
view(0,90)
subplot(1,2,2)
hold on
surf(eks,ers,es_loss)
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$\varepsilon_{\rho}$','Interpreter','latex')
colorbar
shading interp
xlim([0,0.99])
ylim([0,0.99])
view(0,90)

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
for epsilon_rho = ers

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    % Compute the different regimes
    [es_gain(ir), as_gain, e_cons, as_cons, es_loss(ir), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    ir = ir+1;

end

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
for N = Ns

    % Settings for the material parameters
    gamma = 1; delta = gamma*len^2; % small contrast parameter
    vr = 1; % wave speed inside the resonators
    v0 = 1; % wave speed outside the resonators

    % Settings for modulation
    epsilon_kappa = 0.5; % modulation amplitude of kappa
    epsilon_rho = 0.5; % modulation amplitude of rho. It needs to be 0, don't change!
    phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
    phase_rho = ones(1,N).*pi/2; % modulation phases of rho
    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    len = 0.01; li = ones(1,N).*len; % length of the resonator
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
    omega = 0.00024;
    
    % Settings for the material parameters
    mu = omega/len; % operating frequency
    kr = omega/vr; % wave number inside the resonator
    k = omega/v0; % wave number outside of the resonator
    xi = 0.31; Omega = xi*len; % modulation frequency
    O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);

    % Compute the different regimes
    [es_gain(iN), as_gain, e_cons, as_cons, es_loss(iN), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    iN = iN+1;
    
    end

% Create Plot
figure()
hold on
plot(Ns,es_gain,'-*g')
plot(Ns,es_loss,'-*r')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$N$','Interpreter','latex')
grid on


%% Iterate over omega

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
gamma = 0.5; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mus = linspace(0.001,1,20); % scaling coefficient of omega

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
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
for mu = mus
    
    omega = mu*len; % operating frequency
    kr = omega/vr; % wave number inside the resonator
    k = omega/v0; % wave number outside of the resonator
    % Compute the different regimes
    [es_gain(im), as_gain, e_cons, as_cons, es_loss(im), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    im = im+1;
    
end

% Compute the resonant frequency as a reference
C = make_capacitance_finite(N,lij);
w_out = real(get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr));

% Create Plot
figure()
hold on
plot(omegas,es_gain,'-*g')
plot(omegas,es_loss,'-*r')
plot(max(w_out).*ones(2,1),[min(es_loss)-0.01,max(es_gain)+0.01],'--k')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','$\omega_0$','Interpreter','latex')
xlabel('$\omega$','Interpreter','latex')
grid on


%% Iterate over epsilon_rho and omega

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
gamma = 0.5; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mus = linspace(0.001,1,50); omega = mu*len; % operating frequency
kr = omega/vr; % wave number inside the resonator
k = omega/v0; % wave number outside of the resonator

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
eks = linspace(0,0.99,50); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = ones(1,N).*pi/2; % modulation phases of rho

es_loss = zeros(length(eks),length(mus)); es_gain = zeros(length(eks),length(mus));
im = 1;
for mu = mus
    ik = 1;
    for epsilon_kappa = eks

        rs = []; % Fourier coefficients of 1/rho
        ks = []; % Fourier coefficients of 1/kappa
        for j = 1:N
            rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
            ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
            ks = [ks; ks_j];
            rs = [rs; rs_j];
        end
        % Compute the different regimes
        [es_gain(ik,im), as_gain, e_cons, as_cons, es_loss(ik,im), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
        ik = ik+1;

    end
    im = im+1;

end

% Create Plot
figure()
hold on
subplot(1,2,1)
surf(eks,mus,es_gain)
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$\mu$','Interpreter','latex')
colorbar
shading interp
xlim([0,0.99])
ylim([0,0.99])
view(0,90)
subplot(1,2,2)
hold on
surf(eks,mus,es_loss)
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')
ylabel('$\mu$','Interpreter','latex')
colorbar
shading interp
xlim([0,0.99])
ylim([0,0.99])
view(0,90)


%% Iterate over epsilon_kappa and choose omega_0

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
gamma = 0.5; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
eks = linspace(0,0.95,10); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = ones(1,N).*pi/2; % modulation phases of rho

es_loss = zeros(1,length(eks)); es_gain = zeros(1,length(eks)); omegas = zeros(1,length(eks));
ik = 1;
for epsilon_kappa = eks

    % Compute omega0
    if N > 1
        C = make_capacitance_finite(N,lij);
        w_out = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr);
        w_out = w_out(w_out >= 0);
    else
        w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0);
    end
    omega = max(w_out); omegas(ik) = omega;   

    kr = omega/vr; % wave number inside the resonator
    k = omega/v0; % wave number outside of the resonator
    O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    % Compute the different regimes
    [es_gain(ik), as_gain, e_cons, as_cons, es_loss(ik), as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);
    ik = ik+1;

end

% Create Plot
figure()
hold on
plot(eks,es_gain,'-*g')
plot(eks,es_loss,'-*r')
legend('$E_{\mathrm{gain}}$','$E_{\mathrm{loss}}$','Interpreter','latex')
xlabel('$\varepsilon_{\kappa}$','Interpreter','latex')


