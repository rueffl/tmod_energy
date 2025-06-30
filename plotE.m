%% Settings for the structure's geometry

N = 3; % number of resonators
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

%% Iterate over mu and epsilon_kappa

% Settings for the material parameters
gamma = 0.05; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
all_mu = linspace(0,0.99,101); % scaling coefficients for \omega

% incident wave field
a = ones(2*k_tr+1,1); a = a./norm(a,2);

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
all_epsk = linspace(0,0.99,101); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho. It needs to be 0, don't change!
phase_kappa = ones(1,N).*pi/2; % modulation phases of kappa, we assume that \kappa_i(t) is the same accross all resonators
phase_rho = zeros(1,N); % modulation phases of rho

i_k = 1;
E = zeros(length(all_epsk),length(all_mu));

for epsilon_kappa = all_epsk
    i_w = 1;
    for mu = all_mu

        omega = mu*len;
        O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
        rs = []; % Fourier coefficients of 1/rho
        ks = []; % Fourier coefficients of 1/kappa
        for j = 1:N
            rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
            ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
            ks = [ks; ks_j];
            rs = [rs; rs_j];
        end

        A_sum = zeros(2*k_tr+1,2*k_tr+1);
        for i = 1:N
            % Compute the eigenvalues lambda and the eigenvector f of the matrix C_i
            Ci = get_Ci(k_tr,i,omega,Omega,rs,ks,vr); % matrix C_i fo the interior ODE
            [fs,lambdas_square] = eig(Ci); % take the eigenvalues and eigenvectors of C_i
            % fs = fliplr(flip(fs));
            lambdas = sqrt(diag(lambdas_square)); % lambda_j^1 in the exponents of the interior solution
            cjs = lambdas./(len); % constants c_j st \lambda_j=c_j\ell, note that it is the same for each resonator and can therefore be computed for just D_1
            gs = inv(fs); % matrix G which is the inverse of the matrix F containing the eigenvectors f_n
            
            % Construct the matrix g of scattering coefficients
            A = get_Amatrix(gamma,mu,xi,vr,k_tr,cjs,fs,gs,len);
            A_sum = A_sum+A;
        end
        g_sum = inv((eye(2*k_tr+1)-A_sum))*A_sum;

        % Compute the total energy flux for a given incident wave field
        E(i_k,i_w) = get_E(g_sum,a,O);
        i_w = i_w+1;

    end
    i_k = i_k+1;

end

% Create Plot
figure()
hold on
surf(all_epsk,all_mu,E)
xlabel('$\varepsilon_{\kappa}$',Interpreter='latex')
ylabel('$\mu$',Interpreter='latex')
colorbar
shading interp

%% Iterate over a

% Settings for the material parameters
gamma = 0.05; delta = gamma*len^2; % small contrast parameter
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators
mu = 0.9*sqrt(-1); omega = mu*len; % scaling coefficients for \omega
kr = omega/vr; % wave number inside the resonator
k = omega/v0; % wave number outside of the resonator

% Settings for modulation
xi = 0.2; Omega = xi*len; % modulation frequency
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
epsilon_kappa = 0.9; % modulation amplitude of kappa
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

A_sum = zeros(2*k_tr+1,2*k_tr+1);
for i = 1:N
    % Compute the eigenvalues lambda and the eigenvector f of the matrix C_i
    Ci = get_Ci(k_tr,i,omega,Omega,rs,ks,vr); % matrix C_i fo the interior ODE
    [fs,lambdas_square] = eig(Ci); % take the eigenvalues and eigenvectors of C_i
    % fs = fliplr(flip(fs));
    lambdas = sqrt(diag(lambdas_square)); % lambda_j^1 in the exponents of the interior solution
    cjs = lambdas./(len); % constants c_j st \lambda_j=c_j\ell, note that it is the same for each resonator and can therefore be computed for just D_1
    gs = inv(fs); % matrix G which is the inverse of the matrix F containing the eigenvectors f_n
    
    % Construct the matrix g of scattering coefficients
    A = get_Amatrix(gamma,mu,xi,vr,k_tr,cjs,fs,gs,len);
    A_sum = A_sum+A;
end
g_sum = inv((eye(2*k_tr+1)-A_sum))*A_sum;


as0 = linspace(-4,4,1000);
as1 = linspace(-4,4,1000);
i_a0 = 1;
E = zeros(length(as1),length(as0));

for a0 = as1 
    i_a1 = 1;
    for a1 = as0
        % incident wave field
        a = zeros(2*k_tr+1,1); 
        a(k_tr+1) = a0; a(k_tr) = a1; a(k_tr+2) = a1; a = a./norm(O*a,2);
        % Compute the total energy flux for a given incident wave field
        E(i_a1,i_a0) = get_E(g_sum,a,O);
        i_a1 = i_a1+1;
    end
    i_a0 = i_a0+1;
end

% Create Plot
figure()
h_surf = surf(as1,as0,E,'EdgeColor','none');
hold on
set(h_surf, 'HandleVisibility', 'off');
shading interp

% Identify indices where E is within epsilon of 1
epsilon = 0.001;
idx = abs(E - 1) < epsilon;
% Extract the corresponding coordinates
[X, Y] = meshgrid(as1, as0);
x_highlight = X(idx);
y_highlight = Y(idx);
z_highlight = E(idx);
% Overlay red dots at the identified points
h_close_to_one = plot3(x_highlight, y_highlight, z_highlight, 'r.', 'MarkerSize', 10, 'LineWidth', 2);

% Find the global minimum
[minVal, minIdx] = min(E(:));
[rowMin, colMin] = ind2sub(size(E), minIdx);
h_min = plot3(as1(colMin), as0(rowMin), minVal, '*g', 'MarkerSize', 10, 'LineWidth', 2);
% Find the global maximum
[maxVal, maxIdx] = max(E(:));
[rowMax, colMax] = ind2sub(size(E), maxIdx);
h_max = plot3(as1(colMax), as0(rowMax), maxVal, '*b', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('$a_1$',Interpreter='latex')
ylabel('$a_0$',Interpreter='latex')
cb = colorbar;
cb.Label.String = 'E';
legend([h_close_to_one, h_min, h_max], {'E ≈ 1', 'Minimum', 'Maximum'});

view(0,90)
