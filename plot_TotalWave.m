%% Settings for the structure's geometry

N = 10; % number of resonators
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

% Define evaluation points
len_xs = 200;
len_zs = 20;
xs = linspace(xm(1)-spacing,xp(end)+spacing,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end
ts = linspace(0,100,11); % time steps

% Compute the different regimes
[e_gain, as_gain, e_cons, as_cons, e_loss, as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);

% Define incident field
O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
a0 = as_gain(1,1); a1 = as_gain(1,2);
as = zeros(2*k_tr+1,1); as(k_tr+1) = a0; as(k_tr) = a1; as(k_tr+2) = a1; as = as./norm(O*as,2); % vector of coefficients of the incident field
uin = @(x,t) as'*exp(sqrt(-1)*((k).*x+(omega+[-k_tr:k_tr]'.*Omega).*t)); % incident wave
vin = @(x,n) as(n+k_tr+1)*exp(sqrt(-1)*k*x);
dx_vin = @(x,n) as(n+k_tr+1)*sqrt(-1)*k*exp(sqrt(-1)*k*x);

% prepare figure
fig = figure();
c_map = parula(length(ts)+2); ic = 1;

for t = ts

    % Compute solution coefficients
    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,omega,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin, @(x,n) 0); % vector \mathcal{F}
    N_vin = getN_vin(k_tr, N, delta, xm, xp, omega, Omega, v0, lij, vin, @(x,n) 0); % matrix-vector product \mathcal{N}v^{in}
    RHS = MatcalF + N_vin;
    sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}
    
    us_eval_x = zeros(1,len_xs);
    us_eval_z = zeros(1,len_zs);
    usx = @(x) uin(x,t) + get_us(x, t, N, xm, xp, lij, k_tr, v0, omega, Omega, rs, ks, vr, sol, k, vin); % scattered wave field as a function of x for fixed time t, according to formula (31)

    for i = 1:len_xs
        us_eval_x(i) = usx(xs(i));
    end
    
    % plot results
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(xs,real(us_eval_x),'-','DisplayName',strcat('$t=$ ',num2str(t)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(xs,imag(us_eval_x),'-','DisplayName',strcat('$t=$ ',num2str(t)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    for i = 1:N
        for j = 1:len_zs
            us_eval_z(j) = usx(zs(i,j));
        end
        subplot(1,2,1)
        plot(zs(i,:),real(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
        subplot(1,2,2)
        plot(zs(i,:),imag(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
    end

    ic = ic+1;
    
end

subplot(1,2,1)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Re$(u(x,t))$',Interpreter='latex',FontSize=18)

subplot(1,2,2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Im$(u(x,t))$',Interpreter='latex',FontSize=18)

legend('show',interpreter='latex',fontsize=18,location='southoutside')


