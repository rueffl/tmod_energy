function w_out = get_capacitance_approx_spec_hot(epsilon_kappa, phase_kappa, epsilon_rho, phase_rho, Omega, delta, vr, li, C, v0)
%GET_CAPACITANCE_APPROX_SPEC_V Solves the ODE with an additional 1/v * d/dt * psi term using a spectral method
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    modulation phase shift of kappa
%   epsilon_rho:    modulation amplitude of rho
%   phase_rho:      modulation phase shift of rho
%   Omega:          frequency of kappa_i and rho_i
%   delta:          contrast parameter
%   vr:             wave speed inside the resonators
%   li:             length of resonators
%   C:              capacitance matrix for fixed alpha
%   v0:             wave speed outside of the resonators


    alpha = delta.*diag(1./li)*diag(vr)^2;
    GCM = alpha*C;

    M = 1; % Maximal order of Fourier coefficients of 1/\kappa and 1/\rho
    N = size(GCM,1);
    N_fourier = 4; % Length of Fourier series approx
    
    R_mod = zeros(2*M+1,N); % Fourier coefficients of 1/\rho
    K_mod = zeros(2*M+1,N); % Fourier coefficients of 1/\kappa
    
    for i = 1:N        
        R_mod(:,i) = [epsilon_rho/2*exp(-1i*phase_rho(i)); 1; epsilon_rho/2*exp(1i*phase_rho(i))];
        K_mod(:,i) = [epsilon_kappa/2*exp(-1i*phase_kappa(i)); 1; epsilon_kappa/2*exp(1i*phase_kappa(i))]; 
    end

    ns = -N_fourier:N_fourier; % Fourier indicies
    NN = 2*N_fourier+1; % Total amount for fourier indicies
    O = diag(-1i.*ns.'*Omega); % Diagonal matrix consisting of n*Omega
    e = ones(NN,1); 
    INN = eye(NN);
    IN = eye(N);

    iK = zeros(NN*N); % Fourier coefficients of kappa
    iR = zeros(NN*N); % Fourier coefficients of rho
    R = zeros(NN*N); % Fourier coefficients of rho^-1
    
    for i = 1:N
        Ki = zeros(NN,NN);
        Ri = zeros(NN,NN);
        for m = -M:M
            Ki = Ki+diag(e(1:NN-abs(m))*K_mod(m+M+1,i),m);
            Ri = Ri+diag(e(1:NN-abs(m))*R_mod(m+M+1,i),m);
        end
        Ii = (i-1)*NN+1:i*NN;
        iK(Ii,Ii) = inv(Ki); % Fourier coefficients of \kappa
        R(Ii,Ii) = Ri; % Fourier coefficients of 1/\rho
        iR(Ii,Ii) = inv(Ri); % Fourier coefficients of \rho        
    end

    big_O = kron(IN,O);
    big_alpha = kron(alpha./v0,INN);
    big_GCM = kron(GCM,INN);

    upperleft = -1i*R*big_O*iR;
    upperright = 1i*R*iK;
    lowerleft = -1i*(iR*big_GCM-iR*big_alpha*R*(big_O*iR-iR*big_O));
    lowerright = -1i*(iR*big_alpha*R*iK+big_O);

    mat = [upperleft, upperright; lowerleft, lowerright];


    %mat = [-1i.*R*kron(IN,O)*iR+1i.*kron(IN,O), -1i*R*iK; 1i*iR*kron(GCM,INN)+1i.*kron(IN,O)*iR*kron(alpha,INN)./v0, 1i*iK*kron(alpha,INN)./v0-kron(IN,O)];
    
    
    w = sort(eig(mat),'ComparisonMethod','real');
    w_out = w(N*NN-N+1:N*NN+N);
end