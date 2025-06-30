function [e_gain, as_gain, e_cons, as_cons, e_loss, as_loss] = get_Energy_regimes_nonasympt(k_tr,omega,Omega,rs,ks,vr,delta,N,xm,xp)
%GET_ENERGY_REGIMES_NONASYMPT Determines the regime of the incident wave where energy is conserved, gained and lost for large len
%   k_tr:       Truncation parameter
%   omega:      Operating frequency
%   Omega:      Modulation frequency
%   rs:         Fourier coefficients of 1/\rho
%   ks:         Fourier coeddicients of 1/\kappa
%   vr:         Interior wave speed
%   delta:      Contrast parameter
%   N:          Number of resonators
%   xm:         Left boundary points of resonators
%   xp:         Right boundary points of resonators

O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
k = (omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega)./vr; % interior wave number for each considered mode
% Compute the transfer matrix
Stilde_tot = eye(2*(2*k_tr+1));
for i = 1:N
    Ci = get_Ci(k_tr,i,omega,Omega,rs,ks,vr);
    [F,lambdas_square] = eig(Ci); % take the eigenvalues and eigenvectors of C
    lambdas = diag(sqrt(lambdas_square));
    Si = get_tdepStilde(xm(i), xp(i), delta, k, F, lambdas);
%     Si = get_tdepStilde(xm(i), xp(i), delta, k_tr, vr, omega, Omega, F, lambdas);
    Stilde_tot = Si*Stilde_tot;
end

% Construct the left blocks of the scattering matrix
S11 = Stilde_tot(1:(2*k_tr+1), 1:(2*k_tr+1))-Stilde_tot(1:(2*k_tr+1), (2*k_tr+2):end)*inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));
S21 = -inv(Stilde_tot((2*k_tr+2):end, (2*k_tr+2):end))*Stilde_tot((2*k_tr+2):end, 1:(2*k_tr+1));

as0 = linspace(-2,2,100);
as1 = linspace(-2,2,100);
i_a0 = 1;
E = zeros(length(as1),length(as0));

for a0 = as1 
    i_a1 = 1;
    for a1 = as0
        % incident wave field
        a = zeros(2*k_tr+1,1); 
        a(k_tr+1) = a0; a(k_tr) = a1; a(k_tr+2) = a1; a = a./norm(O*a,2);
        % Compute the total energy flux for a given incident wave field
        E(i_a1,i_a0) = get_E_nonasympt(S11,S21,a,O);
        i_a1 = i_a1+1;
    end
    i_a0 = i_a0+1;
end

% Identify indices where E is within epsilon of 1
epsilon = 0.0001;
idx = abs(E - 1) < epsilon;
[X, Y] = meshgrid(as1, as0);
as_cons =[Y(idx),X(idx)];
e_cons = E(idx);
% Find the global minimum
[e_loss, minIdx] = min(E(:));
[rowMin, colMin] = ind2sub(size(E), minIdx);
as_loss = [as0(rowMin), as1(colMin)];
% Find the global maximum
[e_gain, maxIdx] = max(E(:));
[rowMax, colMax] = ind2sub(size(E), maxIdx);
as_gain = [as0(rowMax), as1(colMax)];



end