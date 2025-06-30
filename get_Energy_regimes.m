function [e_gain, as_gain, e_cons, as_cons, e_loss, as_loss] = get_Energy_regimes(k_tr,omega,Omega,rs,ks,vr,gamma,len,N)
%GET_ENERGY_REGIMES Determines the regime of the incident wave where energy is conserved, gained and lost
%   k_tr:       Truncation parameter
%   omega:      Operating frequency
%   Omega:      Modulation frequency
%   rs:         Fourier coefficients of 1/\rho
%   ks:         Fourier coeddicients of 1/\kappa
%   vr:         Interior wave speed
%   gamma:      \delta = \gamma\ell^2
%   len:        Length of each resonator
%   N:          Number of resonators

O = diag(omega.*ones(1,2*k_tr+1)+[-k_tr:k_tr].*Omega);
g_sum = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N);


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
        E(i_a1,i_a0) = get_E(g_sum,a,O);
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