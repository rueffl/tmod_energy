function [G] = get_Gmat(k_tr,omega,Omega,rs,ks,vr,gamma,len,N)
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

mu = omega/len;
xi = Omega/len;
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
G = inv((eye(2*k_tr+1)-A_sum))*A_sum;

end