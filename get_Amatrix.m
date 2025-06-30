function [A] = get_Amatrix(gamma,mu,xi,v,k_tr,cjs,F,G,len)
%GET_Gmatrix Computes the matrix \mathcal{A} which has the coefficients g_{mn}
%   gamma:  constant st delta=gamma*len^2
%   mu:     constant st omega=mu*len
%   xi:     constant st Omega=xi*len
%   v:      wave speed inside the resonators
%   k_tr:   truncation parameter K
%   cjs:    constants c_j st lambda_j=c_j*len
%   F:      matrix of all eigenvectors
%   G:      inverse matrix of F
%   len:    length of resonators

ns = -k_tr:k_tr; % vector of indices
Xinv = diag(v./(mu+ns.*xi));
Z = diag(sqrt(-1).*(cjs.^2));
% compute the matrix \mathcal{A}
A = len/(2*gamma).*Xinv*F*Z*G;

end