function [E] = get_E_nonasympt(S11,S21,a,O)
%GET_E_NONASYMPT Computes the total energy flux E for large l
%   S11:    upper left block of the scattering matrix
%   S21:    lower left block of the scattering matrix
%   a:      vector of coefficients of the incident wave field
%   O:      diagonal matrix with frequencies

E = norm(O*S11*a,2)^2+norm(O*S21*a,2)^2;

end