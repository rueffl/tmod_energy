function [Escat] = get_Escat(g,a,O)
%GET_E Computes the total energy flux E
%   g:  Matrix with scattering coefficients
%   a:  vector of coefficients of the incident wave field
%   O:  diagonal matrix with frequencies

size_ktr = length(a);
Escat = 2*norm(O*g*a,2)^2;

end