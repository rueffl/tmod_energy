function [S] = get_tdepS(xm, xp, delta, k, F, lambdas)
%GET_TDEPS computes the scattering matrix corresponding to one resonator
%   xm:     left boundary point of the resonator
%   xp:     right boundary point of the resonator
%   delta:  contrast parameter
%   k:      vector of interior wave numbers for all considered modes
%   F:      matrix of eigenvectors
%   lambdas:eigenvalues of the resonator

    G = inv(F); 

    % Compute the diagonal matrices E_{-,r}, E_{-,l}, E_{+,r}, E_{+,l}
    E_plus_r = diag(exp(1i*k.*xp));  % E_{+,r}
    E_minus_l = diag(exp(-1i*k.*xm)); % E_{-,l}
    E_plus_l = diag(exp(1i*k.*xm));  % E_{+,l}
    E_minus_r = diag(exp(-1i*k.*xp)); % E_{-,r}
    
    % K matrix (assuming K is diagonal)
    K = diag(k);

    % Compute the matrices T^{\pm} and T^{\pm}_{\lambda}
    T_plus = T_blockdiag(lambdas, xp); T_plus_lambda = T_blockdiag2(lambdas, xp); 

    % Compute the matrices \hat{T} and \hat{T}_{\lambda}
    hat_T = That_blockdiag(lambdas, xm); hat_T_lambda = That_blockdiag2(lambdas, xm);

    % Compute the blocks of the transfer matrix
    Stilde11 = 0.5.*((E_minus_r*F*(T_plus*hat_T)+delta^(-1).*E_minus_r*inv(K)*F*(T_plus_lambda*hat_T))*G*E_plus_l+(delta.*E_minus_r*F*(T_plus*hat_T_lambda)+E_minus_r*inv(K)*F*(T_plus_lambda*hat_T_lambda))*G*K*E_plus_l);
    Stilde12 = 0.5.*((E_minus_r*F*(T_plus*hat_T)+delta^(-1).*E_minus_r*inv(K)*F*(T_plus_lambda*hat_T))*G*E_minus_l-(delta.*E_minus_r*F*(T_plus*hat_T_lambda)+E_minus_r*inv(K)*F*(T_plus_lambda*hat_T_lambda))*G*K*E_minus_l);
    Stilde21 = 0.5.*((E_plus_r*F*(T_plus*hat_T)-delta^(-1).*E_plus_r*inv(K)*F*(T_plus_lambda*hat_T))*G*E_plus_l+(delta.*E_plus_r*F*(T_plus*hat_T_lambda)-E_plus_r*inv(K)*F*(T_plus_lambda*hat_T_lambda))*G*K*E_plus_l);
    Stilde22 = 0.5.*((E_plus_r*F*(T_plus*hat_T)-delta^(-1).*E_plus_r*inv(K)*F*(T_plus_lambda*hat_T))*G*E_minus_l-(delta.*E_plus_r*F*(T_plus*hat_T_lambda)-E_plus_r*inv(K)*F*(T_plus_lambda*hat_T_lambda))*G*K*E_minus_l);

    % Construct the scattering matrix
    S11 = Stilde11-Stilde12*inv(Stilde22)*Stilde21;
    S12 = Stilde12*inv(Stilde22);
    S21 = -inv(Stilde22)*Stilde21;
    S22 = inv(Stilde22);

    S = [S11, S12; S21, S22];

end


function T = T_blockdiag(lambda, x_pm)
% Constructs T^{\pm} as a block-diagonal matrix
% Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]

    N = length(lambda);
    blocks = cell(1, N);  % Preallocate cell array for blocks
    
    for j = 1:N
        block = [exp(1i * lambda(j) * x_pm), exp(-1i * lambda(j) * x_pm)];
        blocks{j} = block;  % Store each 1x2 block
    end
    
    T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
end

function T = T_blockdiag2(lambda, x_pm)
% Constructs T^{\pm} as a block-diagonal matrix
% Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]

    N = length(lambda);
    blocks = cell(1, N);  % Preallocate cell array for blocks
    
    for j = 1:N
        block = [lambda(j)*exp(1i * lambda(j) * x_pm), -lambda(j)*exp(-1i * lambda(j) * x_pm)];
        blocks{j} = block;  % Store each 1x2 block
    end
    
    T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
end

function T = That_blockdiag(lambda, x_pm)
% Constructs T^{\pm} as a block-diagonal matrix
% Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]

    N = length(lambda);
    blocks = cell(1, N);  % Preallocate cell array for blocks
    
    for j = 1:N
        block = [exp(-1i * lambda(j) * x_pm)/2; exp(1i * lambda(j) * x_pm)/2];
        blocks{j} = block;  % Store each 1x2 block
    end
    
    T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
end

function T = That_blockdiag2(lambda, x_pm)
% Constructs T^{\pm} as a block-diagonal matrix
% Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]

    N = length(lambda);
    blocks = cell(1, N);  % Preallocate cell array for blocks
    
    for j = 1:N
        block = [exp(-1i * lambda(j) * x_pm)/(lambda(j)*2); -exp(1i * lambda(j) * x_pm)/(lambda(j)*2)];
        blocks{j} = block;  % Store each 1x2 block
    end
    
    T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
end

