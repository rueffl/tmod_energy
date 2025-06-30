% function [S_tilde] = get_tdepStilde(xm, xp, delta, k_tr, vr, omega, Omega, F, lambdas)
% %GET_STILDE computes the matrix connecting [\alpha_i;\beta_i] to [\alpha_{i+1};\beta_{i+1}]
% %   xm:     left boundary point of the resonator
% %   xp:     right boundary point of the resonator
% %   delta:  contrast parameter
% %   k_tr:   truncation parameter
% %   vr:     wave speed inside the resonators
% %   omega:  operating frequency
% %   Omega:  modulation frequency
% %   F:      matrix of eigenvectors
% %   lambdas:eigenvalues of the resonator
% 
%     kns = (omega+[-k_tr:k_tr].*Omega)./vr;
%     G = inv(F);
% 
%     Eplus_left = diag(exp(1i*kns.*xm)); Eminus_left = diag(exp(-1i*kns.*xm));
%     Eplus_right = diag(exp(1i*kns.*xp)); Eminus_right = diag(exp(-1i*kns.*xp));
%     K = diag(kns);
%     Tplus = T_blockdiag(lambdas, xp); Tminus = T_blockdiag(lambdas, xm);
%     Tpluslambda = T_blockdiag2(lambdas, xp); Tminuslambda = T_blockdiag2(lambdas, xm);
% 
%     Mminus = [G*Eplus_left, G*Eminus_left; delta.*G*K*Eplus_left, -delta.*G*K*Eminus_left];
%     Tauminus = [Tminus; Tminuslambda]; 
%     invTauminus = inv(Tauminus);
%     Mplus = [Eplus_right, Eminus_right; delta.*K*Eplus_right, -delta.*K*Eminus_right]; 
%     invMplus = inv(Mplus);
%     Tauplus = [F*Tplus; F*Tpluslambda];
% 
%     S_tilde = invMplus*Tauplus*invTauminus*Mminus;
% 
% end
% 
% 
% function T = T_blockdiag(lambda, x_pm)
% % Constructs T^{\pm} as a block-diagonal matrix
% % Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]
% 
%     N = length(lambda);
%     blocks = cell(1, N);  % Preallocate cell array for blocks
%     
%     for j = 1:N
%         block = [exp(1i * lambda(j) * x_pm), exp(-1i * lambda(j) * x_pm)];
%         blocks{j} = block;  % Store each 1x2 block
%     end
%     
%     T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
% end
% 
% function T = T_blockdiag2(lambda, x_pm)
% % Constructs T^{\pm} as a block-diagonal matrix
% % Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]
% 
%     N = length(lambda);
%     blocks = cell(1, N);  % Preallocate cell array for blocks
%     
%     for j = 1:N
%         block = [lambda(j)*exp(1i * lambda(j) * x_pm), -lambda(j)*exp(-1i * lambda(j) * x_pm)];
%         blocks{j} = block;  % Store each 1x2 block
%     end
%     
%     T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
% end
% 
% 
% function T = That_blockdiag(lambda, x_pm)
% % Constructs T^{\pm} as a block-diagonal matrix
% % Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]
% 
%     N = length(lambda);
%     blocks = cell(1, N);  % Preallocate cell array for blocks
%     
%     for j = 1:N
%         block = [exp(-1i * lambda(j) * x_pm)/2; exp(1i * lambda(j) * x_pm)/2];
%         blocks{j} = block;  % Store each 1x2 block
%     end
%     
%     T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
% end
% 
% function T = That_blockdiag2(lambda, x_pm)
% % Constructs T^{\pm} as a block-diagonal matrix
% % Each block is a 1x2 row: [exp(i*lambda_j*x_pm), exp(-i*lambda_j*x_pm)]
% 
%     N = length(lambda);
%     blocks = cell(1, N);  % Preallocate cell array for blocks
%     
%     for j = 1:N
%         block = [exp(-1i * lambda(j) * x_pm)/(lambda(j)*2); -exp(1i * lambda(j) * x_pm)/(lambda(j)*2)];
%         blocks{j} = block;  % Store each 1x2 block
%     end
%     
%     T = blkdiag(blocks{:});  % Combine into block-diagonal matrix
% end


function [Stilde] = get_tdepStilde(xm, xp, delta, k, F, lambdas)
%GET_TDEPS computes the transfer matrix corresponding to one resonator
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

    Stilde = [Stilde11, Stilde12; Stilde21, Stilde22];

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

