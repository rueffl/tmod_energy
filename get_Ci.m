function Ci = get_Ci(k_tr,i,w,Omega,rs,ks,vr)
%GET_CI Builds the matrix C_i to compute lambda and f
%   k_tr:   truncation parameter
%   i:      index of the resonator of interest
%   Omega:  modulation frequency
%   w:      subwavelength frequency
%   rs:     Fourier coefficients of 1/\rho_i(t)
%   ks:     Fourier coefficients of 1/\kappa_i(t)
%   vr:     wave speed insid the resonators

%   build -1 subdigonal of B
    b_1 = zeros(2*k_tr,1);
    j_ind = 1;
    for j = -k_tr:k_tr-1
        b_1(j_ind) = make_gamma(j,i,-1,w,Omega,ks,vr);
        j_ind = j_ind + 1;
    end

%   build main diagonal of B
    b0 = zeros(2*k_tr+1,1);
    j_ind = 1;
    for j = -k_tr:k_tr
        b0(j+k_tr+1) = make_gamma(j,i,0,w,Omega,ks,vr);
        j_ind = j_ind + 1;
    end

%   build 1 subdigonal of B
    b1 = zeros(2*k_tr,1);
    j_ind = 1;
    for j = -k_tr+1:k_tr
        b1(j_ind)= make_gamma(j,i,1,w,Omega,ks,vr);
        j_ind = j_ind + 1;
    end

%   build B (eq28)
    Bi = diag(b_1,-1)+diag(b0)+diag(b1,1);

%   build A (eq28)
    a_1 = ones(2*k_tr,1)*rs(i,1);
    a0 = ones(2*k_tr+1,1)*rs(i,2);
    a1 = ones(2*k_tr,1)*rs(i,3);
    Ai = diag(a_1,-1)+diag(a0)+diag(a1,1);

%   compute C
    Ci = Ai\Bi;
    %Ci = inv(Ai)*Bi;

end

function gamma = make_gamma(n,i,m,w,Omega,ks,vr)
    gamma = (w+(n-m)*Omega)*(w+n*Omega)*ks(i,m+2)/(vr^2);
%     gamma = (w+(n-m)*Omega)/(w+n*Omega)*ks(i,m+2)*((w+n*Omega)/vr)^2;
end