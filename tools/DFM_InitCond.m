function [ A, C, Q, R, initZ, initV] = DFM_InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio)
% This function sets initial conditions for the estimation of the DFM
%
% Code from Bańbura, M., and Modugno, M. (2014). “Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%
% Block identification from Delle Chiaie, S., Ferrara, L., and Giannone, D. 
% (2022). “Common factors of commodity prices”, Journal of Applied 
% Econometrics, 37(3), 461–476
%

pC = size(R_mat,2);
ppC = max(p,pC);
n_b = size(blocks,2);

OPTS.disp=0;

[xBal,indNaN] = DFM_remNaNs_spline(xNaN,optNaN);
[T,N] = size(xBal);
NM = N-nQ;

xNaN = xBal;
xNaN(indNaN) = nan;
C = [];
A = [];
Q = [];
initV = [];

res = xBal;
resNaN = xNaN;
indNaN(1:pC-1,:) = true;
for i = 1:n_b
    r_i = r(i);
    %--------------------------------------------------------------------------
    % Observation equation
    %--------------------------------------------------------------------------
    C_i = zeros(N,r_i*ppC);
    idx_i = find(blocks(:,i));
    idx_iM = idx_i(idx_i<NM+1);
    idx_iQ = idx_i(idx_i>NM);
    [ v, d ] = eigs(cov(res(:,idx_iM)),r_i,'lm',OPTS);
    
    if v == 1  % Relevant for local factors that load on only 2 variables 
        v = 0.99;
    end
    
    C_i(idx_iM,1:r_i) = v;
    f = res(:,idx_iM)*v;
    F = [];
    for kk = 0:max(p+1,pC)-1
        F = [F f(pC-kk:end-kk,:)];
    end
    Rcon_i = kron(R_mat,eye(r_i));
    q_i = kron(q,zeros(r_i,1));
    ff = F(:,1:r_i*pC);
    for j = idx_iQ'
        xx_j = resNaN(pC:end,j);
        if sum(~isnan(xx_j)) < size(ff,2)+2
            xx_j = res(pC:end,j);
        end
        ff_j = ff(~isnan(xx_j),:);
        xx_j = xx_j(~isnan(xx_j));
        iff_j = inv(ff_j'*ff_j);
        Cc = iff_j*ff_j'*xx_j;
        Cc = Cc - iff_j*Rcon_i'*inv(Rcon_i*iff_j*Rcon_i')*(Rcon_i*Cc-q_i);
        C_i(j,1:pC*r_i)=Cc';
    end
    ff = [zeros(pC-1,pC*r_i);ff];
    res = res - ff*C_i';
    resNaN = res;
    resNaN(indNaN) = nan;
    C = [C C_i];

    %--------------------------------------------------------------------------
    % Transition equation
    %--------------------------------------------------------------------------
    z = F(:,1:r_i);
    Z = F(:,r_i+1:r_i*(p+1));
    A_i = zeros(r_i*ppC,r_i*ppC)';
    A_temp = inv(Z'*Z)*Z'*z;
    A_i(1:r_i,1:r_i*p) = A_temp';
    A_i(r_i+1:end,1:r_i*(ppC-1)) = eye(r_i*(ppC-1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_i = zeros(ppC*r_i,ppC*r_i);
    e = z  - Z*A_temp;         % VAR residuals
    Q_i(1:r_i,1:r_i) = cov(e); % VAR covariance matrix

    initV_i = reshape(inv(eye((r_i*ppC)^2)-kron(A_i,A_i))*Q_i(:),r_i*ppC,r_i*ppC);

    A = blkdiag(A,A_i);
    Q = blkdiag(Q,Q_i);
    initV = blkdiag(initV,initV_i);


end


R = diag(nanvar(resNaN));

eyeN = eye(N);
eyeN(:,~i_idio) = [];
% Initial conditions
C=[C eyeN];

ii_idio = find(i_idio);
n_idio = length(ii_idio);
B = zeros(n_idio);
S = zeros(n_idio);

for i = 1:n_idio;
    R(ii_idio(i),ii_idio(i)) = 1e-04;

    res_i = resNaN(:,ii_idio(i));
    % number of leading zeros
    leadZero = max( find( (1:T)' == cumsum(isnan(res_i)) ) );
    endZero = max( find( (1:T)' == cumsum(isnan(res_i(end:-1:1))) ) );

    res_i = res(:,ii_idio(i));
    res_i(end-endZero:endZero) = [];
    res_i(1:leadZero) = [];

    BM(i,i) = inv(res_i(1:end-1)'*res_i(1:end-1))*res_i(1:end-1)'*res_i(2:end,:);
    SM(i,i) = cov(res_i(2:end)-res_i(1:end-1)*B(i,i));
end
initViM=[];
if n_idio~=0
    initViM = diag(1./diag(eye(size(BM,1))-BM.^2)).*SM;
else
    BM=[];
    SM=[];
end


C = [C [zeros(NM,5*nQ);kron(eye(nQ),[1 2 3 2 1])]];
Rdiag = diag(R);
sig_e = Rdiag(NM+1:N)/19;
Rdiag(NM+1:N) = 1e-04;
R = diag(Rdiag);


rho0 = 0.1;
BQ = kron(eye(nQ),[[rho0 zeros(1,4)];[eye(4),zeros(4,1)]]);
temp = zeros(5);
temp(1,1) = 1;
SQ = kron(diag((1-rho0^2)*sig_e),temp);

%initViQ = eye((5*NQ));
initViQ = reshape(inv(eye((5*nQ)^2)-kron(BQ,BQ))*SQ(:),5*nQ,5*nQ);

A = blkdiag(A, BM, BQ);
Q = blkdiag(Q, SM, SQ);

% Initial conditions
initZ = zeros(size(A,1),1); 
initV = blkdiag(initV, initViM, initViQ);

end