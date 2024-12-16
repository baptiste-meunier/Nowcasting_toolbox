function  [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = DFM_EMstep(y, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks)
% This function runs the expectation-maximization algorithm
%
% Code from Bańbura, M., and Modugno, M. (2014). “Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%
% Block identification from Delle Chiaie, S., Ferrara, L., and Giannone, D. 
% (2022). “Common factors of commodity prices”, Journal of Applied 
% Econometrics, 37(3), 461–476
%

[n,T] = size(y);
nM = n-nQ;

pC = size(R_mat,2);
ppC = max(p,pC);


n_b = size(blocks,2);

% Compute the (expected) sufficient statistics for a single Kalman filter sequence.

%Running the Kalman filter with the current estimates of the parameters
[Zsmooth, Vsmooth, VVsmooth, loglik] = DFM_runKF(y, A, C, Q, R, Z_0, V_0);

A_new = A;
Q_new = Q;
V_0_new = V_0;

for i = 1:n_b
    r_i = r(i);
    rp = r_i*p;
    rp1 = sum(r(1:i-1))*ppC;
    
    A_i = A(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC);
    Q_i = Q(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC);
    
    EZZ = Zsmooth(rp1+1:rp1+rp,2:end)*Zsmooth(rp1+1:rp1+rp,2:end)'...
        +sum(Vsmooth(rp1+1:rp1+rp,rp1+1:rp1+rp,2:end),3);                        %E(Z'Z)
    EZZ_BB = Zsmooth(rp1+1:rp1+rp,1:end-1)*Zsmooth(rp1+1:rp1+rp,1:end-1)'...
        +sum(Vsmooth(rp1+1:rp1+rp,rp1+1:rp1+rp,1:end-1),3); %E(Z(-1)'Z_(-1))
    EZZ_FB = Zsmooth(rp1+1:rp1+rp,2:end)*Zsmooth(rp1+1:rp1+rp,1:end-1)'...
        +sum(VVsmooth(rp1+1:rp1+rp,rp1+1:rp1+rp,:),3);%E(Z'Z_(-1))

    A_i(1:r_i,1:rp) = EZZ_FB(1:r_i,1:rp) * inv(EZZ_BB(1:rp,1:rp));
    Q_i(1:r_i,1:r_i) = (EZZ(1:r_i,1:r_i) - A_i(1:r_i,1:rp)*EZZ_FB(1:r_i,1:rp)') / T; %d
    
    A_new(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC) = A_i; 
    Q_new(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC) = Q_i;
    V_0_new(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC) = Vsmooth(rp1+1:rp1+r_i*ppC,rp1+1:rp1+r_i*ppC,1);
end

rp1 = sum(r)*ppC;
niM = sum(i_idio(1:nM)); 
% idiosyncratic
EZZ = diag(diag(Zsmooth(rp1+1:end,2:end)*Zsmooth(rp1+1:end,2:end)'))...
    +diag(diag(sum(Vsmooth(rp1+1:end,rp1+1:end,2:end),3)));                        %E(Z'Z)
EZZ_BB = diag(diag(Zsmooth(rp1+1:end,1:end-1)*Zsmooth(rp1+1:end,1:end-1)'))...
    +diag(diag(sum(Vsmooth(rp1+1:end,rp1+1:end,1:end-1),3))); %E(Z(-1)'Z_(-1))
EZZ_FB = diag(diag(Zsmooth(rp1+1:end,2:end)*Zsmooth(rp1+1:end,1:end-1)'))...
    +diag(diag(sum(VVsmooth(rp1+1:end,rp1+1:end,:),3)));%E(Z'Z_(-1)) 

A_i = EZZ_FB * diag(1./diag((EZZ_BB)));
Q_i = (EZZ - A_i*EZZ_FB') / T;

A_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = A_i(1:niM,1:niM); 
Q_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = Q_i(1:niM,1:niM);
V_0_new(rp1+1:rp1+niM,rp1+1:rp1+niM) = diag(diag(Vsmooth(rp1+1:rp1+niM,rp1+1:rp1+niM,1)));

Z_0 = Zsmooth(:,1); %zeros(size(Zsmooth,1),1); %


nanY = isnan(y);
y(nanY) = 0;

% LOADINGS
C_new = C;

% Blocks
bl = unique(blocks,'rows');
n_bl = size(bl,1);
bl_idxM = [];
bl_idxQ = [];
R_con = [];
q_con = [];
for i = 1:n_b
    bl_idxQ = [bl_idxQ repmat(bl(:,i),1,r(i)*ppC)];
    bl_idxM = [bl_idxM repmat(bl(:,i),1,r(i)) zeros(n_bl,r(i)*(ppC-1))];
    R_con = blkdiag(R_con, kron(R_mat,eye(r(i))));
    q_con = [q_con;zeros(r(i)*size(R_mat,1),1)];
end

bl_idxM = logical(bl_idxM);
bl_idxQ = logical(bl_idxQ);

%idio
i_idio_M = i_idio(1:nM);
n_idio_M = length(find(i_idio_M));
c_i_idio = cumsum(i_idio);

for i = 1:n_bl
    bl_i = bl(i,:);
    rs = sum(r(logical(bl_i)));
    idx_i = find(ismember(blocks,bl_i,'rows'));
    
    % MONTHLY
    idx_iM = idx_i(idx_i<nM+1);
    n_i = length(idx_iM);

    denom = zeros(n_i*rs,n_i*rs);
    nom = zeros(n_i,rs);

    i_idio_i = i_idio_M(idx_iM);
    i_idio_ii = c_i_idio(idx_iM);
    i_idio_ii = i_idio_ii(i_idio_i);
    for t=1:T
        nanYt = diag(~nanY(idx_iM,t));
        denom = denom + kron(Zsmooth(bl_idxM(i,:),t+1)*Zsmooth(bl_idxM(i,:),t+1)'...
            +Vsmooth(bl_idxM(i,:),bl_idxM(i,:),t+1),nanYt);
        nom = nom + y(idx_iM,t)*Zsmooth(bl_idxM(i,:),t+1)'...%here's the modification
            -nanYt(:,i_idio_i)*(Zsmooth(rp1+i_idio_ii,t+1)*Zsmooth(bl_idxM(i,:),t+1)'...
            +Vsmooth(rp1+i_idio_ii,bl_idxM(i,:),t+1));
    end
    vec_C = inv(denom)*nom(:);
    C_new(idx_iM,bl_idxM(i,:)) = reshape(vec_C,n_i,rs);

    % QUARTERLY
   idx_iQ = idx_i(idx_i>nM);
   rps = rs*ppC;
   
   R_con_i = R_con(:,bl_idxQ(i,:));
   q_con_i = q_con;
   no_c = ~(any(R_con_i,2));
   R_con_i(no_c,:) = [];
   q_con_i(no_c,:) = [];
    

   for j = idx_iQ'
       denom = zeros(rps,rps);
       nom = zeros(1,rps);
       idx_jQ = j-nM;
       i_idio_jQ = (rp1+n_idio_M+5*(idx_jQ-1)+1:rp1+n_idio_M+5*idx_jQ);
       V_0_new(i_idio_jQ,i_idio_jQ) = Vsmooth(i_idio_jQ,i_idio_jQ,1);
       A_new(i_idio_jQ(1),i_idio_jQ(1)) = A_i(i_idio_jQ(1)-rp1,i_idio_jQ(1)-rp1);
       Q_new(i_idio_jQ(1),i_idio_jQ(1)) = Q_i(i_idio_jQ(1)-rp1,i_idio_jQ(1)-rp1);

       for t=1:T
           nanYt = diag(~nanY(j,t));
           denom = denom + kron(Zsmooth(bl_idxQ(i,:),t+1)*Zsmooth(bl_idxQ(i,:),t+1)'...
                +Vsmooth(bl_idxQ(i,:),bl_idxQ(i,:),t+1),nanYt);
            nom = nom + y(j,t)*Zsmooth(bl_idxQ(i,:),t+1)';
            nom = nom -...
                nanYt*([1 2 3 2 1]*Zsmooth(i_idio_jQ,t+1)*Zsmooth(bl_idxQ(i,:),t+1)'+...
                [1 2 3 2 1]*Vsmooth(i_idio_jQ,bl_idxQ(i,:),t+1));

        end

        C_i = inv(denom)*nom';
        C_i_constr = C_i - inv(denom)*R_con_i'*inv(R_con_i*inv(denom)*R_con_i')*(R_con_i*C_i-q_con_i);
        C_new(j,bl_idxQ(i,:)) = C_i_constr;

    end


end


R_new = zeros(n,n);
for t=1:T
    nanYt = diag(~nanY(:,t));
    R_new = R_new + (y(:,t)-nanYt*C_new*Zsmooth(:,t+1))*(y(:,t)-nanYt*C_new*Zsmooth(:,t+1))'...
        +nanYt*C_new*Vsmooth(:,:,t+1)*C_new'*nanYt...
        +(eye(n)-nanYt)*R*(eye(n)-nanYt);
end

R_new = R_new/T;
RR = diag(R_new); %RR(RR<1e-2) = 1e-2;
RR(i_idio_M) = 1e-04;
RR(nM+1:end) = 1e-04;
R_new = diag(RR);

end