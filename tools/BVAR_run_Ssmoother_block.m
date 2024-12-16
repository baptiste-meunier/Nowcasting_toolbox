function yMCMC = BVAR_run_Ssmoother_block(betaHat, sigmaHat, jump, lags,yQ)
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

% Added line initX(isnan(initX)) = 0; to take care of potential initial missing values

[~,n,nSim] = size(betaHat);

initX = zeros(lags*n,1);

for jt = lags:-1:1
    
%     initX((lags-jt)*n+1:(lags-jt+1)*n,:) = nanmean(yQ(1:lags,:),1)';
    initX((lags-jt)*n+1:(lags-jt+1)*n,:) = yQ(jt,:)';
    
end

initX(isnan(initX)) = 0;

% initV = eye(length(initX))*1e-7;

yInput = yQ(lags+1:end,:);
Tks = size(yInput,1);
yMCMC = NaN(Tks,n,floor(nSim/jump));




% % useful elements of measurement equation
C = zeros(n,n*lags); C(1:n,1:n) = eye(n);%maps first n states (current month) into observables
c1 = zeros(n,1); %constant in ME


parfor jm = 1:floor(size(betaHat,3)/jump)
    
    % transition equation
    A = zeros(n*lags);
    A(n+1:end,1:end-n) = eye(n*(lags-1));
    A(1:n,:) = betaHat(2:end,:,jm*jump)';% autoregressive coefficients
    % constant (in state eq)
    c2 = zeros(n*lags,1);
    c2(1:n,:)=betaHat(1,:,jm*jump)';
    
    % "shock" impact matrix
    Q = zeros(n*lags); 
    sTemp = sigmaHat(:,:,jm*jump);
    [V,E]=eig(sTemp);
    sTemp=V*abs(E)*V';
    Q(1:n,1:n) = sTemp;
    
    aPlus = NaN(n*lags,Tks);
    yPlus = NaN(n,Tks);
    initx = zeros(size(initX));
    
    for tt = 1:Tks
        
        randomBit = [mvnrnd(zeros(1,n),Q(1:n,1:n),1)';zeros(n*(lags-1),1)];
        aPlus(:,tt) =  A*initx+randomBit; 
        initx = aPlus(:,tt);
        yPlus(:,tt) = C*aPlus(:,tt)+c1;
        
    end
        
    yStar = yInput'-yPlus;
    initV = BVAR_lyapunov_symm(A, Q, 1, 1e-6);
    aHatStar = BVAR_runKF_DK(yStar, A, C, Q, zeros(n), initX, initV,c1,c2); %should we add a random draw to initX?
    aTilda = aHatStar+aPlus;
    yMCMC(:,:,jm)=aTilda(1:n,:)';
    
end

end