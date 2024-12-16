function [yHat,logLik] = BVAR_run_Ksmoother_FM(betaHat, sigmaHat, lags,yQ)
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

% Added line initX(isnan(initX)) = 0; to take care of potential initial missing values
[k,n] = size(betaHat);
% transition equation
A = zeros(n*lags);
A(1:n,:) = betaHat(2:end,:)';% autoregressive coefficients
A(n+1:end,1:end-n) = eye(n*(lags-1));
c2 = zeros(n*lags,1);
c2(1:n,:)=betaHat(1,:)';% constant (in state eq)

% % useful elements of measurement equation
C = zeros(n,n*lags); C(1:n,1:n) = eye(n);%maps first n states (current month) into observables
c1 = zeros(n,1); %constant in ME

% "shock" impact matrix 
Q = zeros(n*lags);
Q(1:n,1:n) = sigmaHat;


initX = zeros(lags*n,1);

for jt = lags:-1:1
    
%     initX((lags-jt)*n+1:(lags-jt+1)*n,:) = nanmean(yQ(1:lags,:),1)';
    initX((lags-jt)*n+1:(lags-jt+1)*n,:) = yQ(jt,:)';
    
end

initX(isnan(initX)) = 0;

% initV = eye(length(initX))*1e-7;
initV = BVAR_lyapunov_symm(A, Q, 1, 1e-6);

yInput = yQ(lags+1:end,:);
[xSmooth,logLik] = BVAR_runKF_DK(yInput',A, C, Q, zeros(n), initX, initV,c1,c2);
yHat = xSmooth(1:n,:)';

end