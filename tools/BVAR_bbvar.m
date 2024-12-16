function res = BVAR_bbvar(X_ext,lags,mSeries,stationary,lambda0,theta0,miu0,alpha0,thresh,max_iter,varargin)
% This script estimates a blocked BVAR (BBVAR)
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & 
% Sokol, A. (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

% Preliminaries
if ~isempty(varargin)
    
    doDensity = 1;
    nDraws = varargin{1};
    nSave = varargin{2};
    
else
    
    doDensity = 0;
    nDraws = 10000;
    
end

pos = zeros(1,size(X_ext,2));
pos(stationary) = 1;
[~,pos] = find([pos(mSeries) pos(mSeries) pos]==1);

% Prepare data
month1 = X_ext(1:3:end,mSeries);
month2 = X_ext(2:3:end,mSeries);
month3 = X_ext(3:3:end,:);
endEstimT = find(~isnan(sum(X_ext(3:3:end,:),2)),1,'last');
xQ = [month1 month2 month3];
xQEst = xQ(1:endEstimT,:);
optNaN.method = 1;
optNaN.k = 3; 
xQEst = BVAR_remNaNs_spline(xQEst,optNaN);

% Estimate
estimRes = BVAR_bvarGLP_alt(xQEst,lags,'mcmc',doDensity,'pos',pos,'MCMCfcast',0,'Fcast',0,'MNalpha',0,'Ndraws',2*nDraws,'MCMCconst',.14,'lambda0',lambda0,'theta0',theta0,'miu0',miu0,'alpha0',alpha0,'sur',0,'MNpsi',1,'thresh',thresh,'max_iter',max_iter);
% estimRes = bvarGLP_fixedhyp(xQEst,lags,'mcmc',doDensity,'pos',pos,'MCMCfcast',0,'Fcast',0,'MNalpha',1,'Ndraws',2*nDraws,'MCMCconst',1,'lambda0',lambda0,'theta0',theta0,'miu0',miu0,'alpha0',alpha0,'sur',1,'MNpsi',0,'hyperpriors',0);
betaHat = estimRes.postmax.betahat;
sigmaHat = estimRes.postmax.sigmahat;

% disp('Now running Kalman smoother')
% [X_sm,logLik] = run_Ksmoother_FM(betaHat, sigmaHat, lags,xQ);
[X_sm,logLik] = BVAR_run_Ksmoother_FM(betaHat, sigmaHat, lags,xQ(endEstimT-lags:end,:));
X_sm = [xQ(lags+1:endEstimT-1,:);X_sm];

if doDensity
    
    disp('Now running simulation smoother')
    yMCMC = NaN([size(xQ(lags+1:end,:)) nSave]);
    auxT = size(xQ(endEstimT-lags:end,:),1);
%     yMCMC = run_Ssmoother_block(estimRes.mcmc.beta,estimRes.mcmc.sigma,nDraws/nSave,lags,xQ);
    yMCMC(end-auxT+lags+1:end,:,:) = BVAR_run_Ssmoother_block(estimRes.mcmc.beta,estimRes.mcmc.sigma,nDraws/nSave,lags,xQ(endEstimT-lags:end,:));
    
    for ss = 1:size(yMCMC,3)
        
        yMCMC(1:end-auxT+lags,:,ss) = xQ(lags+1:endEstimT-1,:);
        
    end
    
    res.yMCMC = yMCMC;
    res.beta = estimRes.mcmc.beta;
    res.sigma = estimRes.mcmc.sigma;
    
end

res.X_sm = X_sm;
res.logLik = logLik;
res.xQ = xQ;
res.betaHat = betaHat;
res.sigmaHat = sigmaHat;
res.pos = pos;
res.lambda = estimRes.postmax.lambda;
res.theta = estimRes.postmax.theta;
res.miu = estimRes.postmax.miu;
res.alpha = estimRes.postmax.alpha;
end