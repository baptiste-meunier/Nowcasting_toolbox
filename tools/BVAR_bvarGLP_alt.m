function r = BVAR_bvarGLP_alt(y,lags,varargin)
% This function estimates the BVAR of Giannone, D., Lenza, M., and 
% Primiceri, G. E. (2015). “Prior Selection for Vector Autoregressions”, 
% The Review of Economics and Statistics, 97(2), 436–451
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% y:        data matrix
% lags:     number of lags in the VAR
%
% 


%% set BVAR priors (several options available, see setpriors.m)
%  if varargin=[] --> default settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>2
    for ii=1:2:size(varargin,2)
        eval([varargin{ii},'= varargin{ii+1};'])
    end
end

%------------------------------
if ~exist('hyperpriors','var')
    hyperpriors=1;
end
r.setpriors.hyperpriors=hyperpriors;

%------------------------------
if ~exist('Vc','var')
    Vc=10e6;
end
r.setpriors.Vc=Vc;

%------------------------------
if ~exist('pos','var')
    pos=[];
end
r.setpriors.pos=pos;

%------------------------------
if ~exist('MNalpha','var')
    mn.alpha=0;
else
    mn.alpha=MNalpha;
end
r.setpriors.MNalpha=mn.alpha;

%------------------------------
if ~exist('MNpsi','var')
    mn.psi=1;
else
    mn.psi=MNpsi;
end
r.setpriors.MNpsi=mn.psi;

%------------------------------
if ~exist('sur','var')
    sur=1;
end
r.setpriors.sur=sur;

%------------------------------
if ~exist('noc','var')
    noc=1;
end
r.setpriors.noc=noc;

%------------------------------
if ~exist('Fcast','var')
    Fcast=0;
end
r.setpriors.Fcast=Fcast;

%------------------------------
if ~exist('hz','var')
    hz=1:8; 
else
    hz=1:hz;
end
r.setpriors.hz=hz;

%------------------------------
if ~exist('mcmc','var')
    mcmc=0; 
end
r.setpriors.mcmc=mcmc;

%------------------------------
if ~exist('Ndraws','var')
    M=20000;
else
    M=Ndraws;
end
r.setpriors.Ndraws=M;

%------------------------------
if ~exist('Ndrawsdiscard','var')
    N=round(M/2);
else
    N=Ndrawsdiscard;
end
r.setpriors.Ndrawsdiscard=N;

%------------------------------
if ~exist('MCMCconst','var')
    const=1; 
else
    const=MCMCconst;
end
r.setpriors.MCMCconst=const;

%------------------------------
if ~exist('MCMCfcast','var')
    MCMCfcast=1; 
end
r.setpriors.MCMCfcast=MCMCfcast;

%------------------------------
if ~exist('MCMCstorecoeff','var')
    MCMCstorecoeff=1; 
end
r.setpriors.MCMCstorecoeff=MCMCstorecoeff;
%------------------------------

%% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of the hyperpriors, if choosen
if hyperpriors==1
    
    mode.lambda=.2;     % hyperpriors modes
    mode.miu=1;
    mode.theta=1;
    
    sd.lambda=.4;       % hyperpriors std
    sd.miu=1;
    sd.theta=1;
    
    scalePSI=(0.02^2);    % scale and shape of the IG on psi/(d-n-1). divide by 16 because we are not using 4logs    
    
    priorcoef.lambda = BVAR_GammaCoef(mode.lambda,sd.lambda,0);  % coefficients of hyperpriors
    priorcoef.miu = BVAR_GammaCoef(mode.miu,sd.miu,0);
    priorcoef.theta = BVAR_GammaCoef(mode.theta,sd.theta,0);
    priorcoef.alpha.PSI = scalePSI;
    priorcoef.beta.PSI = scalePSI;
    
else
    priorcoef=[];
end

% bounds for maximization medium_real blocking
MIN.lambda = 0.0001;
MIN.alpha = 0.1;
MIN.theta = 0.0001;
MIN.miu = .0001;%0.0001;
MAX.lambda=5;
MAX.miu=50;%50;
MAX.theta=5;
MAX.alpha=5;
% MIN.lambda = 0.0001;
% MIN.alpha = 0.1;
% MIN.theta = 0.0001;
% MIN.miu = 0.0001;
% MAX.lambda=5;
% MAX.miu=50;
% MAX.theta=50;
% MAX.alpha=5;
% 



%% data matrix manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions
[TT,n]=size(y);
k=n*lags+1;         % # coefficients for each equation


% constructing the matrix of regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(TT,k);
x(:,1)=1;
for i=1:lags
    x(:,1+(i-1)*n+1:1+i*n) = BVAR_lag(y,i);
end

x=x(lags+1:end,:);
y=y(lags+1:end,:);
y0=mean(y(1:lags,:),1);
[T,n]=size(y);

r.postmax.y = y;
r.postmax.x = x;

% MN prior mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(k,n);
diagb=ones(n,1);
diagb(pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
b(2:n+1,:)=diag(diagb);


%% starting values for the minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('lambda0','var')
    lambda0=.2;     % std of MN prior
end
if ~exist('theta0','var')
    theta0=1;       % std of sur prior    
end
if ~exist('miu0','var')
    miu0=1;         % std of noc prior
end
if ~exist('alpha0','var')
   alpha0=2;       % lag-decaying parameter of the MN prior
end


% residual variance of AR(1) for each variable
SS=zeros(n,1);
for i=1:n
    ar1 = BVAR_ols1(y(2:end,i),[ones(T-1,1),y(1:end-1,i)]);
    SS(i)=ar1.sig2hatols;
end
MIN.psi=SS./100;
MAX.psi=SS.*100;
psi0=SS;

if mn.psi==1
    inpsi=-log((MAX.psi-psi0)./(psi0-MIN.psi));
elseif mn.psi==0
    inpsi=[];
end

if mn.alpha==1
    inalpha=-log((MAX.alpha-alpha0)/(alpha0-MIN.alpha));
elseif mn.alpha==0
    inalpha=[];
end

if sur==1
    intheta=-log((MAX.theta-theta0)/(theta0-MIN.theta));
elseif sur==0
    intheta=[];
end

if noc==1
    inmiu=-log((MAX.miu-miu0)/(miu0-MIN.miu));
elseif noc==0
    inmiu=[];
end


x0=[-log((MAX.lambda-lambda0)/(lambda0-MIN.lambda));inpsi;intheta;inmiu;inalpha];

H0=10*eye(length(x0));          % initial guess for the inverse Hessian



%% maximization of the posterior of the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changed max. number of iterations to 1000 to avoid too long convergence
[~,xh,~,~,itct,~,~] = BVAR_csminwel('BVAR_logMLVAR_formin',x0,H0,[],thresh,max_iter,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,hyperpriors,priorcoef);


%% output of the maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR coefficients and residuals at the posterior model
[fh,r.postmax.betahat,r.postmax.sigmahat] = BVAR_logMLVAR_formin(xh,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,hyperpriors,priorcoef);

r.lags = lags;                      % # lags
r.postmax.itct=itct;                % #iteration before reaching maximum
r.postmax.SSar1=SS;                 % residual variance of AR(1) for each variable
r.postmax.logPost=-fh;              % value of the posterior of the hyperparameters at the peak
r.postmax.lambda=MIN.lambda+(MAX.lambda-MIN.lambda)/(1+exp(-xh(1)));    % std of MN prior at the peak

r.postmax.theta=theta0;
r.postmax.miu=miu0;

r.postmax.residuals = r.postmax.y-r.postmax.x*r.postmax.betahat;

if mn.psi==1
    % diagonal elements of the scale matrix of the IW prior on the residual variance
    r.postmax.psi=MIN.psi+(MAX.psi-MIN.psi)./(1+exp(-xh(2:n+1)));
    if sur==1
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-xh(n+2)));
        if noc==1
            % std of noc prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(n+3)));
        end
    elseif sur==0
        if noc==1
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(n+2)));
        end
    end
elseif mn.psi==0
    r.postmax.psi=SS;
    if sur==1
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-xh(2)));
        if noc==1
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(3)));
        end
    elseif sur==0
        if noc==1
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(2)));
        end
    end
end

if mn.alpha==0
    r.postmax.alpha=2;
elseif mn.alpha==1
    % Lag-decaying parameter of the MN prior
    r.postmax.alpha=MIN.alpha+(MAX.alpha-MIN.alpha)/(1+exp(-xh(end)));
end


%% forecasts at the posterior mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Fcast==1
    Y=[y;zeros(hz(end),n)];
    for tau=1:max(hz)
        xT=[1;reshape(Y(T+tau-1:-1:T+tau-lags,:)',k-1,1)]';
        Y(T+tau,:)=xT*r.postmax.betahat;
    end
    r.postmax.forecast=Y(T+hz,:);
end


%% MCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mcmc==1
    
%     % Jacobian of the transformation of the hyperparameters that has been
%     % used for the constrained maximization
%     JJ=exp(xh)./(1+exp(xh)).^2;
%     JJ(1)=(MAX.lambda-MIN.lambda)*JJ(1);
%     
% %     MAX.lambda = 2*MAX.lambda;
% %     MIN.theta = .1*MIN.theta;
%     
%     if mn.psi==1;
%         JJ(2:n+1)=(MAX.psi-MIN.psi).*JJ(2:n+1);
%         if sur==1;
%             JJ(n+2)=(MAX.theta-MIN.theta)*JJ(n+2);
%             if noc==1;
%                 JJ(n+3)=(MAX.miu-MIN.miu)*JJ(n+3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 JJ(n+2)=(MAX.miu-MIN.miu)*JJ(n+2);
%             end
%         end
%     elseif mn.psi==0;
%         if sur==1;
%             JJ(2)=(MAX.theta-MIN.theta)*JJ(2);
%             if noc==1;
%                 JJ(3)=(MAX.miu-MIN.miu)*JJ(3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 JJ(2)=(MAX.miu-MIN.miu)*JJ(2);
%             end
%         end
%     end
%     
%     if mn.alpha==1;
%         JJ(end)=(MAX.alpha-MIN.alpha)*JJ(end);
%     end
%     
%     JJ=diag(JJ);
%     HH=JJ*H*JJ';
%     
%     % regularizing the Hessian (making sure it is positive definite)
%     [V,E]=eig(HH);
%     HH=real(V*abs(E)*V'); %added real
    
    % recovering the posterior mode
    if mn.psi==1
        modepsi=r.postmax.psi;
    elseif mn.psi==0
        modepsi=[];
    end
    
    if mn.alpha==1
        modealpha=r.postmax.alpha;
    elseif mn.alpha==0
        modealpha=[];
    end
    
    if sur==1
        modetheta=r.postmax.theta;
    elseif sur==0
        modetheta=[];
    end
    
    if noc==1
        modemiu=r.postmax.miu;
    elseif noc==0
        modemiu=[];
    end
    
    postmode=[r.postmax.lambda;modepsi;modetheta;modemiu;modealpha];
    
    fun = @(par) BVAR_logMLVAR_formcmc(par,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,0,hyperpriors,priorcoef);%,Tcovid);
    Hess = hessian_new(fun,postmode);
    HH = -inv2(Hess);
    [V,E]=eig(HH);
    E = abs(E);
    HH=real((V*E)*V'); %added real

%     E = diag(E);
%     E(E<0) = 1e-10;
%     HH = real(V*(V'.*E));

    r.postmax.HH = HH;
    

    
    % starting value of the Metropolis algorithm
    disp('Now starting MH algorithm')
    dLength = length(xh);
    P=zeros(dLength,M);
%     cholHess = const*chol(HH,'lower');
    logMLold=-10e15;
    while logMLold==-10e15
        P(:,1)= mvnrnd(postmode',HH*const^2,1)';
        [logMLold,~,~] = BVAR_logMLVAR_formcmc(P(:,1),y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
            max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
    end
    
    % matrix to store the draws of the VAR coefficients if MCMCstorecoeff is on
    if MCMCstorecoeff==1
        r.mcmc.beta  = zeros(k,n,M-N);
        r.mcmc.sigma = zeros(n,n,M-N);
    end
    
    % matrix to store the forecasts if MCMCfcast is on
    if MCMCfcast==1
        r.mcmc.Dforecast=zeros(length(hz),n,M-N);
    end
    
    % Metropolis iterations
    count=0;
    for i=2:M
        if i==1000*floor(.001*i)
                        disp(['Now running the ',num2str(i),'th mcmc iteration (out of ',num2str(M),')'])
        end
        % draw candidate value
        P(:,i)=mvnrnd(P(:,ii-1)',HH*const^2,1)';
        [logMLnew,betadrawnew,sigmadrawnew] = BVAR_logMLVAR_formcmc(P(:,i),y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
            max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
        
        if logMLnew>logMLold
            logMLold=logMLnew;
            count=count+1;
        else
            if rand(1)<exp(logMLnew-logMLold)
                logMLold=logMLnew;
                count=count+1;
            else
                P(:,i)=P(:,i-1);
                % if MCMCfcast is on, take a new draw of the VAR coefficients
                % with the old hyperparameters if have rejected the new ones
                % (the speed of this step could be probably improved)
                if MCMCfcast==1 || MCMCstorecoeff==1
                    [~,betadrawnew,sigmadrawnew] = BVAR_logMLVAR_formcmc(P(:,i),y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
                        max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
                end
            end
        end
        
        % stores draws of VAR coefficients if MCMCstorecoeff is on
        if i>N && MCMCstorecoeff==1     
            r.mcmc.beta(:,:,i-N) = betadrawnew;  
            r.mcmc.sigma(:,:,i-N) = sigmadrawnew;   
        end
        
        % produce and store the forecasts if MCMCfcast is on
        if i>N && MCMCfcast==1
            Y=[y;zeros(hz(end),n)];
            for tau=1:max(hz)
                xT=[1;reshape(Y(T+tau-1:-1:T+tau-lags,:)',k-1,1)]';
                %Y(T+tau,:)=xT*betadrawnew;
                Y(T+tau,:)=xT*betadrawnew+mvnrnd(zeros(1,n),sigmadrawnew);
            end
            r.mcmc.Dforecast(:,:,i-N)=Y(T+hz,:);
            %r.IRF(:,:,:,i) = IRF(betadrawnew,sigmadrawnew);
        end
        
    end
    
    % store the draws of the hyperparameters
    r.mcmc.lambda=P(1,N+1:end);   % std MN prior
    
    if mn.psi==1
        % diagonal elements of the scale matrix of the IW prior on the residual variance
        r.mcmc.PSI=P(2:n+1,N+1:end)';
        if sur==1
            % std of sur prior
            r.mcmc.theta=P(n+2,N+1:end);
            if noc==1
                % std of noc prior
                r.mcmc.miu=P(n+3,N+1:end);
            end
        elseif sur==0
            if noc==1
                % std of noc prior
                r.mcmc.miu=P(n+2,N+1:end);
            end
        end
    elseif mn.psi==0
        if sur==1
            % std of sur prior
            r.mcmc.theta=P(2,N+1:end);
            if noc==1
                % std of noc prior
                r.mcmc.miu=P(3,N+1:end);
            end
        elseif sur==0
            if noc==1
                % std of noc prior
                r.mcmc.miu=P(2,N+1:end);
            end
        end
    end
    
    if mn.alpha==1
        % Lag-decaying parameter of the MN prior
        r.mcmc.alpha=P(end,N+1:end);
    end
    
    % acceptance rate
    r.mcmc.ACCrate=mean(r.mcmc.lambda(2:end)~=r.mcmc.lambda(1:end-1));
     
end