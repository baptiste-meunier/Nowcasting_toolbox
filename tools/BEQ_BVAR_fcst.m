function X_new = BEQ_BVAR_fcst(X_in,Par_BVAR)
% This scripts run a BVAR for forecast
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - X_in [matrix] = input data
% - Par.BVAR [structure] = parameters of the BVAR
%       o nDraw [scalar] = number of draws (0 means a point forecast)
%       o p [scalar] = number of lags
%       o lambda [scalar] = overall shrinkage parameters
%
% OUTPUT
% - X_new [matrix] = interpolated data
%

% Get parameters of the BVAR
lambda = Par_BVAR.lambda;
KK     = Par_BVAR.KK;
iRW    = Par_BVAR.iRW;
p      = Par_BVAR.lags;
nDraws = Par_BVAR.nDraw;

% Standardising the data (for computational reasons)
Tall = size(X_in,1);
Mstd = mean(X_in,1,'omitnan');
Sstd = std(X_in,1,'omitnan');
X_in = (X_in-repmat(Mstd,Tall,1))./repmat(Sstd,Tall,1);

% Position of full data
nNaNsE   =  find(all(isfinite(X_in(end:-1:1,:)),2),1,'first')-1; % Number of NaN at the end of dataset
nNaNsS   = find(all(isfinite(X_in),2),1,'first'); % Position of first NaN

% Correct datasets for NaN
Xub   = X_in(end-nNaNsE+1:end,:);   % unbalanced end data
Xest  = X_in(nNaNsS:end-nNaNsE,:);  % balanced data for estimation
Xinit = X_in(1:end-nNaNsE,:);       % initial data


% =========================================================================
% Matrix of regressors
% =========================================================================

[T,N] = size(Xest);   % The size of the panel

Z = nan(T-p,p*N+1); % matrix of regressors (cutting p lags, and put all variables with their lags)
for jp = 1:p % looping on lags
    Z(:,N*(jp-1)+1:N*jp) = Xest(p+1-jp:end-jp,:);
end
Z(:,end) = ones(size(Z(:,1))); % adding a constant
Y        = Xest(p+1:end,:); % target data


% =========================================================================
% Estimation
% =========================================================================

Np1 = N*p+1; % size of beta

if  isinf(lambda) % OLS case

    Ypr = zeros(0,N);
    Zpr = zeros(0,Np1);
    Tpr = Np1-3;

else
    
    Vc = 10e6;    
    
    SS0 = nan(1,N);

    for jn = 1:N
        [~,~,Car] = BEQ_arfit(Xest(:,jn),p,p);
        SS0(jn) = sqrt(Car);
    end
    
    
    % Dummies
    % Construct dummy for litterman prior
    Yrw = [diag(SS0.*iRW);zeros(N*(p-1),N)]/lambda;
    Zrw = [diag(kron(1:p,SS0)) zeros(N*p,1) ]/lambda;
    
    % Construct dummy for the constant
    Ycs = Vc^-1*zeros(1,N);
    Zcs = Vc^-1*[zeros(1,N*p) 1];
    
    % Construct dummy for the sum of the coefficients 
    if KK>0
        MM0  = mean(Xest(1:p,:));
        mu   = lambda*KK;
        Ynoc = diag(MM0.*iRW)/mu;
        Znoc = [kron(ones(1,p),diag(MM0.*iRW)) zeros(N,1)]/mu;
    else
        Ynoc = zeros(0,N);
        Znoc = zeros(0,N*p+1);
    end
    
    % Construct dummies for prior on covariance matrix of residual;
    Ycv = diag(SS0);
    Zcv = zeros(N,Np1);
    
    % Put together all the information
    Ypr = [Yrw; Ynoc; Ycv; Ycs];
    Zpr = [Zrw; Znoc; Zcv; Zcs];
    
    Tpr = size(Ypr,1);
       
end


% =========================================================================
% Posterior
% =========================================================================

ZZinv = (Zpr'*Zpr + Z'*Z)\eye(Np1);
ZY    = Zpr'*Ypr + Z'*Y;

beta = ZZinv*ZY;

e  = [Y; Ypr]-[Z; Zpr]*beta;
S  = e'*e/(T+Tpr-Np1+2);


%==========================================================================
% State space representation
%==========================================================================

% Transition equation
AA                    = zeros(N*p);
AA(1:N,1:N*p)         = beta(1:end-1,:)'; % constant comes through initial conditions
AA(N+1:N*p,1:N*(p-1)) = eye(N*(p-1));
c2                    = [beta(end,:)';zeros(N*(p-1),1)];

% Measurement equation
CC = zeros(N,N*p); CC(:,1:N) = eye(N);
QQ = zeros(N*p); QQ(1:N,1:N) = S;
c1 = zeros(N,1);

% Initialisation of the Kalman filter
initx = lagmatrix(Xinit,0:p-1);
initx = initx(end,:)';
initV = eye(length(initx))*1e-16;


%==========================================================================
% CONDITIONAL FORECASTS
%==========================================================================

Tub = size(Xub,1);

Par_BVAR.C   = CC;
Par_BVAR.R   = diag(ones(1,N)*1e-16);
Par_BVAR.A   = AA;
Par_BVAR.Q   = QQ;
Par_BVAR.Z_0 = initx;
Par_BVAR.V_0 = initV;
Par_BVAR.c1  = c1;
Par_BVAR.c2  = c2;

if nDraws == 0 % point forecast

    % Kalman filter and smoother
    Res = BEQ_RUN_KF(Xub,Par_BVAR);
    X_new = Res.X_sm;
    X_new = [Xinit;X_new];

else
    
    Par_BVAR.Z_0 = zeros(size(initx));
    Par_BVAR.c1  = zeros(N,1);
    Par_BVAR.c2  = zeros(size(initx));
    
    X_new = nan(Tall,N,nDraws);
    % Durbin and Koopman simulation smoother
    for kg = 1:nDraws
        aplus    = nan(N*p,Tub);
        Xplus    = nan(N,Tub);
        initx_kg = initx;
        for t = 1:Tub
            aplus(:,t) =  AA*initx_kg+...
                [mvnrnd(zeros(N,1),S,1)';zeros(N*(p-1),1)]+c2;
            initx_kg   = aplus(:,t);
            Xplus(:,t) = CC*aplus(:,t)+c1;

        end
        Xstar       = Xub-Xplus';
        Res         = BEQ_RUN_KF(Xstar,Par_BVAR);
        ahatstar    = Res.AmT;
        atilda      = ahatstar+aplus;
        X_new(:,:,kg) = [Xinit;atilda'*CC'+repmat(c1',Tub,1)];
    end
end

% re-attributing mean and variance
X_new = X_new.*repmat(Sstd,[Tall,1,max(nDraws,1)])+repmat(Mstd,[Tall,1,max(nDraws,1)]);
