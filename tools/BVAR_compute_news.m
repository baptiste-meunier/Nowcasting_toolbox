function [y_old,y_new,singlenews,actual,forecast,weight,t_miss,v_miss,innov] = BVAR_compute_news(X_old,X_new,Res_in,Par_in,t_fcst,v_news)
% This code inputs two datasets, BVAR parameters, target time index, and
% target variable index. The function then produces Nowcast updates and
% decomposes the changes into news.
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% INPUTS
% - X_old:  Old data matrix (old vintage)
% - X_new:  New data matrix (new vintage)
% - Res_in: Results from current nowcast
% - Par_in: Parameters from current nowcast
% - t_fcst: Index for target time
% - v_news: Index for target variable
%
% OUTPUTS
% - y_old: Old nowcast
% - y_new: New nowcast
% - single_news: News for each data series
% - actual: Observed series release values
% - fore: Forecasted series values
% - weight: news weight
% - t_miss: time index for data releases
% - v_miss: series index for data releases
% - innov: difference between observed and predicted series values ("innovation")
%


%% Initialize parameters
Beta = Res_in.betaHat;
Sigma = Res_in.sigmaHat;
series = Res_in.series;
lags = Par_in.bvar_lags;

%% Initialize data

% Need to restructure data to conform to BVAR_compute_news function
mSeries = 1:Par_in.nM;

month1 = X_new(1:3:end,mSeries);
month2 = X_new(2:3:end,mSeries);
month3 = X_new(3:3:end,:);
X_new = [month1 month2 month3];

month1 = X_old(1:3:end,mSeries);
month2 = X_old(2:3:end,mSeries);
month3 = X_old(3:3:end,:);
X_old = [month1 month2 month3];

v_news = 3*v_news-2; % Adjust target variable index accordingly


%% Initialize variables
[k,n] = size(Beta);
% transition equation
A = zeros(n*lags);
A(1:n,:) = Beta(2:end,:)';% autoregressive coefficients
A(n+1:end,1:end-n) = eye(n*(lags-1));
 c2 = zeros(n*lags,1);
 c2(1:n,:)=Beta(1,:)';% constant (in state eq)

% % useful elements of measurement equation
C = zeros(n,n*lags); C(1:n,1:n) = eye(n);%maps first n states (current month) into observables
c1 = zeros(n,1);
% "shock" impact matrix
Q = 1e-12*eye(n*lags); %are you sure this shouldn't be proper 0s???
Q(1:n,1:n) = Sigma;

% measurement error
R = zeros(n*lags); 
R(1:n,1:n) = diag(ones(1,n)*1e-12);

initX = zeros(lags*n,1);

initX(1:n,:) = nanmean(X_old(1:lags,:),1)';
initX(isnan(initX)) = 0;

initV = BVAR_lyapunov_symm(A, Q, 1, 1e-6);

Res.A = A;
Res.c1 = c1;
Res.c2 = c2;
Res.C =C;
Res.Q = Q;
Res.R = R;
Res.Z_0 = initX;
Res.V_0 = initV;

r = size(Res.C,2);
[~, N] = size(X_new);
singlenews = zeros(1,N);  % Initialize news vector (will store news for each series)

%% NO FORECAST CASE: Already values for variables v_news at time t_fcst

if ~isnan(X_new(t_fcst,v_news))
    Res_old = para_const(X_old, Res, 0);  % Apply Kalman filter for old data
    
    
    for i=1:size(v_news,2)      % Loop for each target variable
        % (Observed value) - (predicted value)
        singlenews(:,v_news(i)) = X_new(t_fcst,v_news(i)) ...
            - Res_old.X_sm(t_fcst,v_news(i));
        
        % Set predicted and observed y values
        y_old(1,i) = Res_old.X_sm(t_fcst,v_news(i));
        y_new(1,i) = X_new(t_fcst,v_news(i));
    end
    
    % Forecast-related output set to empty
    actual=[];forecast=[];weight=[];t_miss=[];v_miss=[];innov=[];
    
    
else
    %% FORECAST CASE (these are broken down into (A) and (B))
    
%     % Initialize series mean/standard deviation respectively
%     Mx = Res.Mx;
%     Wx = Res.Wx;
    % Calculate indicators for missing values (1 if missing, 0 otherwise)
    miss_old=isnan(X_old);
    miss_new=isnan(X_new);
    
    % Indicator for missing--combine above information to single matrix where:
    % (i) -1: Value is in the old data, but missing in new data
    % (ii) 1: Value is in the new data, but missing in old data
    % (iii) 0: Values are missing from/available in both datasets
    i_miss = miss_old - miss_new;
    
    % Time/variable indicies where case (b) is true
    [t_miss,v_miss]=find(i_miss==1);
    
    %% FORECAST SUBCASE (A): NO NEW INFORMATION
    
    if isempty(v_miss)
        % Fill in missing variables using a Kalman filter
        Res_old = para_const(X_old, Res, 0);
        Res_in = para_const(X_new, Res, 0);
        
        % Set predicted and observed y values. New y value is set to old
        y_old = Res_old.X_sm(t_fcst,v_news);
        y_new = y_old;
        % y_new = Res_new.X_sm(t_fcst,v_news);
        
        % No news, so nothing returned for news-related output
        groupnews=[];singlenews=[];gain=[];gainSer=[];
        actual=[];forecast=[];weight=[];t_miss=[];v_miss=[];innov=[];
    else
        %----------------------------------------------------------------------
        %     v_miss=[1:size(X_new,2)]';
        %     t_miss=t_miss(1)*ones(size(X_new,2),1);
        %----------------------------------------------------------------------
        %% FORECAST SUBCASE (B): NEW INFORMATION
        
        % Difference between forecast time and new data time
        lag = t_fcst-t_miss;
        
        % Gives biggest time interval between forecast and new data
        k = max([abs(lag); max(lag)-min(lag)]);
        
        C = Res.C;   % Observation matrix
        R = Res.R';  % Covariance for observation matrix residuals
        
        % Number of new events
        n_news = size(lag,1);
        
        % Smooth old dataset
        Res_old = para_const(X_old, Res, k);
        Plag = Res_old.Plag;
        
        % Smooth new dataset
        Res_in = para_const(X_new, Res, 0);
        
        % Subset for target variable and forecast time
        y_old = Res_old.X_sm(t_fcst,v_news);
        y_new = Res_in.X_sm(t_fcst,v_news);
        
        
        
        P = Res_old.P(:,:,2:end);
        P1=[];  % Initialize projection onto updates
        
        % Cycle through total number of updates
        for i=1:n_news
            h = abs(t_fcst-t_miss(i));
            m = max([t_miss(i) t_fcst]);
            
            % If location of update is later than the forecasting date
            if t_miss(i)>t_fcst
                Pp=Plag{h+1}(:,:,m);  %P(1:r,h*r+1:h*r+r,m)';
            else
                Pp=Plag{h+1}(:,:,m)';  %P(1:r,h*r+1:h*r+r,m);
            end
            P1=[P1 Pp*C(v_miss(i),1:r)'];  % Projection on updates
        end
        
        for i=1:size(t_miss,1)
%             % Standardize predicted and observed values
%             X_new_norm = (X_new(t_miss(i),v_miss(i)) - Mx(v_miss(i)))./Wx(v_miss(i));
%             X_sm_norm = (Res_old.X_sm(t_miss(i),v_miss(i))- Mx(v_miss(i)))./Wx(v_miss(i));
             X_new_norm = X_new(t_miss(i),v_miss(i)) ;
            X_sm_norm = Res_old.X_sm(t_miss(i),v_miss(i));
           
            % Innovation: Gives [observed] data - [predicted data]
            innov(i)= X_new_norm - X_sm_norm;
     
        end
        
        ins=size(innov,2);
        P2=[];
        p2=[];
        
        % Gives non-standardized series weights
        for i=1:size(lag,1)
            for j=1:size(lag,1)
                h=abs(lag(i)-lag(j));
                m=max([t_miss(i),t_miss(j)]);
                
                if t_miss(j)>t_miss(i)
                    Pp=Plag{h+1}(:,:,m); %P(1:r,h*r+1:(h+1)*r,m)';
                else
                    Pp=Plag{h+1}(:,:,m)'; %P(1:r,h*r+1:(h+1)*r,m);
                end
                if v_miss(i)==v_miss(j) & t_miss(i)~=t_miss(j)
                    WW(v_miss(i),v_miss(j))=0;
                else
                    WW(v_miss(i),v_miss(j))=R(v_miss(i),v_miss(j));
                end
                p2=[p2 C(v_miss(i),1:r)*Pp*C(v_miss(j),1:r)'+WW(v_miss(i),v_miss(j))];
            end
            P2=[P2;p2];
            p2=[];
        end
        
        clear temp
        for i=1:size(v_news,2)      % loop on v_news
            % Convert to real units (unstadardized data)
%             totnews(1,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2)*innov';
%             temp(1,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2).*innov;
%             gain(:,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2);
        totnews(1,i) = C(v_news(i),1:r)*P1*inv(P2)*innov';
            temp(1,:,i) = C(v_news(i),1:r)*P1*inv(P2).*innov;
            gain(:,:,i) = C(v_news(i),1:r)*P1*inv(P2);
 
        end
        
        % Initialize output objects
        singlenews = NaN(max(t_miss)-min(t_miss)+1,N); %,size(v_news,2)
        actual     = NaN(N,1);  % Actual forecasted values
        forecast   = NaN(N,1);  % Forecasted values
        weight     = NaN(N,1,size(v_news,2));
        
        % Fill in output values
        for i=1:size(innov,2)
            actual(v_miss(i),1) = X_new(t_miss(i),v_miss(i));
            forecast(v_miss(i),1) = Res_old.X_sm(t_miss(i),v_miss(i));
            
            for j=1:size(v_news,2)
                singlenews(t_miss(i)-min(t_miss)+1,v_miss(i),j) = temp(1,i,j);
                weight(v_miss(i),:,j) = gain(:,i,j);%/Wx(v_miss(i));
            end
        end
        
        singlenews = sum(singlenews,1);      % Returns total news
        
        
        %[v_miss, ~, ~] = unique(v_miss);

        % Need to convert series indices to original data format
        v_miss_convert = repmat(1:length(series)-1, 1, 3);
        v_miss_convert = [v_miss_convert, length(series)]';
        v_miss_convert = v_miss_convert(~isnan(actual));
        v_miss = v_miss_convert;
        [v_miss, v_miss_order] = sort(v_miss); % Sort in ascending order
        
        % Sort news output according to v_miss
        temp_news = {actual, forecast, weight};

        for i = 1:numel(temp_news)
            temp_news{i} = temp_news{i}(~isnan(temp_news{i}));
            temp_news{i} = temp_news{i}(v_miss_order);
        end

        actual = temp_news{1};
        forecast = temp_news{2};
        weight = temp_news{3};

        
    end
    
end

end

function Res = para_const(X, P, lag)
%para_const()    Implements Kalman filter for "News_BVAR.m"
%
%  Syntax:
%    Res = BVAR(X,P,lag)
%
%  Description:
%    para_const() implements the Kalman filter for the news calculation
%    step. This procedure smooths and fills in missing data for a given
%    data matrix X. In contrast to runKF(), this function is used when
%    model parameters are already estimated.
%
%  Input parameters:
%    X: Data matrix.
%    P: Parameters from the model.
%    lag: Number of lags
%
%  Output parameters:
%    Res [struc]: A structure containing the following:
%      Res.Plag: Smoothed factor covariance for transition matrix
%      Res.P:    Smoothed factor covariance matrix
%      Res.X_sm: Smoothed data matrix
%      Res.F:    Smoothed factors
%


% Kalman filter with specified paramaters
% written for
% "MAXIMUM LIKELIHOOD ESTIMATION OF FACTOR MODELS ON DATA SETS WITH
% ARBITRARY PATTERN OF MISSING DATA."
% by Marta Banbura and Michele Modugno

%% Set model parameters and data preparation

% Set model parameters
Z_0 = P.Z_0;
V_0 = P.V_0;
A = P.A;
C = P.C;
Q = P.Q;
R = P.R;
% Mx = P.Mx;
% Wx = P.Wx;
c1=P.c1;
c2 = P.c2;

% Prepare data
[T,~] = size(X);

% % Standardise x
% Y = ((X-repmat(Mx,T,1))./repmat(Wx,T,1))';
Y = X';
%% Apply Kalman filter and smoother
% See runKF() for details about FIS and SKF

Sf = SKF(Y, A, C, Q, R, Z_0, V_0,c2);  % Kalman filter

Ss = FIS(A,Sf);

%% Calculate parameter output

Vs = Ss.VmT(:,:,2:end);  % Smoothed factor covariance for transition matrix
Vf = Sf.VmU(:,:,2:end);  % Filtered factor posterior covariance
Zsmooth = Ss.ZmT;        % Smoothed factors
Vsmooth = Ss.VmT;        % Smoothed covariance values

Plag{1} = Vs;

for jk = 1:lag
    for jt = size(Plag{1},3):-1:lag+1
        As = Vf(:,:,jt-jk)*A'*pinv(A*Vf(:,:,jt-jk)*A'+Q);
        Plag{jk+1}(:,:,jt) = As*Plag{jk}(:,:,jt);
    end
end

% Prepare data for output
Zsmooth=Zsmooth';
x_sm = Zsmooth(2:end,:)*C';  % Factors to series representation
% X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);  % Standardized to unstandardized
X_sm = x_sm; 

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
Res.Plag = Plag;
Res.P = Vsmooth;
Res.X_sm = X_sm;
Res.F = Zsmooth(2:end,:);


end



%______________________________________________________________________
function S = SKF(Y, A, C, Q, R, Z_0, V_0,c2)
%SKF    Applies Kalman filter
%
%  Syntax:
%    S = SKF(Y, A, C, Q, R, Z_0, V_0)
%
%  Description:
%    SKF() applies a Kalman filter. This is a bayesian algorithm. Looping 
%    forward in time, a 'prior' estimate is calculated from the previous 
%    period. Then, using the observed data, a 'posterior' value is obtained.
%
%  Input parameters:
%    y:   Input data.
%    A:   Transition matrix coefficients. 
%    C:   Observation matrix coefficients.
%    Q:   Covariance matrix (factors and idiosyncratic component)
%    R:   Variance-Covariance for observation equation residuals 
%    Z_0: Initial factor values
%    V_0: Initial factor covariance matrix 
%
%  Output parameters:
%    S.Zm:     Prior/predicted factor state vector (Z_t|t-1)  
%    S.ZmU:    Posterior/updated state vector (Z_t|t)  
%    S.Vm:     Prior/predicted covariance of factor state vector (V_t|t-1)  
%    S.VmU:    Posterior/updated covariance of factor state vector (V_t|t)  
%    S.loglik: Value of likelihood function
%    S.k_t:    Kalman gain
%
%  Model:
%   Y_t = C_t Z_t + e_t,     e_t ~ N(0, R)
%   Z_t = A Z_{t-1} + mu_t,  mu_t ~ N(0, Q)
  
%% INSTANTIATE OUTPUT VALUES ---------------------------------------------
  % Output structure & dimensions of state space matrix
  [~, m] = size(C);
  
  % Outputs time for data matrix. "number of observations"
  nobs  = size(Y,2);
  
  % Initialize output
  S.Zm  = NaN(m, nobs);       % Z_t | t-1 (prior)
  S.Vm  = NaN(m, m, nobs);    % V_t | t-1 (prior)
  S.ZmU = NaN(m, nobs+1);     % Z_t | t (posterior/updated)
  S.VmU = NaN(m, m, nobs+1);  % V_t | t (posterior/updated)
  S.loglik = 0;

%% SET INITIAL VALUES ----------------------------------------------------
  Zu = Z_0;  % Z_0|0 (In below loop, Zu gives Z_t | t)
  Vu = V_0;  % V_0|0 (In below loop, Vu guvse V_t | t)
  
  % Store initial values
  S.ZmU(:,1)    = Zu;
  S.VmU(:,:,1)  = Vu;

%% KALMAN FILTER PROCEDURE ----------------------------------------------
  for t = 1:nobs
      %%% CALCULATING PRIOR DISTIBUTION----------------------------------
      
      % Use transition eqn to create prior estimate for factor
      % i.e. Z = Z_t|t-1
      Z   = A * Zu+c2;
      
      % Prior covariance matrix of Z (i.e. V = V_t|t-1)
      %   Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = 
      %   A*Vu*A' + Q
      V   = A * Vu* A' + Q; 
      V   =  0.5 * (V+V');  % Trick to make symmetric
      
      %%% CALCULATING POSTERIOR DISTRIBUTION ----------------------------
       
      % Removes missing series: These are removed from Y, C, and R
      [Y_t, C_t, R_t, ~] = MissData(Y(:,t), C, R); 

      % Check if y_t contains no data. If so, replace Zu and Vu with prior.
      if isempty(Y_t)
          Zu = Z;
          Vu = V;
      else  
          % Steps for variance and population regression coefficients:
          % Var(c_t*Z_t + e_t) = c_t Var(A) c_t' + Var(u) = c_t*V *c_t' + R
          VC  = V * C_t';  
          iF  = inv(C_t * VC + R_t);
          
          % Matrix of population regression coefficients (QuantEcon eqn #4)
          VCF = VC*iF;  

          % Gives difference between actual and predicted observation
          % matrix values
          innov  = Y_t - C_t*Z;
          
          % Update estimate of factor values (posterior)
          Zu  = Z  + VCF * innov;
          
          % Update covariance matrix (posterior) for time t
          Vu  = V  - VCF * VC';
          Vu   =  0.5 * (Vu+Vu'); % Approximation trick to make symmetric
          
          % Update log likelihood 
          S.loglik = S.loglik + 0.5*(log(det(iF))  - innov'*iF*innov);
      end
      
      %%% STORE OUTPUT----------------------------------------------------
      
      % Store covariance and observation values for t-1 (priors)
      S.Zm(:,t)   = Z;
      S.Vm(:,:,t) = V;

      % Store covariance and state values for t (posteriors)
      % i.e. Zu = Z_t|t   & Vu = V_t|t
      S.ZmU(:,t+1)    = Zu;
      S.VmU(:,:,t+1)  = Vu;
  end 
  
  % Store Kalman gain k_t
  if isempty(Y_t)
      S.k_t = zeros(m,m);
  else
      S.k_t = VCF * C_t;
  end
end


  
%______________________________________________________________________
function S = FIS(A, S)
%FIS()    Applies fixed-interval smoother
%
%  Syntax:
%    S = FIS(A, S)
%
%  Description:
%    SKF() applies a fixed-interval smoother, and is used in conjunction 
%    with SKF(). Starting from the final value, this uses a set of 
%    recursions and works backwards. Smoothers are advantageous in that 
%    future values are used for estimation. See  page 154 of 'Forecasting, 
%    structural time series models and the Kalman filter' for more details 
%    (Harvey, 1990).
%
%  Input parameters:
%    A: Transition matrix coefficients. 
%    S: Output from SKF() step. See SKF() for details.
%
%  Output parameters:
%    S: In addition to the output from SKF(), FIS() adds the following
%       smoothed estimates. Note that t = 1...T gives time:
%    - S.ZmT: Smoothed estimate of factor values (Z_t|T) 
%    - S.VmT: Smoothed estimate of factor covariance matrix (V_t|T = Cov(Z_t|T))
%    - S.VmT_1: Smoothed estimate of Lag 1 factor covariance matrix (Cov(Z_tZ_t-1|T))
%
%  Model:
%   Y_t = C_t Z_t + e_t for e_t ~ N(0, R)
%   Z_t = A Z_{t-1} + mu_t for mu_t ~ N(0, Q)

%% ORGANIZE INPUT ---------------------------------------------------------

% Initialize output matrices
  [m, nobs] = size(S.Zm);
  S.ZmT = zeros(m,nobs+1);
  S.VmT = zeros(m,m,nobs+1);
  
  % Fill the final period of ZmT, VmT with SKF() posterior values
  S.ZmT(:,nobs+1)   = squeeze(S.ZmU(:, nobs+1));
  S.VmT(:,:,nobs+1) = squeeze(S.VmU(:,:, nobs+1));

  % Initialize VmT_1 lag 1 covariance matrix for final period
  S.VmT_1(:,:,nobs) = (eye(m)-S.k_t) *A*squeeze(S.VmU(:,:,nobs));
  
  % Used for recursion process. See companion file for details
  J_2 = squeeze(S.VmU(:,:,nobs)) * A' * pinv(squeeze(S.Vm(:,:,nobs)));

  %% RUN SMOOTHING ALGORITHM ----------------------------------------------
  
  % Loop through time reverse-chronologically (starting at final period nobs)
    for t = nobs:-1:1
                
        % Store posterior and prior factor covariance values 
        VmU = squeeze(S.VmU(:,:,t));
        Vm1 = squeeze(S.Vm(:,:,t));
        
        % Store previous period smoothed factor covariance and lag-1 covariance
        V_T = squeeze(S.VmT(:,:,t+1));
        V_T1 = squeeze(S.VmT_1(:,:,t));
      
        J_1 = J_2;
                
        % Update smoothed factor estimate
        S.ZmT(:,t) = S.ZmU(:,t) + J_1 * (S.ZmT(:,t+1) - A * S.ZmU(:,t)) ; 
        
        % Update smoothed factor covariance matrix
        S.VmT(:,:,t) = VmU + J_1 * (V_T - Vm1) * J_1';   
      
        if t>1
            % Update weight
            J_2 = squeeze(S.VmU(:, :, t-1)) * A' * pinv(squeeze(S.Vm(:,:,t-1)));
            
            % Update lag 1 factor covariance matrix 
            S.VmT_1(:,:,t-1) = VmU * J_2'+J_1 * (V_T1 - A * VmU) * J_2';
        end
    end
    
end

    
function [y,C,R,L]  = MissData(y,C,R)
%______________________________________________________________________
% PROC missdata                                                        
% PURPOSE: eliminates the rows in y & matrices Z, G that correspond to     
%          missing data (NaN) in y                                                                                  
% INPUT    y             vector of observations at time t    
%          S             KF system matrices (structure)
%                        must contain Z & G
% OUTPUT   y             vector of observations (reduced)     
%          Z G           KF system matrices     (reduced)     
%          L             To restore standard dimensions     
%                        where # is the nr of available data in y
%______________________________________________________________________
  
  % Returns 1 for nonmissing series
  ix = ~isnan(y);
  
  % Index for columns with nonmissing variables
  e  = eye(size(y,1));
  L  = e(:,ix);

  % Removes missing series
  y  = y(ix);
  
  % Removes missing series from observation matrix
  C  =  C(ix,:);  
  
  % Removes missing series from transition matrix
  R  =  R(ix,ix);

end

