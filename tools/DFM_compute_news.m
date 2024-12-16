function [now_old,now_new,actual,forecast,weight,v_miss] = DFM_compute_news(X_new,X_old,Res,targetmonth,gdp_index)
% This function computes the impact of news on the nowcast. 
% 
% Code from Bańbura, M., and Modugno, M. (2014). “Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%
% INPUTS 
% - X_new [matrix] = dataset used for current iteration of DFM (with NaN)
% - X_old [matrix] = dataset used for the previous iteration of DFM (with NaN)
% - Res [structure] = results
%       o X_sm [matrix] = Kalman-smoothed dataset (NaN replaced with expectations)
%       o F [matrix] = smoothed states 
%                      Rows give time
%                      Columns are organized according to Res.C
%       o C [matrix] = observation matrix
%                      Rows correspond to each series
%                      Columns give factor loadings
%                      For example, col. 1-5 give loadings for the first block and are organized in reverse-chronological order (f^G_t, f^G_t-1, f^G_t-2, f^G_t-3, f^G_t-4)
%       o R [matrix] = covariance for observation matrix residuals
%       o A [matrix] = transition matrix
%                      Square matrix that follows the same organization scheme as Res.C's columns
%                      Identity matrices are used to account for matching terms on the left and right-hand side.
%                      For example, we place an I4 matrix to account for matching (f_t-1; f_t-2; f_t-3; f_t-4) terms.
%       o Q [matrix] = covariance for transition equation residuals
%       o Mx [vector] = series mean
%       o Wx [vector] = series standard deviation
%       o Z_0 [vector] = initial value of state
%       o V_0 [matrix] = initial value of covariance matrix
%       o r [scalar] = number of estimated factors
%       o p [scalar] = number of lags in transition equation
%       o L [scalar] = log-likelihood
%       o groups [vector] = indicator for group belonging
%       o series [cell vector] = mnemotechnics of series in X_sm (generally their Haver handle)
%       o name_descriptor [cell vector] = full names of series in X_sm
% - targetmonth [scalar] = observation of the targeted month
% - gdp_index [scalar] = column with the target variable
%
% OUTPUTS
% - now_old [scalar] = old nowcast (for targetmonth)
% - now_new [scalar] = old nowcast (for targetmonth)
% - actual [vector] = values from the new data release
% - forecast [vector] = predicted values from the Kalman smoother
% - weight [vector] = weights of new data release
% - v_miss [vector] = indexes variables that miss information relative to new data release
%


%% Case 1: Code was run although GDP data is already released

if ~isnan(X_new(targetmonth,gdp_index))
    
    % The assumption here is that we ran the model just before the release
    % of GDP, so most data was already available, but now also GDP was
    % released.
    Res_old = para_const(X_old,Res,0); % Apply Kalman filter for old data
    now_old = Res_old.X_sm(targetmonth,gdp_index);
    now_new = X_new(targetmonth, gdp_index);
    
    actual = []; forecast = []; weight = []; v_miss = [];   % output of function
else
    %% Check whether there is new data available
    
    % Initialize series mean/standard deviation respectively
    Mx = Res.Mx;
    Wx = Res.Wx;
    
    
    % create an indicator whether for a series there are news
    miss_old=isnan(X_old);
    miss_new=isnan(X_new);
    i_miss = miss_old - miss_new;
    
    % i_miss:
    % case (a) value -1: Data point present in old dataset, but not in new
    % case (b) value 1:  Data point present in new dataset, but missing in
    % old
    % case (c) value 0: Data missing/available in both datasets
    
    [t_miss,v_miss]=find(i_miss==1);        % v_miss: indexes variables that miss information relative to new data release
                                            % t_miss: time index where (b) is true
    %% Case 2: There's no new data and therefore no news (but this step is important for data-revision)
    
    if isempty(v_miss)
       
        % Fill in missing values in the datasets
        Res_old = para_const(X_old,Res,0);      % lags=0 : Then there's no change to runKF
        
        now_old = Res_old.X_sm(targetmonth,gdp_index);
        now_new = now_old;
        actual = []; forecast = []; weight = [];
        
    else
    %% Case 3: There's new data: Start news computation
    
    % Input "targetmonth" gives the row (/date) for which we want to nowcast gdp.
    lag = targetmonth - t_miss;      % lag gives the difference between the month of the nowcast and the month in which new data was released for every series with new information
    % max time interval between nowcast and new data release. Second term
    % takes care of the case that e.g. lag = [1,-1].
    k = max([abs(lag); max(lag)-min(lag)]);
    % Number of new datapoints
    n_news = size(lag,1);
    
    % Initialize matrices : 
    C = Res.C;          % Observation matrix
    R = Res.R;          % VCV of observation equation residuals
    % Smooth new dataset
    Res_new = para_const(X_new,Res,0);
    
    % Smooth old data
    Res_old = para_const(X_old,Res,k);      % Here X_old includes the revisions (input)
    Plag = Res_old.Plag;                    % We use the Plag matrices from Res_old: we want to compare the forecast to the actual values. 
                                            % The forecast involves the Plag output in the computation
    % The nowcasts are:
    now_old = Res_old.X_sm(targetmonth,gdp_index);
    now_new = Res_new.X_sm(targetmonth, gdp_index);
        
    % Compute E[y_k,tk * I_{v+1}] (without premultiplied Lambda_k, multiplied later below) 
    P1 = [];                
    
    for i=1:n_news
        
        
        h = abs(targetmonth-t_miss(i));     
        m = max([t_miss(i) targetmonth]);  
        
        if t_miss(i) > targetmonth
            Pp = Plag{h+1}(:,:,m);          
        else
            Pp = Plag{h+1}(:,:,m)';         
        end
        % P1 has dimension (m x n_news) (recall: m = #of elements in state vector)
        P1 = [P1 Pp*C(v_miss(i),:)'];       % P1 in BM2010 notation: E[ (f_{tk} - E[f_{tk}|Omega_v]) (f_{tj} - E[f_{tj}|Omega_v]) ] * (Lambda_{ij})'
    
    end
    
    for i = 1:size(t_miss,1)
        % standardize predicted and observed values
        X_new_norm = (X_new(t_miss(i),v_miss(i)) - Mx(v_miss(i)))./Wx(v_miss(i));
        X_sm_norm = (Res_old.X_sm(t_miss(i),v_miss(i))- Mx(v_miss(i)))./Wx(v_miss(i));
        innovation(i) = X_new_norm - X_sm_norm;         
    end
        
    % Compute E(I_{v+1} * I_{v+1}')
    P2=[];
    p2=[];
    
     % Gives non-standardized series weights
    for i=1:size(lag,1)
        for j=1:size(lag,1)
            h=abs(lag(i)-lag(j));
            m=max([t_miss(i),t_miss(j)]);
            
            if t_miss(j)>t_miss(i)
                Pp=Plag{h+1}(:,:,m); 
            else
                Pp=Plag{h+1}(:,:,m)'; 
            end
            if v_miss(i)==v_miss(j) & t_miss(i)~=t_miss(j)
                WW(v_miss(i),v_miss(j))=0;
            else
                WW(v_miss(i),v_miss(j))=R(v_miss(i),v_miss(j));
            end
            % BM2010 notation: E(I_{v+1,j} * I_{v+1,l}) = Lambda_{ij} * E[ (f_{tj}-E[f_tj|Omega_{v}]) (f_{tl}-E[f_{tl}|Omega_{v}]) ] * (Lambda_{il}') + 1_{j=l} * R_jj           
            p2=[p2 C(v_miss(i),:)*Pp*C(v_miss(j),:)'+WW(v_miss(i),v_miss(j))];
            
        end
        % P2 has dimension (n_news x n_news) ( size(lag,1) = n_news = size(t_miss,1))
        P2=[P2;p2];
        p2=[];
    end

    % Compute total news (equation 16 and 26 in BM2010, scaled by Wx). This
    % is the Revision: E[y_{k,tk}|I_{v+1}]
    %totalnews = Wx(gdp_index)*C(gdp_index,:)*P1*inv(P2)*innovation';     % (scalar) This is the change in the nowcast of GDP that is due to the news
    gain = Wx(gdp_index)*C(gdp_index,:)*P1*inv(P2);         % Size: 1 x n_news
    
    %% Output:
    % Actual: simply values from the new data release
    actual = NaN(n_news,1);
    % Forecast: predicted values from the Kalman smoother
    forecast = NaN(n_news,1);
    % weight: for each variable (N x 1)
    weight = NaN(n_news,1);

    for i = 1:size(innovation,2)        % innovation is a row vector
        actual(i,1) = X_new(t_miss(i),v_miss(i));                   
        forecast(i,1) = Res_old.X_sm(t_miss(i),v_miss(i));          
        weight(i,1) = gain(i) / Wx(v_miss(i));                      % B_{v+1} vector in BM2010
    end
    
    end     % if no new data was released
    
end         % If there's GDP already    


end

%==========================================================================
%% KALMAN FILTER & SMOOTHER PROCEDURES
%==========================================================================


function Res = para_const(X, P, lag)
%para_const()    Implements Kalman filter for "News_DFM.m"
%
%  Syntax:
%    Res = para_const(X,P,lag)
%
%  Description:
%    para_const() implements the Kalman filter for the news calculation
%    step. This procedure smooths and fills in missing data for a given 
%    data matrix X. In contrast to runKF(), this function is used when
%    model parameters are already estimated.
%
%  Input parameters:
%    X: Data matrix. 
%    P: Parameters from the dynamic factor model.
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
Mx = P.Mx;
Wx = P.Wx;

% Prepare data
[T,~] = size(X);

% Standardise x
Y = ((X-repmat(Mx,T,1))./repmat(Wx,T,1))';

%% Apply Kalman filter and smoother
% See runKF() for details about FIS and SKF

Sf = SKF(Y, A, C, Q, R, Z_0, V_0);  % Kalman filter

Ss = FIS(A, Sf);  % Smoothing step

%% Calculate parameter output
Vs = Ss.VmT(:,:,2:end);  % Smoothed factor covariance for transition matrix
Vf = Sf.VmU(:,:,2:end);  % Filtered factor posterior covariance / updated covariance of factor state vector
Zsmooth = Ss.ZmT;        % Smoothed factors
Vsmooth = Ss.VmT;        % Smoothed covariance values

Plag{1} = Vs;

for jk = 1:lag
    for jt = size(Plag{1},3):-1:lag+1                   % for T(length of time series) until lag+1
        As = Vf(:,:,jt-jk)*A'*pinv(A*Vf(:,:,jt-jk)*A'+Q);   % eq. 3.6.16c in Harvey
        Plag{jk+1}(:,:,jt) = As*Plag{jk}(:,:,jt);
    end
end

% Prepare data for output
Zsmooth=Zsmooth';
x_sm = Zsmooth(2:end,:)*C';  % Factors to series representation
X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);  % Standardized to unstandardized

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
Res.Plag = Plag;                % unsure what this is yet
Res.P = Vsmooth;                % Same output and notation as in runKF
Res.X_sm = X_sm;                % same output and notation: unstandardized smoothed data series
Res.F = Zsmooth(2:end,:);       % same output and notation as in runKF in estimation file

end


%______________________________________________________________________
function S = SKF(Y, A, C, Q, R, Z_0, V_0)
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
      Z   = A * Zu;
      
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
          iF  = inv(C_t * VC + R_t);        % Inverse of equation 3.2.3c in Harvey1989
          
          % Matrix of population regression coefficients (QuantEcon eqn #4)
          VCF = VC*iF;  

          % Gives difference between actual and predicted observation
          % matrix values
          innov  = Y_t - C_t*Z;
          
          % Update estimate of factor values (posterior)
          Zu  = Z  + VCF * innov;           % equation 3.2.3a in Harvey1989
          
          % Update covariance matrix (posterior) for time t
          Vu  = V  - VCF * VC';             % equation 3.2.3b in Harvey1989
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
      S.k_t = VCF * C_t;        % Kalman gain in context of updating equation
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
%  Model: (CHANGE IN NOTATION HERE AGAIN)
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
        S.ZmT(:,t) = S.ZmU(:,t) + J_1 * (S.ZmT(:,t+1) - A * S.ZmU(:,t)) ;       % equation 3.6.16a Harvey1989
        
        % Update smoothed factor covariance matrix
        S.VmT(:,:,t) = VmU + J_1 * (V_T - Vm1) * J_1';                          % equation 3.6.16b 
      
        if t>1
            % Update weight
            J_2 = squeeze(S.VmU(:, :, t-1)) * A' * pinv(squeeze(S.Vm(:,:,t-1)));    % equation 3.6.16c
            
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
  
  % remove missing in VCV
  R  =  R(ix,ix);

end


