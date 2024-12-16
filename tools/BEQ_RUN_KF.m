function S = BEQ_RUN_KF(y, P, lag)
% This script runs a Kalman filter and smoother
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%

x_0 = P.Z_0;
Sig_0 = P.V_0;
A = P.A;
C = P.C;
Q = P.Q;
R = P.R;
c1 = P.c1;
c2 = P.c2;

% run the filter
S = SKF(y',C,R,A,Q, x_0, Sig_0,c1,c2);

% run the smoother
S = FIS(y',C,R,A,S);

X_sm              = S.AmT'*C';
X_sm(isfinite(y)) = y(isfinite(y));
S.X_sm            = X_sm;

if nargin>2 && lag>0
    
    Ps = S.PmT;
    Pf = S.PmU;
    
    
    Plag{1} = Ps;
    
    for jk = 1:lag
        for jt = size(Plag{1},3):-1:jk+1;
            As = Pf(:,:,jt-jk)*A'*pinv(A*Pf(:,:,jt-jk)*A'+Q);
            Plag{jk+1}(:,:,jt) = As*Plag{jk}(:,:,jt);
        end;
    end;
    
    S.Plag = Plag;
end


%______________________________________________________________________
function S = SKF(Y,Z,R,T,Q,A_0,P_0,c1,c2)
%______________________________________________________________________
% Kalman filter for stationary systems with time-varying system matrices
% and missing data.
%
% The model is        y_t   = Z * a_t + eps_t       
%                     a_t+1 = T * a_t + u_t       
%
%______________________________________________________________________
% INPUT  
%        Y         Data                                 (nobs x n)  
% OUTPUT 
%        S.Am       Predicted state vector  A_t|t-1      (nobs x m)  
%        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)  
%        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)  
%        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)  
%        S.loglik   Value of likelihood function
  
% Output structure & dimensions
  [n, m] = size(Z);
  nobs  = size(Y,2);
  
  S.Am  = nan(m,nobs);   S.Pm  = nan(m,m,nobs);
  S.AmU = nan(m,nobs);   S.PmU = nan(m,m,nobs);
  S.ZF = cell(nobs);
  S.V = cell(nobs);
  %______________________________________________________________________
  Au = A_0;  % A_0|0;
  Pu = P_0;  % P_0|0
  S.loglik = 0;

  for t = 1:nobs
%       t
      % A = A_t|t-1   & P = P_t|t-1

      A   = T*Au+c2;
      P   = T*Pu*T' + Q;
      P   =  0.5 * (P+P');
      
      % handling the missing data
      [y_t,Z_t,R_t,c1_t] = MissData(Y(:,t),Z,R,c1);

      if isempty(y_t)
          Au = A;
          Pu = P;
          ZF = zeros(m,0);
          V = zeros(0,1);

      else
          PZ = P*Z_t';
          iF  = eye(length(c1_t))/(Z_t*PZ + R_t);
          ZF = Z_t'*iF;
          PZF = P*ZF;

          V   = y_t - Z_t*A-c1_t;
          Au  = A  + PZF * V;
          Pu  = P  - PZF * PZ';
          Pu   =  0.5 * (Pu+Pu');
          S.loglik = S.loglik + 0.5*(log(det(iF))  - V'*iF*V);
      end
      S.ZF{t} = ZF;
      S.Am(:,t)   = A;
      S.Pm(:,:,t) = P;
      S.V{t} = V;

      S.AmU(:,t)    = Au;
      S.PmU(:,:,t)  = Pu;
  end % t
 

%______________________________________________________________________
function S = FIS(Y,Z,R,T,S)
%______________________________________________________________________
% Fixed intervall smoother (see Harvey, 1989, p. 154)
% FIS returns the smoothed state vector AmT and its covar matrix PmT             
% Use this in conjnuction with function SKF
%______________________________________________________________________
% INPUT  
%        Y         Data                                 (nobs x n)  
%        S Estimates from Kalman filter SKF                                                          
%          S.Am   : Estimates     a_t|t-1                  (nobs x m) 
%          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
%          S.AmU  : Estimates     a_t|t                    (nobs x m) 
%          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
% OUTPUT 
%        S Smoothed estimates added to above
%          S.AmT  : Estimates     a_t|T                    (nobs x m) 
%          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
%        where m is the dim of state vector and t = 1 ...T is time

  [m, nobs] = size(S.Am);
  S.AmT     = zeros(m,nobs);
  S.PmT     = zeros(m,m,nobs);

  r = zeros(m,1);
  N = zeros(m,m);

  for t = nobs:-1:1
	[~,Z_t] = MissData(Y(:,t),Z,R,zeros(length(Y(:,t)),1));

    L = T*(eye(m)-squeeze(S.Pm(:,:,t))*S.ZF{t}*Z_t);
	r = S.ZF{t}*S.V{t}+L'*r;
    N = S.ZF{t}*Z_t+L'*N*L;
    
	S.AmT(:,t)   = S.Am(:,t)+ squeeze(S.Pm(:,:,t))*r;
    S.PmT(:,:,t) = S.Pm(:,:,t)- squeeze(S.Pm(:,:,t))*N*squeeze(S.Pm(:,:,t)); 
  end

%______________________________________________________________________
function [y,C,R,c1]  = MissData(y,C,R,c1)
%______________________________________________________________________
% PROC missdata                                                        
% PURPOSE: eliminates the rows in y & matrices Z, G that correspond to     
%          missing data (NaN) in y                                                                                  
% INPUT    y             vector of observations at time t  (n x 1 )    
%          S             KF system matrices             (structure)
%                        must contain Z & G
% OUTPUT   y             vector of observations (reduced)   (# x 1)     
%          Z G           KF system matrices     (reduced)   (# x ?)     
%          L             To restore standard dimensions     (n x #)     
%                        where # is the nr of available data in y
%______________________________________________________________________
  ix = ~isnan(y);

  y  =  y(ix);
  c1 =  c1(ix);
  C  =  C(ix,:);  
  R  =  R(ix,ix);
