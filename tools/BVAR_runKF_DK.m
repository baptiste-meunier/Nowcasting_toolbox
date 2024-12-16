function [xFS, logLik] = BVAR_runKF_DK(y, A, C, Q, R, x_0, Sig_0,c1,c2)
% This function runs Kalman filter and smoother
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

% run the filter
S = SKF(y,C,R,A,Q, x_0, Sig_0,c1,c2);
% run the smoother
S = FIS(y,C,R,A,S);

xFS = S.AmT;
logLik.logLik = S.logLik;
logLik.T1 = S.T1;
logLik.T2 = S.T2;
logLik.ZF = S.ZF;

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
%        S.logLik   Value of likelihood function

% Output structure & dimensions
[n, m] = size(Z);
nobs  = size(Y,2);

Am  = nan(m,nobs);   Pm  = nan(m,m,nobs);
%AmU = nan(m,nobs);   %PmU = nan(m,m,nobs);
ZF = cell(nobs,1);
V = cell(nobs,1);
logLik = 0;

%______________________________________________________________________
Au = A_0;  % A_0|0;
Pu = P_0;  % P_0|0

T1 = zeros(nobs,1);
T2 = T1;
Tt = T';


for t = 1:nobs
    %       t
    % A = A_t|t-1   & P = P_t|t-1
    
    A   = T*Au+c2;
    P   = T*Pu*Tt + Q;
    P   =  (P+P')/2;
    
    % handling the missing data
    [y_t,Z_t,R_t,c1_t] = MissData(Y(:,t),Z,R,c1);
    
    if isempty(y_t)
        Au = A;
        Pu = P;
        ZF_t = zeros(m,0);
        V_t = zeros(0,1);
        
    else
        Z_tt = Z_t';
        PZ = P*Z_tt;
        F  = (Z_t*PZ + R_t);
        Finv = BVAR_inv2(F);
        ZF_t = Z_tt*Finv;
        PZF = P*ZF_t; %Kalman gain
        V_t   = y_t - Z_t*A-c1_t;
        Au  = A  + PZF * V_t;
        Pu  = P  - PZF * PZ';
        Pu   =  (Pu+Pu')/2;
        T1(t) = .5*log(det(Finv));
        T2(t) = -.5*V_t'*Finv*V_t;
        logLik = logLik + T1(t)+T2(t)-.5*n*log(2*pi);

    end
    ZF{t} = ZF_t;
    Am(:,t)   = A;
    Pm(:,:,t) = P;
    V{t} = V_t;
    
%    AmU(:,t)    = Au;
%    PmU(:,:,t)  = Pu;
    

    
end % t

S.Am = Am;
S.AmU = Au;
S.ZF = ZF;
S.V = V;
S.Pm = Pm;
%S.PmU = PmU;
S.logLik = logLik;
S.T1 = T1;
S.T2 = T2;
end

%______________________________________________________________________
function SS = FIS(Y,Z,R,T,S)
%______________________________________________________________________
% Fixed interval smoother (see Durbin and Koopman, 2001, p. 64-71)
% FIS returns the smoothed state vector AmT and its covar matrix PmT
% Use this in conjnuction with function SKF
%______________________________________________________________________
% INPUT
%        Y         Data                                 (nobs x n)
%        S Estimates from Kalman filter SKF
%          S.Am   : Estimates     a_t|t-1                  (nobs x m)
%          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
% OUTPUT
%        S Smoothed estimates added to above
%          S.AmT  : Estimates     a_t|T                    (nobs x m)
%        where m is the dim of state vector and t = 1 ...T is time


SS = S;
[m, nobs]        = size(S.Am);
AmT           = zeros(m,nobs);
% PmT           = zeros(m,m,nobs);
AmT(:,nobs)   = S.AmU  ;

r = zeros(m,1);
ZF = S.ZF;
V = S.V;
Am = S.Am;
Pm = S.Pm;

ee = eye(m);

for t = nobs:-1:1
    
    [~,Z_t] = MissData(Y(:,t),Z,R,zeros(length(Y(:,t)),1));
    
    Pmt = squeeze(Pm(:,:,t));
    ZFt = ZF{t};
    r = ZFt*V{t}+(T*(ee-Pmt*ZFt*Z_t))'*r;
    AmT(:,t) = Am(:,t)+ Pmt*r;
    
end

SS.AmT = AmT;

end

%______________________________________________________________________
function [y,C,R,c1]  = MissData(y,C,R,c1)
%______________________________________________________________________
% PROC missdata
% PURPOSE: eliminates the rows in y & matrices C, R and vector c1 that correspond to
%          missing data (NaN) in y
%______________________________________________________________________
ix = ~isnan(y);

y  =  y(ix);
c1 =  c1(ix);
C  =  C(ix,:);
R  =  R(ix,ix);

end
