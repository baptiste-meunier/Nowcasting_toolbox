function Res = DFM_estimate(X,Par,Res_old)
% This function estimates the results of the DFM
%
% Code from Bańbura, M., and Modugno, M. (2014). "Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%
% Block identification from Delle Chiaie, S., Ferrara, L., and Giannone, D. 
% (2022). “Common factors of commodity prices”, Journal of Applied 
% Econometrics, 37(3), 461–476
%
% INPUTS 
% - X [matrix] = input dataset
% - Par [structure] = parameters of models
%       o yearq [scalar] = year for which latest GDP is available
%       o qend [scalar] = quarter for which latest GDP is available
%       o nM [scalar] = number of monthly variables
%       o nQ [scalar] = number of quarterly variables
%       o blocks [vector/matrix] = block identification (DFM)
%       o r [scalar] = number of estimated factors (DFM)
%       o block_factors [scalar] = switch on whether to do block factors (=1) or not (=0) (DFM)
%       o p [scalar] = number of lags (DFM)
%       o idio [scalar] = idiosyncratic component specification: 1 = AR(1), 0 = iid (DFM)
%       o thresh [scalar] = threshold for convergence of EM algorithm (DFM)
%       o max_iter [scalar] = number of iterations for the EM algorithm (DFM)
%       o lagM [scalar] = number of lags for monthly regressor(s) (in quarterly terms) (BEQ)
%       o lagQ [scalar] = number of lags for quarterly regressor(s) (in quarterly terms) (BEQ) 
%       o lagY [scalar] = number of lags for the endogenous variable (in quarterly terms) (BEQ)
%       o Dum [matrix] = dates of the dummies (year, month). Format should be a k x 2 matrix with k = number of dummies (each dummy on a row) and year / month in columns (BEQ)
%       o bvar_lags [scalar] = number of lags (BVAR)
%       o blocks_full [vector/matrix] = block identification (for do_loop=2)
%       o size_data_in [scalar] = lengths of raw data set
% - Res_old [structure] = results from a previous iteration
%
% OUTPUT
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
%

nQ = Par.nQ;

nM = size(X,2)-nQ;

thresh = Par.thresh;
r = Par.r;
if max(r) > nM
    warning("Number of factors exceeding number of monthly variables. Setting number of factors to smaller value")
    r(r>nM) = nM;
end
p = Par.p;
max_iter = Par.max_iter;

% - Idiosyncratic component: AR(1) choice via Par.idio
if Par.idio == 1
    i_idio = logical([ones(nM,1);zeros(nQ,1)]); 
else
    i_idio = logical([zeros(nM,1);zeros(nQ,1)]); 
end


R_mat = [2 -1 0 0 0;...
    3 0 -1 0 0;...
    2 0 0 -1 0;...
    1 0 0 0 -1]; % matrix of constraints on the loading of quarterly data

q = zeros(4,1);

if p > 5
    p = 5;
    warning('This routine cannot handle more than 5 lags; p set to 5.')
end


blocks = Par.blocks;

%--------------------------------------------------------------------------
% Preparation of the data
%--------------------------------------------------------------------------
[T,N] = size(X);

if nargin<3
    % Standardise x
    Mx = nanmean(X);
    Wx = (nanstd(X));
    
else
    Mx = Res_old.Mx;
    Wx = Res_old.Wx;
end

xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);

%--------------------------------------------------------------------------
% Initial Conditions
%--------------------------------------------------------------------------

%Removing missing values (for initial estimators)
optNaN.method = 2; % Remove leading and closing zeros
optNaN.k = 3;

if nargin<3
    [A, C, Q, R, Z_0, V_0] = DFM_InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
else
    A = Res_old.A;
    C = Res_old.C;
    Q = Res_old.Q;
    R = Res_old.R;
    Z_0 = Res_old.Z_0;
    V_0 = Res_old.V_0;
end
    
    
    
% some auxiliary variables for the iterations
previous_loglik = -inf;
num_iter = 0;
LL = -inf;
converged = 0;

% y for the estimation is WITH missing data
y = xNaN';


%--------------------------------------------------------------------------
%THE EM LOOP
%--------------------------------------------------------------------------

%The model can be written as
%y = C*Z + e;
%Z = A*Z(-1) + v
%where y is NxT, Z is (pr)xT, etc

%remove the leading and ending nans for the estimation
optNaN.method = 3;
y_est = DFM_remNaNs_spline(xNaN,optNaN)';

% disp('Change in Loglik:')

while (num_iter < max_iter) & ~converged
    [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = ...
        DFM_EMstep(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks);
    
    C = C_new;
    R = R_new;
    A = A_new;
    Q = Q_new;

    % Checking convergence
    if num_iter>2
    [converged,decrease(num_iter+1)] = DFM_em_converged(loglik, previous_loglik, thresh,1);
    end
    
    LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
    
    if mod(num_iter,10)==0
        disp(['Running the ',num2str(num_iter),'th iteration'])
        % disp(diff(LL(end-1:end))/mean(LL(end-1:end))*200)
    end
end

%----------------------------------------------
%final run of the Kalman filter
%----------------------------------------------
Zsmooth = DFM_runKF(y, A, C, Q, R, Z_0, V_0)';
x_sm = Zsmooth(2:end,:)*C';


Res.X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);
Res.F = Zsmooth(2:end,:);

%--------------------------------------------------------------------------
%   Loading the structure with the results
%--------------------------------------------------------------------------
Res.C = C;
Res.R = R;
Res.A = A;
Res.Q = Q;
Res.Mx = Mx;
Res.Wx = Wx;
Res.Z_0 = Z_0;
Res.V_0 = V_0;
Res.r = r;
Res.p = p;
Res.L = loglik;

end