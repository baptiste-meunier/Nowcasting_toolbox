function Res = BVAR_estimate(yest_out,Par,datet_in)
% This script runs a mixed-frequency BVAR following the block BVAR of 
% Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. (2022). 
% "Nowcasting with large Bayesian vector autoregressions", Journal of
% Econometrics, 231, 500-519
%
% Code adapted from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & 
% Sokol, A. (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
%
% INPUTS
% - yest_out [matrix] = input data (monthly and quarterly, quarterly with NaN on 1st and 2nd months)
% - Par [structure] = parameters of the model
%       o nM [scalar] = number of monthly variables to be used (always the first ones)
%       o nQ [scalar] = number of quarterly variables to be used (always the last ones)
%       o bvar_lags [scalar] = number of lags in the BVAR
% - datet_in [matrix] = dates (year / month)
%
% OUTPUTS
% - Res [structure] = results of the model
%       o X_sm [matrix] = dataset with interpolated data for target variable
%

%--------------------------------------------------------------------------
%% Step 0. Set parameters and run preliminary checks
%--------------------------------------------------------------------------

warning off

% Set starting values for optimisation (same as Cimadomo et al., 2022)
lambda0 = 0.2;
miu0 = 1;
theta0 = 1;
alpha0 = 2;

% Ensure data extends from month 1 of a quarter until month 3 of another (needed in bbvar code)
[X_ext,datet_ext] = BVAR_filldata(yest_out,datet_in);

% Check that the number of observations is sufficient
% Otherwise the code risks getting stuck in 'bad gradient' error
% For monthly data
min_obs = 3*2*Par.bvar_lags + 1; % Same criteria as in BEQ_estimate (for consistency)
iMax = find(~any(isnan(X_ext),2),1,'last'); % position of the last non-NaN
idx_keep_Xm = [];
for cc = 1:Par.nM        
    iX = find(~isnan(X_ext(:,cc)),1,'first');
    if iMax - iX + 1 > min_obs
        idx_keep_Xm = [idx_keep_Xm, cc];
    end
end

% Same for quarterly data
min_obs = 2*Par.bvar_lags + 1; % Same criteria as in BEQ_estimate (for consistency)
idx_keep_Xq = [];
if Par.nQ > 1
    for cc = (Par.nM+1):(Par.nM+Par.nQ-1)        
        iX = find(~isnan(X_ext(:,cc)),1,'first');
        if iMax - iX + 1 > min_obs
            idx_keep_Xq = [idx_keep_Xq, cc];
        end
    end
end

X_ext = X_ext(:,[idx_keep_Xm,idx_keep_Xq,end]);
Par.nM = length(idx_keep_Xm); % changing nM
Par.nQ = length(idx_keep_Xq) + 1; % changing nQ still +1 for the target series

% Add return if all is NaN
if Par.nM == 0
    disp('Not enough variables satisfying criteria, returning the initial dataset with NaN')
    Res.X_sm = yest_out;
    return
end

% Set parameters for B-BVAR
stationary = 1:size(X_ext,2);
mSeries = 1:Par.nM;


%--------------------------------------------------------------------------
%% Step 1. Run the estimations
%--------------------------------------------------------------------------

thresh = Par.bvar_thresh; % threshold for convergence
max_iter = Par.bvar_max_iter; % number of iterations


try
    Res = BVAR_bbvar(X_ext,Par.bvar_lags, ...
                     mSeries,stationary, ...
                     lambda0,theta0,miu0,alpha0,thresh,max_iter);
catch
    disp('BVAR not converging, returning the initial dataset with NaN')
    Res.X_sm = yest_out;
    return
end

%--------------------------------------------------------------------------
%% Step 2. Convert forecasts back into X_sm
%--------------------------------------------------------------------------

% Store the previous X_sm
Res.X_sm_init = Res.X_sm;
Res.X_sm = [];

% Decompose by months (see structure of input data in bbvar code)
month1 = Res.X_sm_init(:,1:Par.nM);
month1 = [month1,nan(size(month1,1),Par.nQ)];
month2 = Res.X_sm_init(:,(Par.nM+1):(2*Par.nM));
month2 = [month2,nan(size(month2,1),Par.nQ)];
month3 = Res.X_sm_init(:,(2*Par.nM+1):end);
nobs = size(month3,1);

% Recreate X_sm in monthly frequency
n_init = 3*Par.bvar_lags;
Res.X_sm = NaN(n_init,Par.nM+Par.nQ);
for ll = 1:nobs
    n_line = n_init + 3*(ll-1);
    Res.X_sm(n_line+1,:) = month1(ll,:);
    Res.X_sm(n_line+2,:) = month2(ll,:);
    Res.X_sm(n_line+3,:) = month3(ll,:);
end

% Check dimensions match
if size(Res.X_sm) ~= size(X_ext)
    error('In BVAR, sizes of re-created output data after BVAR do not match')
end

% Cut back to original dimension (of xest)
idx1 = find(datet_ext(:,1)==datet_in(1,1) & datet_ext(:,2)==datet_in(1,2));
idx2 = find(datet_ext(:,1)==datet_in(end,1) & datet_ext(:,2)==datet_in(end,2));
Res.X_sm = Res.X_sm(idx1:idx2,:);

% Check dimensions match
if size(Res.X_sm) ~= size(yest_out)
    error('In BVAR, sizes of output and input do not match')
end

end