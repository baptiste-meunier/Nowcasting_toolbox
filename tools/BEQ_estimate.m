function Res = BEQ_estimate(xest_out,Par,datet_in,nameseries,re_estimate,coeffs_in)
% This script runs a combination of bridge equations in the spirit of
% Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B. (2023). "Nowcasting
% employment in the euro area", Working Paper Series, No 2815, European
% Central Bank
%
% Code adapted from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B. 
% (2023). "Nowcasting employment in the euro area", Working Paper Series, 
% No 2815, European Central Bank
%
% INPUTS
% - xest [matrix] = input data (monthly and quarterly, quarterly with NaN on 1st and 2nd months)
% - Par [structure] = parameters of the model
%       o nM [scalar] = number of monthly variables to be used (always the first ones)
%       o nQ [scalar] = number of quarterly variables to be used (always the last ones)
%                       with target variable (generally GDP) at the very end
%       o lagM [scalar] = number of lags for monthly regressor(s) (in quarterly terms) (BEQ)
%       o lagQ [scalar] = number of lags for quarterly regressor(s) (in quarterly terms) (BEQ) 
%       o lagY [scalar] = number of lags for the endogenous variable (in quarterly terms) (BEQ)
%       o type [scalar] = type of interpolation performed (see BEQ_estimate for details)
%       o Dum [matrix] = dates of the dummies (year, month). Format should be a k x 2 matrix with k = number of dummies (each dummy on a row) and year / month in columns (BEQ)
% - datet [matrix] = dates (year / month)
% - nameseries [cell array] = name of the series (short version)
% - type [scalar] = switch on how to interpolate data
%       o 901 = BVAR with all data
%       o 902 = BVAR on selected data
%       o 903 = univariate BVAR
%       o 904 = all three above
% - re_estimate [sclar] = whether to re-estimate OLS or use the coefficients passed to the functions
% - coeffs_in [matrix] = matrix of coefficients (used if the OLS is not re-estimated)
%
% OUTPUTS
% - Res [structure] = results of the model
%       o Y_fcst [vector] = predictions of target variable (averaged over all bridge equations)
%       o Y_fcst_indiv [matrix] = predictions for each individual model
%       o beq_combinations [matrix] = specifications of all bridge equations
%           > col. 1 = type of interpolation
%               + 901 = BVAR with all data
%               + 902 = BVAR on selected data
%               + 903 = univariate BVAR
%           > cols. 2 & 3 = monthly regressors
%           > col. 4 = quarterly regressor (if any)
%       o Date_fcst [matrix] = dates (year / month)
%       o X_sm [matrix] = dataset with interpolated data for target variable
%       o contrib_X_sm [matrix] = contributions of individual variables to the forecast
%           > NB: does not sum to the actual data (given contributions from the error term)
%                 but sums for forecasts
%       o contrib [3D-matrix] = contributions of individual variables to the individual forecasts
%           > NB: does not sum to the actual data (given contributions from the error term)
%                 but sums for forecasts
%           > col. 1 = time
%           > col. 2 = variables (including constant and dummies on the right)
%           > col. 3 = bridge equation
%       o contrib_names [cell array] = names of the variables that correspond to contributions
%       o coeffs [matrix] = coefficients for the bridge equations
%           > NB: match the beq_combinations
%

%--------------------------------------------------------------------------
%% Step 0. Set parameters and run preliminary checks
%--------------------------------------------------------------------------

warning off

% Complement BEQ parameters
Par.nDraw = 0; % number of draws (0 = point forecast)

% BVAR parameters (for interpolation)
% We take the same parameters as in Banbura et al. (2023)
Par_BVAR.nDraw = 0; % number of draws
Par_BVAR.KK = 0;
Par_BVAR.iRW = 0;
Par_BVAR.lags = 6; % number of lags

% Check dimensions
if (size(datet_in,1) ~= size(xest_out,1))
    error('In bridge equations, number of observations for time (datet) and data (xest) do not match')
end

% Check type
if ~ismember(Par.type,901:904)
    error('Type for bridge equations (in Par.type) should be 901, 902, 903, or 904. Please refer to function BEQ_estimate')
end


%--------------------------------------------------------------------------
%% Step 1. Prepare data
%--------------------------------------------------------------------------

% Get data
Xm  = xest_out(:,1:Par.nM);                         % monthly data
XqY = xest_out(:,(size(xest_out,2)-Par.nQ+1):end);  % quarterly data
XqY = XqY(mod(datet_in(:,2),3)==0,:);               % selecting only 3rd months (otherwise is NaN)
dateQ = datet_in(mod(datet_in(:,2),3)==0,:);        % selecting also 3rd months for dates
Y = XqY(:,end);                                     % only target series
Xq = XqY(:,1:end-1);                                % only regressors

% Ensuring there is no NaN between two non-missing values - otherwise it cause cause the BVAR interpolation to fail
% Missing data in-between is filled via linear interpolation
% Not a frequent case (so far only one series in UK nowcast) but can be the case with do_Covid = 3
Xm = BEQ_inbetween_nan(Xm);
Xq = BEQ_inbetween_nan(Xq);
Y = BEQ_inbetween_nan(Y);

% Check that the number of observations is sufficient for interpolation
% For monthly data
min_obs = 2*Par_BVAR.lags + 1; % see BEQ_arfit (lines 86 to 88) for explanation on minimum number of observations
iMax = find(~any(isnan(Xm),2),1,'last'); % position of the last non-NaN
idx_keep_Xm = [];
for cc = 1:Par.nM        
    iX = find(~isnan(Xm(:,cc)),1,'first');
    if iMax - iX + 1 > min_obs
        idx_keep_Xm = [idx_keep_Xm, cc];
    end
end
Xm = Xm(:,idx_keep_Xm);
Par.nM = length(idx_keep_Xm); % changing number of variables
Par.trf.transf_m = Par.trf.transf_m(idx_keep_Xm); % adjusting transformations

% Same for quarterly data
if Par.nQ > 1
    iMax = find(~any(isnan(Xq),2),1,'last'); % position of the last non-NaN
    idx_keep_Xq = [];
    for cc = 1:(Par.nQ-1)        
        iX = find(~isnan(Xq(:,cc)),1,'first');
        if iMax - iX + 1 > min_obs
            idx_keep_Xq = [idx_keep_Xq, cc];
        end
    end
    Xq = Xq(:,idx_keep_Xq);
    Par.nQ = length(idx_keep_Xq) + 1; % still +1 for the target series
    Par.trf.transf_q = Par.trf.transf_q([idx_keep_Xq,end]); % adjusting transformations
end

% Check that the number of observations is sufficient for OLS regression
% For monthly series
XmQ = Xm(mod(datet_in(:,2),3)==0,:);        % selecting only 3rd months
iY = find(~isnan(Y),1,'last');
idx_keep_Xm = [];
for cc = 1:Par.nM        
    iX = find(~isnan(XmQ(:,cc)),1,'first');
    if iX + Par.lagM + Par.lagY + 2 < iY
        idx_keep_Xm = [idx_keep_Xm, cc];
    end
end
Xm = Xm(:,idx_keep_Xm);
Par.nM = length(idx_keep_Xm);
Par.trf.transf_m = Par.trf.transf_m(idx_keep_Xm); % adjusting transformations

% Same for quarterly series
if Par.nQ > 1
    idx_keep_Xq = [];
    for cc = 1:(Par.nQ-1)        
        iX = find(~isnan(Xq(:,cc)),1,'first');
        if iX + Par.lagQ + Par.lagY + 2 < iY
            idx_keep_Xq = [idx_keep_Xq, cc];
        end
    end
    Xq = Xq(:,idx_keep_Xq);
    Par.nQ = length(idx_keep_Xq) + 1; % still +1 for the target series
    Par.trf.transf_q = Par.trf.transf_q([idx_keep_Xq,end]); % adjusting transformations
end


% Add return if all is NaN
if Par.nM == 0
    disp('Not enough variables satisfying criteria, returning the initial dataset with NaN')
    Res.beq_combinations = [];
    Res.Y_fcst = [];
    Res.Y_fcst_indiv = NaN;
    Res.Date_fcst = [];
    Res.X_sm = xest_out;
    return
end


%--------------------------------------------------------------------------
%% Step 2. Prepare combinations of bridge equations
%--------------------------------------------------------------------------
% These are all regressions with either:
% - one monthly regressor (Par.nM)
% - two monthly regressors (binomial coefficient of 2 among Par.nM)
%   -> both of above with each possible quarterly regressor and none (above times Par.nQ)
%       -> then all are run with three different interpolation techniques (above times 3) identified as
%           o BVAR on all variables (type 901)
%           o BVAR only for selected variables (type 902)
%           o Univariate BVAR (type 903)
% Number of specifications is then 3 * Par.nQ * (Par.nM + (2 in Par.nM))

% Combinations of monthly regressors
beq_comb = nchoosek(1:Par.nM,2); % all combinations with 2 monthly regressors among Par.nM (bivariate)
temp = 1:Par.nM;
temp = [nan(Par.nM,1),temp']; % univariate regressions with 1 monthly regressor
beq_comb = [beq_comb;temp]; % merging

% Including quarterly regressors (if any) in the combinations
beq_comb_q = [beq_comb,nan(size(beq_comb,1),1)]; % adding a first column with NaN (meaning no quarterly regressor)
if Par.nQ > 1 % NB: by default Par.nQ is at least one (the quarterly target variable)
    for nn = 1:(Par.nQ-1)
        temp = [beq_comb,repmat(nn,size(beq_comb,1),1)];
        beq_comb_q = [beq_comb_q;temp];
    end
end

% Indexing on the types of interpolations (901 to 903)
% Structure of beq_comb_q
% col. 1 = type of interpolation
% cols. 2 & 3 = monthly regressors
% col. 4 = quarterly regressor (if any)
switch Par.type
    case 901
        beq_comb_q = [repmat(901,size(beq_comb_q,1),1),beq_comb_q];
    case 902
        beq_comb_q = [repmat(902,size(beq_comb_q,1),1),beq_comb_q];
    case 903
        beq_comb_q = [repmat(903,size(beq_comb_q,1),1),beq_comb_q];
    case 904
        beq_comb_q = [repmat(901,size(beq_comb_q,1),1),beq_comb_q; ...
                      repmat(902,size(beq_comb_q,1),1),beq_comb_q; ...
                      repmat(903,size(beq_comb_q,1),1),beq_comb_q];
end

% Add to Res structure
nEqs = size(beq_comb_q,1);
Res.beq_combinations = beq_comb_q;

% Prepare container for coefficients
ncol_max = 1 + 2*(Par.lagM + 1) + (Par.lagQ + 1) + Par.lagY + size(Par.Dum,1);
Res.coeffs = nan(nEqs,ncol_max);

% Allocate coefficients if they are not re-estimated
if re_estimate == 0
    if size(coeffs_in,1) ~= nEqs || size(coeffs_in,2) ~= ncol_max
        error('Not enough coefficients passed into the function to compute the forecasts.')
    end 
elseif re_estimate == 1
    coeffs_in = zeros(nEqs,ncol_max);    
else
    error('Parameter re_estimate should take values 0 or 1.')
end


%--------------------------------------------------------------------------
%% Step 3. Running each bridge equation
%--------------------------------------------------------------------------

% Loop over bridge equations
for i = 1:nEqs

    if mod(i,50)==0
        disp(['Doing iteration ',num2str(i),' out of ',num2str(nEqs)])
    end

    %% A. Interpolate data

    if beq_comb_q(i,1) == 901 % interpolate with BVAR on all data

        if i == 1 % do once, only for first i (otherwise take previous)

            % Interpolate monthly data
            if size(Xm,2)>1
                Par_BVAR.lambda = 0.2; % shrinkage parameters for multivariate BVAR (as in Banbura et al., 2023)
            else
                Par_BVAR.lambda = 0.5; % shrinkage parameters for univariate BVAR (as in Banbura et al., 2023)
            end
            Xm_ext = BEQ_Run_extrapolation_BVAR(Xm,Par_BVAR);

            % Interpolate quarterly data
            if Par.nQ > 1 % if we have more than the quarterly target variable
                if size(Xq,2)>1
                    Par_BVAR.lambda = 0.2; % shrinkage parameters for multivariate BVAR (as in Banbura et al., 2023)
                else
                    Par_BVAR.lambda = 0.5; % shrinkage parameters for univariate BVAR (as in Banbura et al., 2023)
                end
                Xq_ext = BEQ_Run_extrapolation_BVAR(Xq,Par_BVAR);
            end
        end
    
    elseif beq_comb_q(i,1) == 902 % only selection here, interpolation is done below

        Xm_ext = Xm;
        if Par.nQ > 1 % if we have more than the quarterly target variable
            Xq_ext = Xq;
        end

    elseif beq_comb_q(i,1) == 903 % interpolate with BVAR on all data

        Par_BVAR.lambda = 0.5; % shrinkage parameters for univariate BVAR (as in Banbura et al., 2023)

        if i == 1 || beq_comb_q(i-1,1) ~= 903 % do once, when going to the univariate BVAR regression (otherwise, take previous)

            Xm_ext = NaN(size(Xm)); % container
            for cc = 1:size(Xm,2)
                Xm_ext(:,cc) = BEQ_Run_extrapolation_BVAR(Xm(:,cc),Par_BVAR); % interpolate one by one
            end

            if Par.nQ > 1 % if we have more than the quarterly target variable
                Xq_ext = NaN(size(Xq)); % container for quarterly variables
                for cc = 1:size(Xq,2)
                    Xq_ext(:,cc) = BEQ_Run_extrapolation_BVAR(Xq(:,cc),Par_BVAR);
                end
            end
            
        end

    end


    %% B. Prepare data

    % Select quarterly variables
    sel_q = beq_comb_q(i,4);
    sel_q = sel_q(~isnan(sel_q));
    if ~isempty(sel_q)
        Xq_ext_i = Xq_ext(:,sel_q);
    else
        Xq_ext_i = [];
    end

    % Select monthly variables
    sel_m = beq_comb_q(i,2:3);
    sel_m = sel_m(~isnan(sel_m));
    Xm_ext_i = Xm_ext(:,sel_m);
    transf_m_beq = Par.trf.transf_m(sel_m); % also do for transformations (for transform_data function)

    % If 902 (interpolation with selected variables), interpolate here
    if beq_comb_q(i,1) == 902

        % Monthly data
        if size(Xm,2)>1
            Par_BVAR.lambda = 0.2; % shrinkage parameters for multivariate BVAR (as in Banbura et al., 2023)
        else
            Par_BVAR.lambda = 0.5; % shrinkage parameters for univariate BVAR (as in Banbura et al., 2023)
        end
        Xm_ext_i = BEQ_Run_extrapolation_BVAR(Xm_ext_i,Par_BVAR);

        % Quarterly data (if any)
        if ~isempty(sel_q)
            Par_BVAR.lambda = 0.5; % shrinkage parameters for univariate BVAR (as in Banbura et al., 2023)
            Xq_ext_i = BEQ_Run_extrapolation_BVAR(Xq_ext_i,Par_BVAR);
        end
        
    end

    
    %% C. Run forecast

    % Get the coeffs
    % NB: are used only if re_estimate == 0
    coeffs_in_i = coeffs_in(i,:);
    coeffs_in_i = coeffs_in_i(~isnan(coeffs_in_i));

    % Do the forecasts
    [Res.Y_fcst_indiv(:,i),Res.Date_fcst,ctb_i,coeffs] = ...
        BEQ_forecast(Xm_ext_i,datet_in,Xq_ext_i,Y,dateQ,Par,transf_m_beq,re_estimate,coeffs_in_i);

    % Add the coefficients
    Res.coeffs(i,1:length(coeffs)) = coeffs';

    % Check coefficients (if re_estimate == 0)
    if re_estimate == 0
        if ~isequal(coeffs',coeffs_in_i)
            error('Even if re_estimate was set to 0, coefficients in and out of BEQ_forecast are different.')
        end
    end


    %% D. Organize contributions

    % Prepare container for contributions
    % NB: needs to be here to get the number of observations
    % Dim. 1 = time
    % Dim. 2 = variables (including constant and dummies on the right)
    % Dim. 3 = equation
    if i == 1
        Res.contrib = nan(size(ctb_i,1),(Par.nM + Par.nQ + 1 + 1*(size(Par.Dum,1)>0)),nEqs); % for 2nd dimension +1 for the constant
    end

    % Add a container with names for the contributions
    if i == 1

        % Container
        Res.contrib_names = repmat({''},1,size(Res.contrib,2));
        
        % Monthly and quarterly variables (incl. target at the end)
        for cc = 1:Par.nM+Par.nQ
            Res.contrib_names{cc} = nameseries{cc};
        end
        
        % Constant
        Res.contrib_names{Par.nM+Par.nQ+1} = 'Constant';
        
        % Dummies (if any)
        if size(Par.Dum,1) > 0
            Res.contrib_names{Par.nM+Par.nQ+1+1} = 'Dummies';
        end

    end

    % Allocate contribution of the constant (cflag = 1 in BEQ_forecast)
    Res.contrib(:,(Par.nM + Par.nQ + 1),i) = ctb_i(:,1);
    
    % Allocate contributions of monthly variables
    idx_col = 2;
    if ~isnan(beq_comb_q(i,2))
        Res.contrib(:,beq_comb_q(i,2),i) = ctb_i(:,idx_col);
        idx_col = idx_col + 1;
    end
    if ~isnan(beq_comb_q(i,3))
        Res.contrib(:,beq_comb_q(i,3),i) = ctb_i(:,idx_col);
        idx_col = idx_col + 1;
    end

    % Allocate contribution of quarterly variable (if any)
    if ~isnan(beq_comb_q(i,4)) % this is a quarterly variable
        Res.contrib(:,(Par.nM + beq_comb_q(i,4)),i) = ctb_i(:,idx_col);
        idx_col = idx_col + 1;
    end

    % Allocate contribution of lagged endogenous variable (if any)
    if Par.lagY > 0
        Res.contrib(:,(Par.nM + Par.nQ),i) = ctb_i(:,idx_col);
        idx_col = idx_col + 1;
    end

    % Allocate contributions of dummies (if any)
    if size(Par.Dum,1) > 0
        Res.contrib(:,(Par.nM + Par.nQ + 1 + 1),i) = ctb_i(:,idx_col);
    end

end % end of loop on nEqs


%--------------------------------------------------------------------------
%% Step 4. Obtain nowcasts from combined bridge equations
%--------------------------------------------------------------------------

% Prepare results (average over all bridge equations)
% Res.Y_fcst = mean(Res.Y_fcst_indiv,2,'omitnan');
Res.contrib(isnan(Res.contrib)) = 0;
ctb_tot = mean(Res.contrib,3,'omitnan');

% Prepare results (median over all bridge equations)
Res.Y_fcst = median(Res.Y_fcst_indiv,2,'omitnan');

% Rescale contributions so they match the median
% To be commented out if the averaging method is not the median
for ii = 1:length(Res.Date_fcst)

    % Get index in datet_in and check if missing target variable
    idx_datet_in = Res.Date_fcst(ii,1)==datet_in(:,1) & Res.Date_fcst(ii,2)==datet_in(:,2);
    is_missing_Y = isnan(xest_out(idx_datet_in,end));

    % Rescale contributions if Y (the target variable) is missing
    if is_missing_Y

        % Get the difference between contributions and predictions
        sum_ctb = sum(ctb_tot(ii,:),2);
        disc = Res.Y_fcst(ii) - sum_ctb;

        % Create rescaled contriutions
        ctb_add_on = ctb_tot(ii,:)/sum_ctb*disc;
        ctb_tot(ii,:) = ctb_tot(ii,:) + ctb_add_on;

        % Check that it worked
        if round(Res.Y_fcst(ii),2) ~= round(sum(ctb_tot(ii,:)),2) && ~isnan(Res.Y_fcst(ii)) && ~isnan(sum(ctb_tot(ii,:)))
            error(['Contributions do not match forecast for observation ',num2str(ii),' in datet_in.'])
        end

    end

end

% Write in Res.X_sm (needed for the results) and Res.contrib_X_sm
Res.X_sm = xest_out;
Res.contrib_X_sm = nan(size(xest_out,1),size(ctb_tot,2));
for ii = 1:length(Res.Date_fcst)
    jj = find(datet_in(:,1)==Res.Date_fcst(ii,1) & datet_in(:,2)==Res.Date_fcst(ii,2));
    Res.X_sm(jj,end) = Res.Y_fcst(ii);
    Res.contrib_X_sm(jj,:) = ctb_tot(ii,:);
end  

end