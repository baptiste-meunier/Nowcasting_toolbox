function [fcstnew, Date_fcst,ctb_new,coeffs] = BEQ_forecast(Xm_ext_i,datet_in,Xq_ext_i,Y,dateQ,Par,transf_m_beq,re_estimate,coeffs_in_i)
% This function estimates a bridge equation
%
% Code adapted from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - Xm_ext_i [matrix] = monthly regressors (including extrapolations)
% - dateM [matrix] = monthly dates
% - Xq [matrix] = quarterly regressors at quarterly frequency (including extrapolations)
% - Y [matrix] = dependent variable at quarterly frequency 
% - dateQ [matrix] = quarterly dates (year / month)
% - Par [structure] = settings of bridge equations
%       o lagM [scalar] = number of lags for monthly regressor(s) (in quarterly terms)
%       o lagQ [scalar] = number of lags for quarterly regressor(s) (in quarterly terms)
%       o lagY [scalar] = number of lags for the endogenous variable (in quarterly terms)
%       o Dum [matrix] = dates of the dummies (year / month)
%       o nDraw [scalar] = number of draws (for density forecasts)
% - transf_m_beq [vector] = transformations for monthly regressors
% - re_estimate [sclar] = whether to re-estimate OLS or use the coefficients passed to the functions
% - coeffs_in_i [matrix] = matrix of coefficients (used if the OLS is not re-estimated)
%
% OUTPUTS
% - fcstnew [vector] = forecasts of the dependent variable
% - Date_fcst [matrix] = quarterly dates (year / month)
% - ctb_new [matrix] = contributions to forecast
% - coeffs [matrix] = coefficients of OLS
% 

% -------------------------------------------------------------------------
%% Step 1: Aggregation of data at different frequencies
% -------------------------------------------------------------------------

% Convert monthly data back to levels (before aggregation)
Xm_levels = common_transform_data(Xm_ext_i,transf_m_beq,-1);

% Aggregate monthly data to quarterly
AggInd = repmat({'A'},1,size(Xm_levels,2)); % letter "A" means we take the average
if ~isempty(Xm_ext_i)
    [XmQ,dateMQ] = BEQ_m2q(Xm_levels,datet_in,AggInd);
else % where there are no monthly indicators
    XmQ = [];
end

% Re-transform monthly data (back to stationary data)
if ~isempty(XmQ)

    XmQ = common_transform_data(XmQ,transf_m_beq,1);
    dateMQ = dateMQ(2:end,:); % deleting first row (which is also deleted in transformation)

    % Check dimensions of XmQ and dateMQ
    if size(dateMQ,1)~=size(XmQ,1)
        error('Problem in the dimensions for the bridge equations (monthly dates and data aggregated to quarterly)') 
    end

    % Check that the first date of dateMQ corresponds to second of dateQ
    if dateMQ(1,:)~=dateQ(2,:)
        error('Problem in the dimensions for the bridge equations (dates)') 
    else

        % Delete first observation of XmQ (to avoid issues with first aggregated value)
        % Meaning we need to delete the first two observations of Xq and Y and dateQ (to match lengths)
        XmQ = XmQ(2:end,:);
        if ~isempty(Xq_ext_i)
            Xq_ext_i = Xq_ext_i(3:end,:);
        end
        Y = Y(3:end,:);
        dateQ = dateQ(3:end,:);

    end

    % Adjust if we have excess observations in XmQ
    % This can be the case due to the extension of monthly data in BEQ_m2q
    if size(dateMQ,1) > size(dateQ,1)
        XmQ = XmQ(1:size(dateQ,1),:);
    end

    % Check dimension with Y (target variable) and, when relevant, Xq (quarterly regressors)   
    if size(XmQ,1)~=size(Y,1)
        error('Problem in the dimensions for the bridge equations (size of XmQ and Y)') 
    end
    if ~isempty(Xq_ext_i)
        if size(XmQ,1)~=size(Xq_ext_i,1)
            error('Problem in the dimensions for the bridge equations (size of XmQ and Xq)')
        end
    end

end


% -------------------------------------------------------------------------
%% Step 2: Obtain regressors with lags
% -------------------------------------------------------------------------

% Container for regressors (with lags)
Xlags = [];

% Monthly indicators
if ~isempty(XmQ) % add monthly predictors (and lags) to the regressor matrix
    for i = 1:size(XmQ,2)
        Xlags = [Xlags,lagmatrix(XmQ(:,i),0:Par.lagM)];
    end
end

% Quarterly indicators
if ~isempty(Xq_ext_i) % add quarterly predictors (and lags) to the regressor matrix
    for i = 1:size(Xq_ext_i,2)
        Xlags = [Xlags,lagmatrix(Xq_ext_i(:,i),0:Par.lagQ)];         
    end
end 

% Lags of endogenous variable
if Par.lagY > 0  
    Xlags = [Xlags,lagmatrix(Y,1:Par.lagY)];     
end

nobs  = size(Xlags,1);

% Add dummies
if ~isempty(Par.Dum)
    for idums = 1:size(Par.Dum,1)

        % Create dummy
        Dummy = zeros(nobs,1);

        % Identify position of the 1
        % NB: if the date to be dummied is not yet passed, the series will
        %     be all 0 and will get discarded in the check below
        idx_dummy = find(dateQ(:,1)==Par.Dum(idums,1) & dateQ(:,2)==Par.Dum(idums,2));        
        if ~isempty(idx_dummy) % this controls for dummies set for future horizons            
            Dummy(idx_dummy,1) = 1;            
        end

        % Add dummy to Xlags
        Xlags = [Xlags Dummy];

    end
end


% -------------------------------------------------------------------------
%% Step 3: Forecast the quarterly variable using bridge equation
% -------------------------------------------------------------------------

% Check coefficients (if re_estimate == 0)
if re_estimate == 0
    if length(coeffs_in_i) ~= (size(Xlags,2) + 1) % +1 accomodates for the constant
        error('Not enough coefficients')
    end
end

% Check Xlags and delete series that are all NaN or all 0 for estimation sample
chk_Xlags = Xlags(~isnan(Y),:);
chk_Xlags(chk_Xlags == 0) = NaN; % replacing 0 with NaN to make it easier to check
Xlags = Xlags(:,~all(isnan(chk_Xlags))); % removing columns that are all NaN
idx_delete = find(all(isnan(chk_Xlags)));

% Correct coeffs_in_i (if re_estimate == 0)
if re_estimate == 0
    coeffs_chk = coeffs_in_i(2:end);
    coeffs_chk = coeffs_chk(~all(isnan(chk_Xlags))); % removing coefficients for columns that are deleted
    coeffs_in_i = [coeffs_in_i(1),coeffs_chk];
end

% Estimate bridge equation
sflag = 1;
cflag = 1;
if re_estimate == 0
    beta = coeffs_in_i';
    sig_beq = 1;
elseif re_estimate == 1   
    [beta,sig_beq] = BEQ_OLS(Y,Xlags,cflag,sflag);
else
    error('Parameter re_estimate should only take values in 0 and 1.')
end

if Par.nDraw == 0
    e = zeros(nobs,1);
else
    % random draw of the residuals
    e = randn(nobs,1)*sig_beq;
end

% Forecast
fcst             = [ones(nobs,1) Xlags] * beta + e;
fcst(~isnan(Y))  = Y(~isnan(Y)); % replacing with actual data (when available)

% Create contributions
ctb = [ones(nobs,1) Xlags].*repmat(beta',nobs,1);
ctb(isnan(ctb)) = 0; % replacing NaN with 0

% Get coefficients
coeffs = beta;


% -------------------------------------------------------------------------
%% Step 4: Address NaN for lagged dependent variable
% -------------------------------------------------------------------------

% In case lagged dependent variable is included in the equation
if Par.lagY > 0

    nL = find(isnan(Xlags(end,:)),1,'first'); % position of last non-lag of Y in Xlags (this will be the first to be NaN)
    nLD = Par.lagY;
    iLast = find(~isnan(fcst),1,'last'); % position for last non-NaN

    % Add lags of the dependent variable sequentially (using the previous period forecast)
    for j = iLast+1:nobs

        tmp = nan(nobs,Par.lagY); % temp file for prediction
        i = 0; % counter

        for jlag = 1:Par.lagY

            % Add to temp
            i = i+1;
            tmp(1:j,i) = [nan(jlag,1);fcst(1:j-jlag)]; % lagmatrix does not work properly here as it loses last observations

            % Add to contribution (only if data is missing)
            if isnan(Xlags(j,nL+jlag-1))
                ctb_add = ctb(j-jlag,:) * beta(nL+jlag);
                ctb(j,:) = ctb(j,:) + ctb_add;
            end

        end

        % Compute new forecast
        Xlags(j,nL:(nL+nLD-1)) = tmp(j,:); % replacing in Xlags by forecast
        fcst(j,:) = [1 Xlags(j,:)] * beta + e(j); % adding a forecast      

    end
    
end

% Check contributions and forecast
% Only when there is no data for Y - otherwise, there is the contribution of the error term
for rr = find(~isnan(Y),1,'last')+1:nobs
    if (round(fcst(rr),2) ~= round(sum(ctb(rr,:)),2)) && ~isnan(sum(ctb(rr,:))) && ~isnan(fcst(rr)) && (fcst(rr) < 100) && (fcst(rr) > -100)
        error(['Contributions do not match forecast for observation ',num2str(rr)])
    end
end

% Output
fcstnew = fcst;
Date_fcst = dateQ;


% -------------------------------------------------------------------------
%% Step 5: Arrange contributions and coeffs with missing variable
% -------------------------------------------------------------------------
% See above chk_Xlags can end up deleting some variables that are 0 / NaN

% Correct for deleted variables
for ii = 1:length(idx_delete)
    idx_ii = idx_delete(ii);
    if idx_ii == 1
        coeffs = [0;coeffs];
        ctb = [zeros(size(ctb,1),1),ctb];
    elseif idx_ii == length(coeffs)+1
        coeffs = [coeffs;0];
        ctb = [ctb,zeros(size(ctb,1),1)];
    elseif idx_ii > length(coeffs)+1
        error('Size of coefficients not appropriate. Please check.')
    else
        coeffs = [coeffs(1:idx_ii-1);0;coeffs(idx_ii:end)];
        ctb = [ctb(:,1:idx_ii-1),zeros(size(ctb,1),1),ctb(:,idx_ii:end)];
    end
end

% Check size for coeffs
ncols = cflag + size(Xm_ext_i,2)*(Par.lagM+1) + size(Xq_ext_i,2)*(Par.lagQ+1) + Par.lagY + size(Par.Dum,1);
if length(coeffs) ~= ncols
    error('Problem for sizes in recovering coefficients. Please check.')
end

% Check size for ctb
if size(ctb,2) ~= ncols
    error('Problem for sizes in recovering contributions. Please check.')
end


% -------------------------------------------------------------------------
%% Step 6: Arrange contributions by variable
% -------------------------------------------------------------------------

% Sum contributions by variable
nvars = cflag + size(Xm_ext_i,2) + size(Xq_ext_i,2) + 1*(Par.lagY > 0) + 1*(size(Par.Dum,1) > 0);
ctb_new = nan(size(ctb,1),nvars);
new_col = 1;
old_col = 1;

% Contribution from constant (where applicable)
% NB: constant will be in the first column
if cflag == 1
    ctb_new(:,new_col) = ctb(:,old_col);
    new_col = new_col + 1;
    old_col = old_col + 1;
end

% Contributions from monthly variables
if size(Xm_ext_i,2) > 0

    for ii = 1:size(Xm_ext_i,2)

        ctb_new(:,new_col) = sum(ctb(:,old_col:(old_col+Par.lagM)),2,"omitnan");
        new_col = new_col + 1;
        old_col = old_col + 1 + Par.lagM;

    end
end

% Contributions from quarterly variables
if size(Xq_ext_i,2) > 0

    for ii = 1:size(Xq_ext_i,2)

        ctb_new(:,new_col) = sum(ctb(:,old_col:(old_col+Par.lagQ)),2,"omitnan");
        new_col = new_col + 1;
        old_col = old_col + Par.lagQ;

    end
end

% Contributions from lagged endogenous
if Par.lagY > 0

    ctb_new(:,new_col) = sum(ctb(:,old_col:(old_col+Par.lagY-1)),2,"omitnan");
    new_col = new_col + 1;
    old_col = old_col + Par.lagY;

end

% Contributions from dummies
% NB: all contributions from dummies are put into a common contribution
if size(Par.Dum,1) > 0

    for ii = 1:size(Par.Dum,1)
            
        ctb_new(isnan(ctb_new(:,new_col)),new_col) = 0;
        ctb_new(:,new_col) = ctb_new(:,new_col) + ctb(:,old_col);
        old_col = old_col + 1;

    end
end


end