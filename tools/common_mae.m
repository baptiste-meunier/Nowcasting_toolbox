function [MAE] = common_mae(xest,Par,t_m,m,datet,do_Covid,country,MAE,do_mae,nameseries)
% Calcuates the mean absolute error (MAE) for the time the nowcast is run based
% on last 10 years. E.g. when run in the first month of a quarter, will calculate error
% based on past predictions for every first month of each quarter for past 10 years.

% INPUTS 
% - xest [matrix] = values of the input series
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
%       o type [scalar] = type of interpolation performed (see BEQ_estimate for details)
%       o Dum [matrix] = dates of the dummies (year, month). Format should be a k x 2 matrix with k = number of dummies (each dummy on a row) and year / month in columns (BEQ)
%       o bvar_lags [scalar] = number of lags (BVAR)
%       o blocks_full [vector/matrix] = block identification (for do_loop=2)
%       o size_data_in [scalar] = lengths of raw data set
% - t_m [matrix] = dates (year / month)
% - m [scalar] = number of months ahead
% - datet [matrix] = dates (year / month). Difference with t_m is that datet relates to the data when adding NaN for months ahead
% - do_Covid [scalar] = method for Covid correction 
% - country [structure] = country and model parameters
%       o name [string] = name of the country
%       o model [string] = model type
% - MAE [structure] mean absolute errors
%       o Bac [structure] MAE for backcasting
%           > mae_1st [scalar] = MAE for predictions in the 1st month of the quarter (value set by the user)
%           > mae_2nd [scalar] = MAE for predictions in the 2nd month of the quarter (value set by the user)
%           > mae_3rd [scalar] = MAE for predictions in the 3rd month of the quarter (value set by the user)
%       o Now [structure] MAE for nowcasting - same content as Bac
%       o For [structure] MAE for forecasting - same content as Bac
% - do_mae [scalar] = switch on how to compute the MAE
% - nameseries [cell array] = name of the series (short version)
%
% OUTPUT
% - MAE [structure] mean absolute errors
%       o Bac [structure] MAE for backcasting
%           > Pred [matrix] out-of-sample real-time predictions
%               + col. 1 = Real-time year
%               + col. 2 = Real-time month
%               + col. 3 = Predicted year
%               + col. 4 = Predicted month
%               + col. 5 = Prediction (out-of-sample, real-time)
%               + col. 6 = Actual
%           > mae [scalar] = Mean Abolute Error, adjusted for outliers
%           > mae_unadjusted [scalar] = Mean Abolute Error (unadjusted)
%           > mae_1st [scalar] = MAE for predictions in the 1st month of the quarter (value set by the user)
%           > mae_2nd [scalar] = MAE for predictions in the 2nd month of the quarter (value set by the user)
%           > mae_3rd [scalar] = MAE for predictions in the 3rd month of the quarter (value set by the user)
%       o Now [structure] MAE for nowcasting - same content as Bac
%       o For [structure] MAE for forecasting - same content as Bac
%

if do_mae == 1

    disp('Section 7: Error evaluation starts');

    % ---------------------------------------------------------------------
    %% A. INITIALIZATION
    % ---------------------------------------------------------------------

    % Take current month as month of nowcast
    MAE.nowcast_date = [year(datetime("today")), month(datetime("today"))];

    % Set start and end of evaluation 
    end_eval = find(t_m(:,1)==year(datetime("today")) & t_m(:,2)==month(datetime("today")))-3; % end_eval is three months before
    if end_eval > 132
        start_eval = end_eval - 120; % 10 years before end_eval
    else
        error("At least 11 years of data are needed to run the MAE")
    end

    % Structure of missing values (NaN)
    xact = xest(1:end-m,:); % removing last m NaN
    [nb_obs,nb_var] = size(xact);
    MAE.delay = zeros(1,nb_var);
    for nn = 1:nb_var
        last_non_nan = find(~isnan(xact(:,nn)), 1, 'last');
        MAE.delay(1,nn) = nb_obs - last_non_nan;
    end
    
    % Save initial parameters (for Covid correct and NaN correction)
    nM_init = Par.nM_init;
    blocks_init = Par.blocks;
    r_init = Par.r;
    transf_m_init = Par.trf.transf_m;
    transf_q_init = Par.trf.transf_q;

    C = {'Bac','Now','For'};

    % Initiate containers for predictions and directions
    for hh = 1:numel(C)
        hor = C{hh};
        MAE.(hor).Pred = [];
        MAE.(hor).Dirc = [];
    end

    % ---------------------------------------------------------------------
    %% B. LOOP OVER OBSERVATIONS
    % ---------------------------------------------------------------------

    % Loops run to estimate model at respective month of each quarter for past 10 years
    % Steps of 3 to only estimate the relevant months
    for i=start_eval:3:end_eval

        disp(['Out-of-sample predictions for year ', num2str(t_m(i,1)), ' month ', num2str(t_m(i,2))])

        % Re-create data that would have been available to a real-time forecaster (with NaN)
        xin = xest(1:i,:); % data until month i
        datet_in = datet(1:(i+m),:);
        for j = 1:size(xest,2)
            if MAE.delay(j) ~= 0
                xin(end-MAE.delay(j)+1:end,j) = NaN; % delete last k points to reproduce missing observations
            end
        end
        yest = [xin; nan(m,size(xin,2))]; % adding the m NaN at the end for prediction

        % Compute number of months till end of quarter 
        % Should be i+2 in first month, i+1 in second month, and i+0 in third month
        iQ = i + mod(3-mod(t_m(i,2),3),3);

        % Create indicator if target is available
        % NB: This is in case this is deleted by do_Covid == 2 or 3
        is_target_bac = ~isnan(yest(iQ-3,end));

        % Correct for Covid and NaN      
        [yest_out,nM_out,blocks_out,r_out,~,~,~,~,transf_m_out,transf_q_out] = ...
            common_NaN_Covid_correct(yest,datet_in,do_Covid,nM_init,blocks_init,r_init, ...
                                     ones(1,size(yest,2)),{'Placeholder'},repmat({'name'},1,size(yest,2)), ...
                                     repmat({'name'},1,size(yest,2)), ...
                                     transf_m_init,transf_q_init); % because in this case we don't need the groups name and ID
        Par.nM = nM_out;
        Par.nQ = size(yest_out,2) - Par.nM;
        Par.blocks = blocks_out;
        Par.r = r_out;
        Par.trf.transf_m = transf_m_out;
        Par.trf.transf_q = transf_q_out;

        % Estimate model
        switch country.model
            case 'DFM'
                Res = DFM_estimate(yest_out,Par);
            case 'BEQ'
                Res = BEQ_estimate(yest_out,Par,datet_in,nameseries,1,[]);
            case 'BVAR'
                Res = BVAR_estimate(yest_out,Par,datet_in);
        end
      
        % Write results for point forecasts
        % 1st column = year of real-time
        % 2nd column = month of real-time
        % 3rd column = year of back/now/fore-cast
        % 4th column = month of back/now/fore-cast
        % 5th column = DFM prediction
        % 6th column = actual value
        MAE.Bac.Pred = [MAE.Bac.Pred; datet_in(i,1:2) datet_in(iQ-3,1:2) Res.X_sm(iQ-3,end) xest(iQ-3,end)];
        if is_target_bac % manual fix for backcasting if actual data is deleted (the reason is that it can be deleted with do_Covid == 2 or 3)
            MAE.Bac.Pred(end,5) = xest(iQ-3,end);
        end
        MAE.Now.Pred = [MAE.Now.Pred; datet_in(i,1:2) datet_in(iQ,1:2) Res.X_sm(iQ,end) xest(iQ,end)];
        MAE.For.Pred = [MAE.For.Pred; datet_in(i,1:2) datet_in(iQ+3,1:2) Res.X_sm(iQ+3,end) xest(iQ+3,end)];

        % Write results for directional accuracy
        % 1st column = year of real-time
        % 2nd column = month of real-time
        % 3rd column = year of back/now/fore-cast
        % 4th column = month of back/now/fore-cast
        % 5th column = DFM direction
        % 6th column = actual direction
        MAE.Bac.Dirc = [MAE.Bac.Dirc; datet_in(i,1:2) datet_in(iQ-3,1:2) Res.X_sm(iQ-3,end)>Res.X_sm(iQ-6,end) xest(iQ-3,end)>xest(iQ-6,end)];
        if is_target_bac % manual fix for backcasting if actual data is already there (the reason is that it can be deleted with do_Covid == 2 or 3)
            MAE.Bac.Dirc(end,5) = xest(iQ-3,end)>xest(iQ-6,end);
        end
        MAE.Now.Dirc = [MAE.Now.Dirc; datet_in(i,1:2) datet_in(iQ,1:2) Res.X_sm(iQ,end)>Res.X_sm(iQ-3,end) xest(iQ,end)>xest(iQ-3,end)];
        MAE.For.Dirc = [MAE.For.Dirc; datet_in(i,1:2) datet_in(iQ+3,1:2) Res.X_sm(iQ+3,end)>Res.X_sm(iQ,end) xest(iQ+3,end)>xest(iQ,end)];

    end % end of loop on observations

    % ---------------------------------------------------------------------
    %% C. CALCULATE ERRORS
    % ---------------------------------------------------------------------

    % Calculate mean absolute error for all horizons
    for hh = 1:numel(C)
        hor = C{hh};
        if hor == "For"
            MAE.(hor).Pred=MAE.(hor).Pred(sum(isnan(MAE.(hor).Pred),2)==0,:); % gets rid of NaN for forecast
        end
        errors_unadjusted = abs(MAE.(hor).Pred(:,6) - MAE.(hor).Pred(:,5));
        [errors_adjusted,~] = common_outliers(errors_unadjusted,0); % produces NaN for outliers
        MAE.(hor).mae = mean(errors_adjusted,"omitnan"); % omit outliers in calculation
        MAE.(hor).mae_unadjusted = mean(errors_unadjusted,"omitnan");
    end

    % Calculate forecast directional accuracy for all horizons
    for hh = 1:numel(C)
        hor = C{hh};
        MAE.(hor).fda = 1-mean(abs(MAE.(hor).Dirc(:,5) - MAE.(hor).Dirc(:,6)));
    end


    disp('Section 7: Error evaluation completed');

elseif do_mae == 0

    % Get the month of the quarter
    mth_qtr = mod(month(datetime("today")),3);  

    % Use errors pre-set by user in main file depending on month of quarter
    if mth_qtr == 1

        MAE.Bac.mae = MAE.Bac.mae_1st;
        MAE.Now.mae = MAE.Now.mae_1st;
        MAE.For.mae = MAE.For.mae_1st;
        MAE.Bac.fda = MAE.Bac.fda_1st;
        MAE.Now.fda = MAE.Now.fda_1st;
        MAE.For.fda = MAE.For.fda_1st;

    elseif mth_qtr == 2

        MAE.Bac.mae = MAE.Bac.mae_2nd;
        MAE.Now.mae = MAE.Now.mae_2nd;
        MAE.For.mae = MAE.For.mae_2nd;
        MAE.Bac.fda = MAE.Bac.fda_2nd;
        MAE.Now.fda = MAE.Now.fda_2nd;
        MAE.For.fda = MAE.For.fda_2nd;

    elseif mth_qtr == 0

        MAE.Bac.mae = MAE.Bac.mae_3rd;
        MAE.Now.mae = MAE.Now.mae_3rd;
        MAE.For.mae = MAE.For.mae_3rd;
        MAE.Bac.fda = MAE.Bac.fda_3rd;
        MAE.Now.fda = MAE.Now.fda_3rd;
        MAE.For.fda = MAE.For.fda_3rd;

    end

    disp('Section 7: Error evaluation completed, using values obtained in the 2023 review of models');

else

    error('do_mae should be either 0 or 1')

end % end of condition on do_error_eval

end