function[Loop,Eval] = common_eval_models(do_loop,Loop,Eval,xest,Par,t_m,m,country,datet,do_Covid,groups)
% This script runs the evaluation of models 
% 
% INPUTS 
% - do_loop [scalar 0/1] = switch on whether to loop on models
% - Loop [structure] = user parameters of the loop
%       o n_iter [scalar] = number of random models to be tested      
%       o min_p [scalar] = minimum number of lags
%       o max_p [scalar] = maximum number of lags (NB: p>5 won't work)
%       o min_r [scalar] = minimum number of factors
%       o max_r [scalar] = maximum number of factors
%       o min_startyear [scalar] = minimum start year of the estimation sample
%       o max_startyear [scalar] = maximum start year of the estimation sample
%       o startmonth [scalar] = start month (same for any start year)
%       o min_var [scalar] = minimum number of variables
%       o max_var [scalar] = maximum number of variables (should not exceed the actual number of variables in load_data)
%       o name_loop [string] = name for the loop
%       o excel_loopfile [string] = name for the output Excel file
%       o do_random [scalar 0/1] = switch on whether to randomize the models in the 
%       o mdl_list [table] = list of pre-defined models for evaluation
% - Eval [structure] = parameters of the evaluation
%       o eval_startyear [scalar] = 2022;
%       o eval_startmonth [scalar] = 1;
%       o eval_endyear [scalar] = 2022;
%       o eval_endmonth [scalar] = 6;
%       o gdp_rel [scalar] = month of quarter at which GDP is available 
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
%       o Dum [matrix] = dates of the dummies (year, month). Format should be a k x 2 matrix with k = number of dummies (each dummy on a row) and year / month in columns (BEQ)
%       o bvar_lags [scalar] = number of lags (BVAR)
%       o blocks_full [vector/matrix] = block identification (for do_loop=2)
%       o size_data_in [scalar] = lengths of raw data set
% - t_m [matrix] = dates (year / month)
% - m [scalar] = number of months ahead
% - country [structure] = country parameters
%       o name [string] = name of the country
%       o stat [string] = transformation of target data (qoq or qoqAR)
%       o model [string] = model type
% - datet [matrix] = dates (year / month). Difference with t_m is that datet relates to the data when adding NaN for months ahead
% - do_Covid [scalar] = method for Covid correction 
% - groups [vector] = group ID of each variable
%
% OUTPUTS
% - Loop [structure] = results of the loop
% - Eval [structure] = results of the evaluation


% ---------------------------------------------------------------------
%% 0. CHECK HYPER-PARAMETERS OF LOOP
% ---------------------------------------------------------------------

if do_loop == 1

    % Check 1 - min_start date
    first_year = datet(1,1);
    first_month = datet(1,2);
    if (first_year > Loop.min_startyear) || (first_year == Loop.min_startyear && first_month > Loop.startmonth)
        error(['Loop: the minimum start date (Loop.min_startyear and Loop.startmonth) is before the starting date of the dataset = ',num2str(first_year),' M',num2str(first_month),'. Either set an earlier date for loading data (Par.startyear) or a later date for minimum estimation start date (Loop.min_startyear)'])
    end
    
    % Check 2 - number of variables
    if Loop.max_var > (Par.nM + Par.nQ)
        error(['Loop: the maximum number of variables (Loop.max_var) is higher than the number of variables (Par.nM + Par.nQ) in the dataset = ',num2str(Par.nM+Par.nQ)])
    end
    
    % Check 3 - number of lags in DFM
    if strcmp(country.model,'DFM') && Loop.max_p > 5
        error('Loop: the maximum number of lags (Loop.max_p) is higher than 5, code might not work. Please set it lower.')
    end

    % Check 4 - min start data and Eval data
    if (Eval.eval_startyear < Loop.min_startyear) || (Eval.eval_startyear == Loop.min_startyear && Eval.eval_startmonth < Loop.startmonth)
        error(['Loop: the minimum start date (Loop.min_startyear and Loop.startmonth) is before the starting date of the evaluation = ',num2str(Eval.eval_startyear),' M',num2str(Eval.eval_startmonth),'. Either set an earlier date for minimum estimation start date (Loop.min_startyear) or a later date for evaluation sample (Eval.eval_startyear)'])
    end

    % Check 5 - sufficient timespan for BEQ and BVAR
    if strcmp(country.model,'BEQ') || strcmp(country.model,'BVAR')
        if (Eval.eval_startyear - Loop.max_startyear < 3)
            error(['Loop: the maximum start date (Loop.max_startyear) is less than three years before starting year of the evaluation = ',num2str(Eval.eval_startyear),'. This would lead to a high proporition of non-feasible models. It is advised to set either an earlier date for maximum estimation start date (Loop.max_startyear) or a later date for evaluation sample (Eval.eval_startyear) with at least 3 years between the two.'])
        end
    end

end


% ---------------------------------------------------------------------
%% 1. PREPARE EVAL STRUCTURE
% ---------------------------------------------------------------------

% Get the names (horizons) of the eval structure
C = {'Bac','Now','For'};

if do_loop == 0

    Loop.n_iter = 1;
    excel_evalfile = strcat('./eval/',country.name,'/',country.name,'_',country.model,'_evaluation.xlsx'); % Excel file for evaluation metrics

elseif do_loop == 1

    % Store initial parameters (as they are modified at each iteration)
    Loop.xest = xest;
    Loop.nM = Par.nM;
    Loop.nQ = Par.nQ;
    Loop.t_m = t_m;
    Loop.datet = datet;
    Loop.n_var_xest = size(Loop.xest,2) - 1; % number of regressors. -1 because the last variable is the target (quarterly GDP growth)
    Loop.Dum = Par.Dum;

    % Save the structure of r and blocks, and transformations
    % NB: r is modified in random loops, there it is to save r for blocks
    % NB: transformations are used only for bridge equations
    Loop.blocks = Par.blocks;
    Loop.r = Par.r;
    Loop.trf.transf_m = Par.trf.transf_m;
    Loop.trf.transf_q = Par.trf.transf_q;
    
    % Containers for model parameters and results
    Loop.var_sel = cell(Loop.n_iter,1);
    Loop.groups_sel = cell(Loop.n_iter,1);
    [Loop.parameters,Loop.parameters_name,Loop.results_names] = common_initiate_loop(Loop.n_iter,country.model);
    for hh = 1:numel(C)
        hor = C{hh};
        Loop.(hor).results = nan(Loop.n_iter,40); % number of columns should correspond to number of RMSE / FDA entered (section 2.D)
    end

elseif do_loop == 2

    % Store initial parameters (as they are modified at each iteration)
    Loop.xest = xest;
    Loop.nM = Par.nM;
    Loop.nQ = Par.nQ;
    Loop.t_m = t_m;
    Loop.datet = datet;
    Loop.n_var_xest = size(Loop.xest,2) - 1;
    Loop.Dum = Par.Dum;

    % Save the structure of r and blocks, and transformations
    % NB: r is modified in random loops, there it is to save r for blocks
    % NB: transformations are used only for bridge equations
    Loop.blocks = Par.blocks_full;
    Loop.r = Par.r;
    Loop.trf.transf_m = Par.trf.transf_m;
    Loop.trf.transf_q = Par.trf.transf_q;

    % Set number of iterations to number of models in predefined list
    Loop.n_iter = height(Loop.mdl_list); 
    
    % Containers for model parameters and results
    Loop.var_sel = cell(Loop.n_iter,1);
    Loop.groups_sel = cell(Loop.n_iter,1);
    [Loop.parameters,Loop.parameters_name,Loop.results_names] = common_initiate_loop(Loop.n_iter,country.model);
    for hh = 1:numel(C)
        hor = C{hh};
        Loop.(hor).results = nan(Loop.n_iter,40); % number of columns should correspond to number of RMSE / FDA entered (section 2.D)
    end

end


% ---------------------------------------------------------------------
%% 2. RUN THE LOOP
% ---------------------------------------------------------------------

for n_iter_mod = 1:Loop.n_iter

    % ---------------------------------------------------------------------
    %% A. INITIALISE MODEL (RANDOM IF DO_LOOP == 1)
    % ---------------------------------------------------------------------

    if do_loop == 1 % then we randomize the hyper-parameters of the model

        disp(['LOOP: DOING ITERATION ', num2str(n_iter_mod), ' OUT OF ', num2str(Loop.n_iter)])

        % Randomize hyper-parameters
        [rnd_parameters,sel_var] = common_randomize_model(n_iter_mod,Loop,country.model,do_Covid,Par);
        Loop.parameters(n_iter_mod,:) = rnd_parameters;
        Loop.n_var_iter = Loop.parameters(n_iter_mod,4);

        % Adjust model hyper-parameters
        switch country.model
            case 'DFM'
                Par.r = Loop.r;
                Par.r(1) = rnd_parameters(5);
                Par.p = rnd_parameters(6);
            case 'BEQ'
                Par.r = Loop.r; % even though there are no factors, this is to avoid bug in NaN_Covid_correct function
                Par.lagM = rnd_parameters(5);
                Par.lagQ = rnd_parameters(6);
                Par.lagY = rnd_parameters(7);
                Par.Dum = Loop.Dum;
            case 'BVAR'
                Par.r = Loop.r; % even though there are no factors, this is to avoid bug in NaN_Covid_correct function
                Par.bvar_lags = rnd_parameters(5);
        end

        % Adjust the sample dates with new startyear
        startyear = rnd_parameters(3);
        startmonth = Loop.startmonth;
        idx_keep = find(Loop.t_m(:,1)==startyear & Loop.t_m(:,2)==startmonth); 
        t_m = Loop.t_m(idx_keep:end,:); % cutting dates
        datet = Loop.datet(idx_keep:end,:); % cutting dates
        xest = Loop.xest(idx_keep:end,:); % cutting the data

        % Adjust the sample variables to sel_var
        Loop.var_sel{n_iter_mod,:} = num2str(sel_var);
        xest = xest(:,[sel_var,end]); % adjusting the sample
                                      % adding end to have the quarterly variable (always last variable)

        % Adjust other parameters
        Par.blocks = Loop.blocks([sel_var,end]',:); % blocks
        Par.nM = sum(sel_var<=Loop.nM); % changing number of monthly variables in parameters
                                        % all variables that are before Loop.nM
        Par.nQ = Loop.n_var_iter + 1 - Par.nM; % changing number of quarterly variables in parameters

        % Adjust transformations (for bridge equations)
        sel_var_m = sel_var(sel_var<=Loop.nM); % monthly variables
        Par.trf.transf_m = Loop.trf.transf_m(sel_var_m);
        sel_var_q = sel_var(sel_var>Loop.nM); % quarterly variables
        sel_var_q = [sel_var_q,(Loop.nM+Loop.nQ)];
        sel_var_q = sel_var_q - Loop.nM;
        Par.trf.transf_q = Loop.trf.transf_q(sel_var_q);

        % Adjust groups of variables
        Loop.groups_sel{n_iter_mod,:} = num2str(groups(sel_var));

        % Preparing name for Excel file
        excel_evalfile = strcat('./eval/',country.name,'/',country.name,'_',country.model,'_evaluation_',Loop.name_loop,'_',num2str(n_iter_mod),'.xlsx'); % Excel file for evaluation metrics


    elseif do_loop == 2 % then we take hyperparameters from pre-defined list

        disp(['LOOP: DOING ITERATION ', num2str(n_iter_mod), ' OUT OF ', num2str(Loop.n_iter)])

        % Load parameters from list of models
        sel_var = sort(str2num(cell2mat(Loop.mdl_list{n_iter_mod,"variables"})));
        do_Covid = Loop.mdl_list{n_iter_mod,"do_Covid"};
        Loop.n_var_iter = Loop.mdl_list{n_iter_mod,"nb_variables"};

        switch country.model
            case 'DFM'
                Loop.parameters(n_iter_mod,:) = Loop.mdl_list{n_iter_mod,2:8};
                if Loop.mdl_list{n_iter_mod,"blocks"} == 1 % Then we add block factors
                    % if block factors are included use:
                    Par.blocks = Loop.blocks;
                    Par.r = ones(size(Loop.blocks,2),1); % Number of block factors
                    Par.r(1) = Loop.mdl_list{n_iter_mod,"nb_factors"};
                else
                    Par.blocks = ones(Par.size_data_in,1); % Then we continue without block factors
                    Par.r = ones(1,1); 
                    Par.r(1) = Loop.mdl_list{n_iter_mod,"nb_factors"};
                end
                Par.p = Loop.mdl_list{n_iter_mod,"nb_lags"};
            case 'BEQ'
                Loop.parameters(n_iter_mod,:) = Loop.mdl_list{n_iter_mod,2:8};
                Par.lagM = Loop.mdl_list{n_iter_mod,"nb_lagM"};
                Par.lagQ = Loop.mdl_list{n_iter_mod,"nb_lagQ"};
                Par.lagY = Loop.mdl_list{n_iter_mod,"nb_lagY"};
                Par.blocks = Loop.blocks; % still put blocks and r to not cause code to stop
                Par.r = ones(size(Loop.blocks,2),1);
                Par.Dum = Loop.Dum;
            case 'BVAR'
                Loop.parameters(n_iter_mod,:) = Loop.mdl_list{n_iter_mod,2:6};
                Par.bvar_lags = Loop.mdl_list{n_iter_mod,"nb_lags"};
                Par.blocks = Loop.blocks; % still put blocks and r to not cause code to stop
                Par.r = ones(size(Loop.blocks,2),1);
        end

        % Adjust the sample dates with new startyear
        startyear = Loop.mdl_list{n_iter_mod,"startyear"};
        startmonth = Loop.startmonth;
        idx_keep = find(Loop.t_m(:,1)==startyear & Loop.t_m(:,2)==startmonth); 
        t_m = Loop.t_m(idx_keep:end,:); % cutting dates
        datet = Loop.datet(idx_keep:end,:); % cutting dates
        xest = Loop.xest(idx_keep:end,:); % cutting the data

        % Adjust the sample variables to sel_var
        Loop.var_sel{n_iter_mod,:} = num2str(sel_var);
        xest = xest(:,[sel_var,end]); % adjusting the sample
                                      % adding end to have the quarterly variable (always last variable)

        % Adjust other parameters
        Par.nM = sum(sel_var<=Loop.nM); % changing number of monthly variables in parameters
                                        % all variables that are before Loop.nM
        Par.nQ = Loop.n_var_iter + 1 - Par.nM; % changing number of quarterly variables in parameters

        
        % Adjust groups of variables
        Loop.groups_sel{n_iter_mod,:} = num2str(groups(sel_var));

        % Adjust transformations (for bridge equations)
        sel_var_m = sel_var(sel_var<=Loop.nM); % monthly variables
        Par.trf.transf_m = Loop.trf.transf_m(sel_var_m);
        sel_var_q = sel_var(sel_var>Loop.nM); % quarterly variables
        sel_var_q = [sel_var_q,(Loop.nM+Loop.nQ)];
        sel_var_q = sel_var_q - Loop.nM;
        Par.trf.transf_q = Loop.trf.transf_q(sel_var_q);

        % Preparing name for Excel file
        excel_evalfile = strcat('./eval/',country.name,'/',country.name,'_',country.model,'_evaluation_',Loop.name_loop,'_',num2str(n_iter_mod),'.xlsx'); % Excel file for evaluation metrics

    end

    % Initial and end date of the out-of-sample
    start_eval = find(t_m(:,1)==Eval.eval_startyear & t_m(:,2)==Eval.eval_startmonth); 
    end_eval = find(t_m(:,1)==Eval.eval_endyear & t_m(:,2)==Eval.eval_endmonth);

    % Structure of missing values (NaN)
    xact = xest(1:end-m,:); % removing last m NaN
    [nb_obs,nb_var] = size(xact);
    nb_var = nb_var - 1; % to account for the fact that the last series is GDP
    Eval.delay = zeros(1,nb_var);
    for nn = 1:nb_var
        last_non_nan = find(~isnan(xact(:,nn)), 1, 'last');
        Eval.delay(1,nn) = nb_obs - last_non_nan;
    end
    
    % Save initial parameters
    % NB: Saved here because can change dynamically in common_NaN_Covid_correct
    nM_init = Par.nM;
    nQ_init = Par.nQ;
    blocks_init = Par.blocks;
    r_init = Par.r;
    transf_m_init = Par.trf.transf_m;
    transf_q_init = Par.trf.transf_q;

    % Initiate containers for predictions and directions
    for hh = 1:numel(C)
        hor = C{hh};
        Eval.(hor).Pred = [];
        Eval.(hor).Dirc = [];
    end

    % Check 1 - do_Covid and dummies for bridge equations
    % If do_Covid = 1 (= dummy correction for Covid), put dummies in Par.Dum (because otherwise dummies will not be selected in combinations)
    if strcmp(country.model,'BEQ')
        if do_Covid == 1
    
            % Changing the value of do_Covid (to avoid doing the corrections below)
            do_Covid = 0; 
    
            % Adding a dummy for 2020Q2 (if not already the case)
            if isempty(Par.Dum)
                Par.Dum = [2020, 6];
            else
                if isempty(find(Par.Dum(:,1)==2020 & Par.Dum(:,2)==6, 1))
                    Par.Dum = [Par.Dum; 2020, 6];
                end
            end
    
            % Adding a dummy for 2020Q3 (if not already the case)
            if isempty(find(Par.Dum(:,1)==2020 & Par.Dum(:,2)==9, 1))
                Par.Dum = [Par.Dum; 2020, 9];
            end
    
        elseif do_Covid ==4
    
            % Changing the value of do_Covid (to avoid doing the corrections below)
            do_Covid = 0; 
    
            % Adding a dummy for 2020Q1 (if not already the case)
            if isempty(Par.Dum)
                Par.Dum = [2020, 3];
            else
                if isempty(find(Par.Dum(:,1)==2020 & Par.Dum(:,2)==3, 1))
                    Par.Dum = [Par.Dum; 2020, 3];
                end
            end
    
            % Adding a dummy for 2020Q2 (if not already the case)
            if isempty(find(Par.Dum(:,1)==2020 & Par.Dum(:,2)==6, 1))
                Par.Dum = [Par.Dum; 2020, 6];
            end    
        
        end
    end

   
    % ---------------------------------------------------------------------
    %% B. LOOP OVER OBSERVATIONS
    % ---------------------------------------------------------------------

    % Loops run to estimate model at each month 
    for i=start_eval:end_eval

        disp(['Out-of-sample predictions for year ', num2str(t_m(i,1)), ' month ', num2str(t_m(i,2))])

        % Re-creates data that would have been available to a real-time forecaster (with NaN)
        xin = xest(1:i,:); % data until month i
        datet_in = datet(1:(i+m),:);
        for j = 1:nM_init
            if Eval.delay(j) ~= 0
                xin(end-Eval.delay(j)+1:end,j) = NaN; % delete last k points to reproduce missing observations
            end
        end
        for j = 1:nQ_init
            if j == nQ_init
                xin(end-Eval.gdp_rel+1:end,nM_init+j) = NaN; % delete last observations 
                                                             % done for GDP series
            elseif Eval.delay(nM_init+j) > 0
                xin(end-Eval.gdp_rel+1:end,nM_init+j) = NaN; % also suppose that the delay is the same as for GDP
                                                             % i.e. this is released in Eval.gdp_rel month of next quarter
            else % in this case it means Eval.delay(Par.nM+j) == 0, i.e. the varibale is available by end of quarter
                 % it can be the case for some US data
                 % in this case, we suppose that data becomes available on last month
                 % meaning - do nothing
            end
        end
        yest = [xin; nan(m,size(xin,2))]; % adding the m NaN at the end for prediction

        % Compute number of months till end of quarter 
        % Should be i+2 in first month, i+1 in second month, and i+0 in third month
        iQ = i + mod(3-mod(t_m(i,2),3),3);

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
                re_estimate = 1;
                coeffs_in = [];
                nameseries = repmat({'name'},1,size(yest_out,2));
                Res = BEQ_estimate(yest_out,Par,datet_in,nameseries,re_estimate,coeffs_in);
            case 'BVAR'
                Res = BVAR_estimate(yest_out,Par,datet_in);
        end

        % Run an AR benchmark
        temp = isfinite(xin(:,end)); % takes automatically last gdp release. 
                                     % usually in the last month gdp of previous quarter is always available
        DQ = xin(temp,end);
        reg = [ones(size(DQ(1:end-1,:),1),1) DQ(1:end-1,:)]; % preparing reg for AR(1)
        beta = inv(reg'*reg)*reg'*DQ(2:end,:); % AR(1) coeff
        num_months = mod(t_m(i,2),3); % number of months in quarter
        if num_months==0 
            num_months=3; % replacing 0 by 3 to have the number of months in quarter
        end
        if num_months < Eval.gdp_rel
            AR_back = [1 DQ(end,:)]*beta;
            AR_back_prev = [DQ(end,:)];
        else
            AR_back = [DQ(end,:)];
            AR_back_prev = [DQ(end-1,:)];
        end
        AR_now = [1 AR_back]*beta;
        AR_1 = [1 AR_now]*beta;

        % Write results for point forecasts
        % 1st column = year of real-time
        % 2nd column = month of real-time
        % 3rd column = year of back/now/fore-cast
        % 4th column = month of back/now/fore-cast
        % 5th column = AR prediction
        % 6th column = DFM prediction
        % 7th column = actual value
        Eval.Bac.Pred = [Eval.Bac.Pred; t_m(i,1:2) t_m(iQ-3,1:2) AR_back Res.X_sm(iQ-3,end) xest(iQ-3,end)];
        if BEQ_monOfQ(t_m(i,2))>= Eval.gdp_rel % manual fix for backcasting if actual data is already there (the reason is that it can be deleted with do_Covid == 2 or 3)
            Eval.Bac.Pred(end,6) = xest(iQ-3,end);
        end
        Eval.Now.Pred = [Eval.Now.Pred; t_m(i,1:2) t_m(iQ,1:2) AR_now Res.X_sm(iQ,end) xest(iQ,end)];
        Eval.For.Pred = [Eval.For.Pred; t_m(i,1:2) t_m(iQ+3,1:2) AR_1 Res.X_sm(iQ+3,end) xest(iQ+3,end)];

        % Write results for directional accuracy
        % 1st column = year of real-time
        % 2nd column = month of real-time
        % 3rd column = year of back/now/fore-cast
        % 4th column = month of back/now/fore-cast
        % 5th column = AR direction (direction means 1 = increase, 0 = decrease)
        % 6th column = DFM direction
        % 7th column = actual direction
        Eval.Bac.Dirc = [Eval.Bac.Dirc; t_m(i,1:2) t_m(iQ-3,1:2) AR_back>AR_back_prev Res.X_sm(iQ-3,end)>Res.X_sm(iQ-6,end) xest(iQ-3,end)>xest(iQ-6,end)];
        if BEQ_monOfQ(t_m(i,2))>= Eval.gdp_rel % manual fix for backcasting if actual data is already there (the reason is that it can be deleted with do_Covid == 2 or 3)
            Eval.Bac.Dirc(end,6) = xest(iQ-3,end)>xest(iQ-6,end);
        end
        Eval.Now.Dirc = [Eval.Now.Dirc; t_m(i,1:2) t_m(iQ,1:2) AR_now>AR_back Res.X_sm(iQ,end)>Res.X_sm(iQ-3,end) xest(iQ,end)>xest(iQ-3,end)];
        Eval.For.Dirc = [Eval.For.Dirc; t_m(i,1:2) t_m(iQ+3,1:2) AR_1>AR_now Res.X_sm(iQ+3,end)>Res.X_sm(iQ,end) xest(iQ+3,end)>xest(iQ,end)];

    end % end of loop on observations


    % ---------------------------------------------------------------------
    %% C. WRITE RESULTS IN EXCEL (FOR ONE MODEL)
    % ---------------------------------------------------------------------

    % Delete past Excel file
    delete(excel_evalfile);

    for hh = 1:numel(C)

        hor = C{hh};

        % Write raw predictions in Excel
        writecell({'year RT','month RT','year','month','AR',country.model,'true'}, ...
                  excel_evalfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  'A1') % names of the columns
        writematrix(Eval.(hor).Pred, ...
                    excel_evalfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    'A2') % actual predictions

        % Evaluate the RMSFE
        Eval.(hor).RMSFE_AR1_all = sqrt(mean((Eval.(hor).Pred(:,5) - Eval.(hor).Pred(:,7)).^2));
        Eval.(hor).RMSFE_mod_all = sqrt(mean((Eval.(hor).Pred(:,6) - Eval.(hor).Pred(:,7)).^2));
        temp_postCovid = Eval.(hor).Pred(:,3)>2020;
        Eval.(hor).RMSFE_AR1_postCovid = sqrt(mean((Eval.(hor).Pred(temp_postCovid,5) - Eval.(hor).Pred(temp_postCovid,7)).^2));
        Eval.(hor).RMSFE_mod_postCovid = sqrt(mean((Eval.(hor).Pred(temp_postCovid,6) - Eval.(hor).Pred(temp_postCovid,7)).^2));  
        temp_preCovid = Eval.(hor).Pred(:,3)<2020;
        Eval.(hor).RMSFE_AR1_preCovid = sqrt(mean((Eval.(hor).Pred(temp_preCovid,5) - Eval.(hor).Pred(temp_preCovid,7)).^2));
        Eval.(hor).RMSFE_mod_preCovid = sqrt(mean((Eval.(hor).Pred(temp_preCovid,6) - Eval.(hor).Pred(temp_preCovid,7)).^2));
        temp_Covid = Eval.(hor).Pred(:,3)==2020;
        Eval.(hor).RMSFE_AR1_Covid = sqrt(mean((Eval.(hor).Pred(temp_Covid,5) - Eval.(hor).Pred(temp_Covid,7)).^2));
        Eval.(hor).RMSFE_mod_Covid = sqrt(mean((Eval.(hor).Pred(temp_Covid,6) - Eval.(hor).Pred(temp_Covid,7)).^2));
        temp_noCovid = Eval.(hor).Pred(:,3)~=2020;
        Eval.(hor).RMSFE_AR1_noCovid = sqrt(mean((Eval.(hor).Pred(temp_noCovid,5) - Eval.(hor).Pred(temp_noCovid,7)).^2));
        Eval.(hor).RMSFE_mod_noCovid = sqrt(mean((Eval.(hor).Pred(temp_noCovid,6) - Eval.(hor).Pred(temp_noCovid,7)).^2));


        % Evaluate the RMSFE at different months
        D = {'first_month','second_month','third_month'};
        for dd = 1:numel(D)
            mth = D{dd};
            temp_D = Eval.(hor).Pred(mod(Eval.(hor).Pred(:,2),3)==mod(dd,3),:); % Select only values corresponding to 1st, 2nd, or 3rd month

            Eval.(hor).(mth).RMSFE_AR1_all = sqrt(mean((temp_D(:,5) - temp_D(:,7)).^2));
            Eval.(hor).(mth).RMSFE_mod_all = sqrt(mean((temp_D(:,6) - temp_D(:,7)).^2));
            
            temp_postCovid = temp_D(:,3)>2020;
            Eval.(hor).(mth).RMSFE_AR1_postCovid = sqrt(mean((temp_D(temp_postCovid,5) - temp_D(temp_postCovid,7)).^2));
            Eval.(hor).(mth).RMSFE_mod_postCovid = sqrt(mean((temp_D(temp_postCovid,6) - temp_D(temp_postCovid,7)).^2));  
            
            temp_preCovid = temp_D(:,3)<2020;
            Eval.(hor).(mth).RMSFE_AR1_preCovid = sqrt(mean((temp_D(temp_preCovid,5) - temp_D(temp_preCovid,7)).^2));
            Eval.(hor).(mth).RMSFE_mod_preCovid = sqrt(mean((temp_D(temp_preCovid,6) - temp_D(temp_preCovid,7)).^2));
            
            temp_Covid = temp_D(:,3)==2020;
            Eval.(hor).(mth).RMSFE_AR1_Covid = sqrt(mean((temp_D(temp_Covid,5) - temp_D(temp_Covid,7)).^2));
            Eval.(hor).(mth).RMSFE_mod_Covid = sqrt(mean((temp_D(temp_Covid,6) - temp_D(temp_Covid,7)).^2));

            temp_noCovid = temp_D(:,3)~=2020;
            Eval.(hor).(mth).RMSFE_AR1_noCovid = sqrt(mean((temp_D(temp_noCovid,5) - temp_D(temp_noCovid,7)).^2));
            Eval.(hor).(mth).RMSFE_mod_noCovid = sqrt(mean((temp_D(temp_noCovid,6) - temp_D(temp_noCovid,7)).^2));
        end

        % Write accuracy metrics (RMSFE) in Excel
        writematrix([Eval.(hor).RMSFE_AR1_all,Eval.(hor).RMSFE_mod_all; ...
                     Eval.(hor).RMSFE_AR1_preCovid,Eval.(hor).RMSFE_mod_preCovid; ...
                     Eval.(hor).RMSFE_AR1_Covid,Eval.(hor).RMSFE_mod_Covid; ...
                     Eval.(hor).RMSFE_AR1_postCovid,Eval.(hor).RMSFE_mod_postCovid; ...
                     Eval.(hor).RMSFE_AR1_noCovid,Eval.(hor).RMSFE_mod_noCovid; ...
                     Eval.(hor).first_month.RMSFE_AR1_all,Eval.(hor).first_month.RMSFE_mod_all; ...
                     Eval.(hor).first_month.RMSFE_AR1_preCovid,Eval.(hor).first_month.RMSFE_mod_preCovid; ...
                     Eval.(hor).first_month.RMSFE_AR1_Covid,Eval.(hor).first_month.RMSFE_mod_Covid; ...
                     Eval.(hor).first_month.RMSFE_AR1_postCovid,Eval.(hor).first_month.RMSFE_mod_postCovid; ...
                     Eval.(hor).first_month.RMSFE_AR1_noCovid,Eval.(hor).first_month.RMSFE_mod_noCovid; ...
                     Eval.(hor).second_month.RMSFE_AR1_all,Eval.(hor).second_month.RMSFE_mod_all; ...
                     Eval.(hor).second_month.RMSFE_AR1_preCovid,Eval.(hor).second_month.RMSFE_mod_preCovid; ...
                     Eval.(hor).second_month.RMSFE_AR1_Covid,Eval.(hor).second_month.RMSFE_mod_Covid; ...
                     Eval.(hor).second_month.RMSFE_AR1_postCovid,Eval.(hor).second_month.RMSFE_mod_postCovid; ...
                     Eval.(hor).second_month.RMSFE_AR1_noCovid,Eval.(hor).second_month.RMSFE_mod_noCovid; ...
                     Eval.(hor).third_month.RMSFE_AR1_all,Eval.(hor).third_month.RMSFE_mod_all; ...
                     Eval.(hor).third_month.RMSFE_AR1_preCovid,Eval.(hor).third_month.RMSFE_mod_preCovid; ...
                     Eval.(hor).third_month.RMSFE_AR1_Covid,Eval.(hor).third_month.RMSFE_mod_Covid; ...
                     Eval.(hor).third_month.RMSFE_AR1_postCovid,Eval.(hor).third_month.RMSFE_mod_postCovid; ...
                     Eval.(hor).third_month.RMSFE_AR1_noCovid,Eval.(hor).third_month.RMSFE_mod_noCovid], ...
                    excel_evalfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    strcat('E',num2str(size(Eval.(hor).Pred,1)+3))) % RMSE

        % Write lines names in Excel
        writecell({'All';'pre-Covid (<2020)';'Covid (2020)';'post-Covid (>2020)'; 'no Covid'; ...
                   '1m - all';'1m - pre-Covid (<2020)';'1m - Covid (2020)';'1m - post-Covid (>2020)'; '1m - no Covid'; ...
                   '2m - all';'2m - pre-Covid (<2020)';'2m - Covid (2020)';'2m - post-Covid (>2020)'; '2m - no Covid'; ...
                   '3m - all';'3m - pre-Covid (<2020)';'3m - Covid (2020)';'3m - post-Covid (>2020)'; '3m - no Covid'; ...
                   'NB: Refers to year of prediction (for example, a one-quarter-ahead forecast made in Q4 2019 will enter the Covid (2020) period)'}, ...
                  excel_evalfile, ...
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  strcat('A',num2str(size(Eval.(hor).Pred,1)+3))) % Names

        % Write raw directions in Excel
        writecell({'year RT','month RT','year','month','AR',country.model,'true'}, ...
                  excel_evalfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  'I1') % names of the columns
        writematrix(Eval.(hor).Dirc, ...
                    excel_evalfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    'I2') % actual predictions
        
        % Evaluate the forecast directional accuracy (FDA)
        Eval.(hor).FDA_AR1_all = 1-mean(abs(Eval.(hor).Dirc(:,5)-Eval.(hor).Dirc(:,7)));
        Eval.(hor).FDA_mod_all = 1-mean(abs(Eval.(hor).Dirc(:,6)-Eval.(hor).Dirc(:,7)));
        temp_postCovid = Eval.(hor).Dirc(:,3)>2020;
        Eval.(hor).FDA_AR1_postCovid = 1-mean(abs(Eval.(hor).Dirc(temp_postCovid,5)-Eval.(hor).Dirc(temp_postCovid,7)));
        Eval.(hor).FDA_mod_postCovid = 1-mean(abs(Eval.(hor).Dirc(temp_postCovid,6)-Eval.(hor).Dirc(temp_postCovid,7)));
        temp_preCovid = Eval.(hor).Dirc(:,3)<2020;
        Eval.(hor).FDA_AR1_preCovid = 1-mean(abs(Eval.(hor).Dirc(temp_preCovid,5)-Eval.(hor).Dirc(temp_preCovid,7)));
        Eval.(hor).FDA_mod_preCovid = 1-mean(abs(Eval.(hor).Dirc(temp_preCovid,6)-Eval.(hor).Dirc(temp_preCovid,7)));
        temp_Covid = Eval.(hor).Dirc(:,3)==2020;
        Eval.(hor).FDA_AR1_Covid = 1-mean(abs(Eval.(hor).Dirc(temp_Covid,5)-Eval.(hor).Dirc(temp_Covid,7)));
        Eval.(hor).FDA_mod_Covid = 1-mean(abs(Eval.(hor).Dirc(temp_Covid,6)-Eval.(hor).Dirc(temp_Covid,7)));
        temp_noCovid = Eval.(hor).Dirc(:,3)~=2020;
        Eval.(hor).FDA_AR1_noCovid = 1-mean(abs(Eval.(hor).Dirc(temp_noCovid,5)-Eval.(hor).Dirc(temp_noCovid,7)));
        Eval.(hor).FDA_mod_noCovid = 1-mean(abs(Eval.(hor).Dirc(temp_noCovid,6)-Eval.(hor).Dirc(temp_noCovid,7)));

        % Evaluate the FDA at different months
        for dd = 1:numel(D)
            mth = D{dd};
            temp_D = Eval.(hor).Dirc(mod(Eval.(hor).Dirc(:,2),3)==mod(dd,3),:); % Select only values corresponding to 1st, 2nd, or 3rd month

            Eval.(hor).(mth).FDA_AR1_all = 1-mean(abs(temp_D(:,5)-temp_D(:,7)));
            Eval.(hor).(mth).FDA_mod_all = 1-mean(abs(temp_D(:,6)-temp_D(:,7)));
            
            temp_postCovid = temp_D(:,3)>2020;
            Eval.(hor).(mth).FDA_AR1_postCovid = 1-mean(abs(temp_D(temp_postCovid,5)-temp_D(temp_postCovid,7)));
            Eval.(hor).(mth).FDA_mod_postCovid = 1-mean(abs(temp_D(temp_postCovid,6)-temp_D(temp_postCovid,7))); 
            
            temp_preCovid = temp_D(:,3)<2020;
            Eval.(hor).(mth).FDA_AR1_preCovid = 1-mean(abs(temp_D(temp_preCovid,5)-temp_D(temp_preCovid,7)));
            Eval.(hor).(mth).FDA_mod_preCovid = 1-mean(abs(temp_D(temp_preCovid,6)-temp_D(temp_preCovid,7)));
            
            temp_Covid = temp_D(:,3)==2020;
            Eval.(hor).(mth).FDA_AR1_Covid = 1-mean(abs(temp_D(temp_Covid,5)-temp_D(temp_Covid,7)));
            Eval.(hor).(mth).FDA_mod_Covid = 1-mean(abs(temp_D(temp_Covid,6)-temp_D(temp_Covid,7)));

            temp_noCovid = temp_D(:,3)~=2020;
            Eval.(hor).(mth).FDA_AR1_noCovid = 1-mean(abs(temp_D(temp_noCovid,5)-temp_D(temp_noCovid,7)));
            Eval.(hor).(mth).FDA_mod_noCovid = 1-mean(abs(temp_D(temp_noCovid,6)-temp_D(temp_noCovid,7))); 
        end

        % Write accuracy metrics (FDA) in Excel
        writematrix([Eval.(hor).FDA_AR1_all,Eval.(hor).FDA_mod_all; ...
                     Eval.(hor).FDA_AR1_preCovid,Eval.(hor).FDA_mod_preCovid; ...
                     Eval.(hor).FDA_AR1_Covid,Eval.(hor).FDA_mod_Covid; ...
                     Eval.(hor).FDA_AR1_postCovid,Eval.(hor).FDA_mod_postCovid; ...
                     Eval.(hor).FDA_AR1_noCovid,Eval.(hor).FDA_mod_noCovid; ...
                     Eval.(hor).first_month.FDA_AR1_all,Eval.(hor).first_month.FDA_mod_all; ...
                     Eval.(hor).first_month.FDA_AR1_preCovid,Eval.(hor).first_month.FDA_mod_preCovid; ...
                     Eval.(hor).first_month.FDA_AR1_Covid,Eval.(hor).first_month.FDA_mod_Covid; ...
                     Eval.(hor).first_month.FDA_AR1_postCovid,Eval.(hor).first_month.FDA_mod_postCovid; ...
                     Eval.(hor).first_month.FDA_AR1_noCovid,Eval.(hor).first_month.FDA_mod_noCovid; ...
                     Eval.(hor).second_month.FDA_AR1_all,Eval.(hor).second_month.FDA_mod_all; ...
                     Eval.(hor).second_month.FDA_AR1_preCovid,Eval.(hor).second_month.FDA_mod_preCovid; ...
                     Eval.(hor).second_month.FDA_AR1_Covid,Eval.(hor).second_month.FDA_mod_Covid; ...
                     Eval.(hor).second_month.FDA_AR1_postCovid,Eval.(hor).second_month.FDA_mod_postCovid; ...
                     Eval.(hor).second_month.FDA_AR1_noCovid,Eval.(hor).second_month.FDA_mod_noCovid; ...
                     Eval.(hor).third_month.FDA_AR1_all,Eval.(hor).third_month.FDA_mod_all; ...
                     Eval.(hor).third_month.FDA_AR1_preCovid,Eval.(hor).third_month.FDA_mod_preCovid; ...
                     Eval.(hor).third_month.FDA_AR1_Covid,Eval.(hor).third_month.FDA_mod_Covid; ...
                     Eval.(hor).third_month.FDA_AR1_postCovid,Eval.(hor).third_month.FDA_mod_postCovid; ...
                     Eval.(hor).third_month.FDA_AR1_noCovid,Eval.(hor).third_month.FDA_mod_noCovid], ...
                    excel_evalfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    strcat('M',num2str(size(Eval.(hor).Dirc,1)+3))) % FDA


        % ---------------------------------------------------------------------
        %% D. WRITE RESULTS IN LOOP.RESULTS (IF DO_LOOP==1 or 2)
        % ---------------------------------------------------------------------

        % Write in loop (if relevant)
        if do_loop == 1 || do_loop == 2

            % Write accuracy metrics (RMSE and FDA)
            Loop.(hor).results(n_iter_mod,1) = Eval.(hor).RMSFE_mod_all;
            Loop.(hor).results(n_iter_mod,2) = Eval.(hor).RMSFE_mod_preCovid;
            Loop.(hor).results(n_iter_mod,3) = Eval.(hor).RMSFE_mod_Covid;
            Loop.(hor).results(n_iter_mod,4) = Eval.(hor).RMSFE_mod_postCovid;
            Loop.(hor).results(n_iter_mod,5) = Eval.(hor).RMSFE_mod_noCovid;
            Loop.(hor).results(n_iter_mod,6) = Eval.(hor).FDA_mod_all;
            Loop.(hor).results(n_iter_mod,7) = Eval.(hor).FDA_mod_preCovid;
            Loop.(hor).results(n_iter_mod,8) = Eval.(hor).FDA_mod_Covid;
            Loop.(hor).results(n_iter_mod,9) = Eval.(hor).FDA_mod_postCovid;
            Loop.(hor).results(n_iter_mod,10) = Eval.(hor).FDA_mod_noCovid;
            Loop.(hor).results(n_iter_mod,11) = Eval.(hor).first_month.RMSFE_mod_all;
            Loop.(hor).results(n_iter_mod,12) = Eval.(hor).first_month.RMSFE_mod_preCovid;
            Loop.(hor).results(n_iter_mod,13) = Eval.(hor).first_month.RMSFE_mod_Covid;
            Loop.(hor).results(n_iter_mod,14) = Eval.(hor).first_month.RMSFE_mod_postCovid;
            Loop.(hor).results(n_iter_mod,15) = Eval.(hor).first_month.RMSFE_mod_noCovid;
            Loop.(hor).results(n_iter_mod,16) = Eval.(hor).first_month.FDA_mod_all;
            Loop.(hor).results(n_iter_mod,17) = Eval.(hor).first_month.FDA_mod_preCovid;
            Loop.(hor).results(n_iter_mod,18) = Eval.(hor).first_month.FDA_mod_Covid;
            Loop.(hor).results(n_iter_mod,19) = Eval.(hor).first_month.FDA_mod_postCovid;
            Loop.(hor).results(n_iter_mod,20) = Eval.(hor).first_month.FDA_mod_noCovid;
            Loop.(hor).results(n_iter_mod,21) = Eval.(hor).second_month.RMSFE_mod_all;
            Loop.(hor).results(n_iter_mod,22) = Eval.(hor).second_month.RMSFE_mod_preCovid;
            Loop.(hor).results(n_iter_mod,23) = Eval.(hor).second_month.RMSFE_mod_Covid;
            Loop.(hor).results(n_iter_mod,24) = Eval.(hor).second_month.RMSFE_mod_postCovid;
            Loop.(hor).results(n_iter_mod,25) = Eval.(hor).second_month.RMSFE_mod_noCovid;
            Loop.(hor).results(n_iter_mod,26) = Eval.(hor).second_month.FDA_mod_all;
            Loop.(hor).results(n_iter_mod,27) = Eval.(hor).second_month.FDA_mod_preCovid;
            Loop.(hor).results(n_iter_mod,28) = Eval.(hor).second_month.FDA_mod_Covid;
            Loop.(hor).results(n_iter_mod,29) = Eval.(hor).second_month.FDA_mod_postCovid;
            Loop.(hor).results(n_iter_mod,30) = Eval.(hor).second_month.FDA_mod_noCovid;
            Loop.(hor).results(n_iter_mod,31) = Eval.(hor).third_month.RMSFE_mod_all;
            Loop.(hor).results(n_iter_mod,32) = Eval.(hor).third_month.RMSFE_mod_preCovid;
            Loop.(hor).results(n_iter_mod,33) = Eval.(hor).third_month.RMSFE_mod_Covid;
            Loop.(hor).results(n_iter_mod,34) = Eval.(hor).third_month.RMSFE_mod_postCovid;
            Loop.(hor).results(n_iter_mod,35) = Eval.(hor).third_month.RMSFE_mod_noCovid;
            Loop.(hor).results(n_iter_mod,36) = Eval.(hor).third_month.FDA_mod_all;
            Loop.(hor).results(n_iter_mod,37) = Eval.(hor).third_month.FDA_mod_preCovid;
            Loop.(hor).results(n_iter_mod,38) = Eval.(hor).third_month.FDA_mod_Covid;
            Loop.(hor).results(n_iter_mod,39) = Eval.(hor).third_month.FDA_mod_postCovid;
            Loop.(hor).results(n_iter_mod,40) = Eval.(hor).third_month.FDA_mod_noCovid;

            % Write parameters and variables
            writematrix(Loop.parameters(n_iter_mod,:), ...
                        excel_evalfile, ...
                        'Sheet', ...
                        'Parameters', ...
                        'Range', ...
                        'A2') % Parameters
            str_col = common_num2xlcol(1 + size(Loop.parameters,2));
            col_var = strcat(str_col,num2str(2));
            writecell(Loop.var_sel(n_iter_mod,:), ...
                      excel_evalfile, ...
                      'Sheet', ...
                      'Parameters', ...
                      'Range', ...
                      col_var) % Variables selected
            str_col = common_num2xlcol(1 + size(Loop.parameters,2) + 1);
            col_var = strcat(str_col,num2str(2));
            writecell(Loop.groups_sel(n_iter_mod,:), ...
                      excel_evalfile, ...
                      'Sheet', ...
                      'Parameters', ...
                      'Range', ...
                      col_var) % Groups of variables
            writecell(Loop.parameters_name, ...
                      excel_evalfile, ...
                      'Sheet', ...
                      'Parameters', ...
                      'Range', ...
                      'A1') % Names of variables


        end

    end % end of loop on horizons

end % end of loop on iterations


% ---------------------------------------------------------------------
%% 3. WRITE RESULTS ACROSS ALL MODELS (IF DO_LOOP==1 or 2)
% ---------------------------------------------------------------------

if do_loop == 1 || do_loop == 2

    % Delete past Excel file
    delete(Loop.excel_loopfile);

    for hh = 1:numel(C)

        hor = C{hh};
        
        % Write loop name
        Loop_ID = cell(Loop.n_iter,1);
        if do_loop == 1
            Loop_ID(:) = {Loop.name_loop};
        elseif do_loop == 2 %then we write batch names from model list
            Loop_ID(:) = Loop.mdl_list{:,"batch"};
        end
        writecell(Loop_ID, ...
                  Loop.excel_loopfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  'A2') % names of the loop              

        % Write results of the loop
        writecell(Loop.results_names, ...
                  Loop.excel_loopfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  'A1') % names of the columns
        writematrix(Loop.parameters, ...
                    Loop.excel_loopfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    'B2') % metrics
        str_col = common_num2xlcol(1 + size(Loop.parameters,2) + 1);
        col_var = strcat(str_col,num2str(2));
        writematrix(Loop.(hor).results, ...
                    Loop.excel_loopfile, ...
                    'Sheet', ...
                    hor, ...
                    'Range', ...
                    col_var) % metrics
        str_col = common_num2xlcol(1 + size(Loop.parameters,2) + 1 + size(Loop.(hor).results,2));
        col_var = strcat(str_col,num2str(2));
        writecell(Loop.var_sel, ...
                  Loop.excel_loopfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  col_var) % variables selected
        str_col = common_num2xlcol(1 + size(Loop.parameters,2) + 1 + size(Loop.(hor).results,2) + 1);
        col_var = strcat(str_col,num2str(2));
        writecell(Loop.groups_sel, ...
                  Loop.excel_loopfile, ...                  
                  'Sheet', ...
                  hor, ...
                  'Range', ...
                  col_var) % groups selected

    end % end of loop on horizons

end % end of condition on do_loop

end % end of function