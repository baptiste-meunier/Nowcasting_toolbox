%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%                            NOWCASTING TOOLBOX                            %
%                                                                          %
%    When using models in the toolbox, users are kindly requested          %
%    to cite the original papers:                                          %
%   - For DFM: Bańbura, M., and Modugno, M. (2014). "Maximum likelihood    %
%              estimation of factor models on datasets with arbitrary      %
%              pattern of missing data", Journal of Applied Econometrics,  % 
%              29(11), 133–160                                             %
%   - For DFM and block structure (on top of above): Delle Chiaie, S.,     %
%             Ferrara, L., and Giannone, D. (2022). "Common factors of     %
%             commodity prices", Journal of Applied Econometrics, 37(3),   % 
%             461–476                                                      %
%   - For BVAR: Cimadomo, J., Giannone, D., Lenza, M., Monti, F., and      %
%               Sokol, A. (2022). "Nowcasting with large Bayesian vector   %
%               autoregressions", Journal of Econometrics, 231(2), 500–519 %
%   - For bridge equations: Bańbura, M., Belousova, I., Bodnár, K., and    %
%                           Tóth, M. B. (2023). "Nowcasting employment in  %
%                           the euro area", Working Paper Series, No 2815, % 
%                           European Central Bank                          %
%                                                                          %
%    Approximate running times in nowcasting mode (do_eval = 0 and         %
%    do_loop = 0) for DFM (NB: approximate timing, depends on model        %
%    specifications and machine running the code                           %
%   - Only nowcast       (do_range = 0 / do_mae = 0) --> around 1 minute   %
%   - Nowcast with range (do_range = 1 / do_mae = 0) --> around 15 minutes %
%   - Nowcast with MAE   (do_range = 0 / do_mae = 1) --> around 20 minutes %
%   - Nowcast with both  (do_range = 1 / do_mae = 1) --> around 35 minutes %
%                                                                          %
%    The toolbox builds on long-standing effort in the External            %
%    Developments division of the ECB to develop a nowcasting tool. This   %
%    involved a number of colleagues notably Simona Delle Chiaie, Frederik % 
%    Kurcz, and Gabriel Perez-Quiros who conducted previous nowcasting     %
%    tools in the division, as well as Davide Brignone, Alistair Dieppe,   %
%    Julia Doleschel, Rinalds Gerinovics, and Roberta de Stefani.          %
%                                                                          %
%    These programmes are the responsibilities of the authors and not of   %
%    the ECB and all errors and ommissions remain those of the authors.    %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear


% -------------------------------------------------------------------------
%% 0. TOOLBOX SETTINGS (to be changed by the user)
% -------------------------------------------------------------------------

do_eval = 1;   % switch on the use of the toolbox
               % 0 = nowcast
               % 1 = model evaluation      

do_loop = 0;   % switch on whether to loop over different models
               % 0 = single model (user-defined in code below)
               % 1 = automatic loop over random models (selected within bounds set by user in code below)      
               % 2 = custom loop over user-defined list of models (in Eval_list_mod.xlsx)
               % NB: if do_loop is 1 or 2 then do_eval is automatically set to 1          

do_range = 0;  % switch on whether to run alternative models (obtained by disconnecting 1 or 2 groups of variables)
               % 0 = no
               % 1 = yes (takes significantly longer to estimate)

do_mae = 0;    % switch on how to compute the Mean Absolute Error (MAE) and Forecast Directional Accuracy (FDA) from past forecast errors
               % 0 = take user-specified values (in code below)
               % 1 = compute based on past 10 years (takes significantly longer to estimate)

do_subset = 0; % switch on whether to take a sub-set of the input data
               % 0 = take full input data (in DFM_data_xxx.xlsx)
               % 1 = take only the sample var_keep 


% -------------------------------------------------------------------------
%% 1. MODEL INPUTS (to be changed by the user)
% -------------------------------------------------------------------------

% Country hyper-parameters
country.name = 'Example1';
country.model = 'DFM'; % either 'DFM' or 'BEQ' or 'BVAR'

% Model specifications
if strcmp(country.name,'Example1')

    % common
    Par.startyear = 2005;       % starting year for estimation
    Par.startmonth = 1;         % starting month for estimation
    do_Covid = 0;               % switch on correction for Covid observations
                                % 0 = do nothing
                                % 1 = add dummies (one for June 2020 and one for Sept. 2020)
                                % 2 = put to NaN (from Feb. 2020 to Sep. 2020 included)
                                % 3 = outlier-correction (with outliers replaces by NaN)
                                % 4 = add dummies (one for Mar. 2020 and one for June 2020)

    % DFM
    Par.p = 4;                  % number of lags 
    Par.r = 5;                  % number of factors
    Par.idio = 1;               % idiosyncratic component specification. 0 = iid, 1 = AR(1).
    Par.thresh = 1e-4;          % threshold for convergence of EM algorithm
    Par.max_iter = 100;         % number of iterations of EM algorithm
    Par.block_factors = 0;      % switch to include block factors. NB: if included, will be one factor per block (on top of global factors). Can be changed in code.

    % Bridge equation (BEQ)
    Par.lagM = 1;               % number of lags for monthly regressor(s) (in quarterly terms). 0 = only contemporaneous.
    Par.lagQ = 1;               % number of lags for quarterly regressor(s) (in quarterly terms). 0 = only contemporaneous.
    Par.lagY = 1;               % number of lags for the endogenous variable (in quarterly terms). 0 = no lags of endogenous.
    Par.type = 901;             % type of interpolation (see BEQ_estimate for details)
    Par.Dum = [2020,3; ...
               2020,6;
               2020,9];         % dates of the dummies (year, month). Format should be a k x 2 matrix with k = number of dummies (each dummy on a row) and year / month in columns.

    % BVAR
    Par.bvar_lags = 5;          % number of lags
    Par.bvar_thresh = 1e-6;     % threshold for convergence
    Par.bvar_max_iter = 200;    % number of iterations

    % Additional input for mean absolute error
    MAE.Bac.mae_1st = 0.15;     % MAE, backcasting, 1st month of quarter
    MAE.Bac.mae_2nd = 0.15;     % MAE, backcasting, 2nd month of quarter
    MAE.Bac.mae_3rd = 0.15;     % MAE, backcasting, 3rd month of quarter
    MAE.Now.mae_1st = 0.24;     % MAE, nowcasting, 1st month of quarter
    MAE.Now.mae_2nd = 0.24;     % MAE, nowcasting, 2nd month of quarter
    MAE.Now.mae_3rd = 0.21;     % MAE, nowcasting, 3rd month of quarter
    MAE.For.mae_1st = 0.48;     % MAE, forecasting, 1st month of quarter
    MAE.For.mae_2nd = 0.30;     % MAE, forecasting, 2nd month of quarter
    MAE.For.mae_3rd = 0.29;     % MAE, forecasting, 3rd month of quarter
    MAE.Bac.fda_1st = 0.88;     % FDA, backcasting, 1st month of quarter
    MAE.Bac.fda_2nd = 0.88;     % FDA, backcasting, 2nd month of quarter
    MAE.Bac.fda_3rd = 0.88;     % FDA, backcasting, 3rd month of quarter
    MAE.Now.fda_1st = 0.73;     % FDA, nowcasting, 1st month of quarter
    MAE.Now.fda_2nd = 0.85;     % FDA, nowcasting, 2nd month of quarter
    MAE.Now.fda_3rd = 0.88;     % FDA, nowcasting, 3rd month of quarter
    MAE.For.fda_1st = 0.56;     % FDA, forecasting, 1st month of quarter
    MAE.For.fda_2nd = 0.68;     % FDA, forecasting, 2nd month of quarter
    MAE.For.fda_3rd = 0.70;     % FDA, forecasting, 3rd month of quarter                                                          

end

% Additional inputs for model evaluation
% NB: if do_eval=1
Eval.data_update_lastyear = 2023;      % date at which the data for evaluation was updated
Eval.data_update_lastmonth = 10;
Eval.eval_startyear = 2020;            % start date for out-of-sample evaluation
Eval.eval_startmonth = 10;
Eval.eval_endyear = 2022;              % end date for out-of-sample evaluation
Eval.eval_endmonth = 10;
Eval.gdp_rel = 2;                      % month of quarter at which GDP (for previous quarter) is available. Should be at most 3. 
                                       % for example, 2 means that GDP is available on the second month (but not on the first)

% Additional inputs for loop across models
% NB: if do_loop=1
Loop.n_iter = 10;                      % number of random models to be tested
Loop.name_loop = 'b1';                 % name for the loop
Loop.min_startyear = 2008;             % minimum start year of estimation
Loop.max_startyear = 2014;             % maximum start year of estimation
Loop.startmonth = 1;                   % start month of estimation (same for any start year)
Loop.min_var = 5;                      % minimum number of variables
Loop.max_var = 10;                     % maximum number of variables
Loop.min_p = 1;                        % DFM: minimum number of lags
Loop.max_p = 4;                        % DFM: maximum number of lags (NB: p>5 won't work)
Loop.min_r = 2;                        % DFM: minimum number of factors
Loop.max_r = 6;                        % DFM: maximum number of factors
Loop.min_lagM = 1;                     % BEQ: minimum number of lags (in quarterly terms) for monthly regressors
Loop.max_lagM = 4;                     % BEQ: maximum number of lags (in quarterly terms) for monthly regressors
Loop.min_lagQ = 1;                     % BEQ: minimum number of lags for quarterly regressors
Loop.max_lagQ = 4;                     % BEQ: maximum number of lags for quarterly regressors
Loop.min_lagY = 1;                     % BEQ: minimum number of lags for the endogenous variable
Loop.max_lagY = 2;                     % BEQ: maximum number of lags for the endogenous variable
Loop.min_bvar_lags = 2;                % BVAR: minimum number of lags
Loop.max_bvar_lags = 4;                % BVAR: maximum number of lags
Loop.do_random = 1;                    % switch on whether to randomize models in the loop
                                       % 0 = different runs of loop will select same models (can be useful to test the effect of one change on same set of models, e.g. Covid corrections)
                                       % 1 (default) = different runs of the loop will select different models
                                                            
% Additional inputs for custom loop across models 
% NB: if do_loop=2
Loop.list_name = 'Eval_list_DFM.xlsx'; % Name of excel with list of models (has to be located in folder 'eval/country.name/')
Loop.name_customloop = 'customloop';   % Name of the loop
Loop.alter_covid = 1;                  % Switch to 1 to alter between Covid corrections (0 tests the list of models)

% Additional inputs for subsetting input data
%NB: if do_subset=1
var_keep = [1   3   5   7   9  10  11  13  14  15  19  21  22  24  25  26  27  28  30  31];        % variables to keep


% -------------------------------------------------------------------------
% -------------- BELOW THIS LINE CODE SHOULD NOT BE MODIFIED --------------
% -------------------------------------------------------------------------

disp ('Section 1: Model input loaded');


% -------------------------------------------------------------------------
%% 2. ADD FOLDERS AND FILES - ALSO RUN SOME CHECKS
% -------------------------------------------------------------------------

% Generic parameters
m = 6; % number of months ahead 

% Get (or simulate for evaluation) today's date
if do_eval == 0
    date_today = [year(datetime("today")),month(datetime("today"))]; % today's date
elseif do_eval == 1
    date_today = [Eval.data_update_lastyear,Eval.data_update_lastmonth]; % date of the data "freeze" for out-of-sample evaluation
end

% Add paths
addpath('./tools');
addpath('./dataset');
namesave = strcat('sav_',date); % current date
outputfolder = strcat('./output/',country.name,'/');
rootfolder = cd;

% Update name of loop if do_loop == 2
if do_loop == 2
    Loop.name_loop = Loop.name_customloop;
end

% Prepare file names
excel_datafile = strcat('data_',country.name); % Excel file containing data (if users use exceldata =1)
excel_outputfile = strcat('./output/',country.name,'/',country.name,'_tracking.xlsx'); % Excel file containing tracking and news decomposition
Loop.excel_loopfile = strcat('./eval/',country.name,'/',country.name,'_',country.model,'_loop_',Loop.name_loop,'.xlsx'); % Excel file for loop over random models
newsfile = 'cur_nowcast.mat'; % compare news relative to this run

% Check 1 - Matlab version compatibility
if verLessThan('matlab','9.12')
    error('Running EXT-Now requires Matlab version R2022a or higher.')
end

% Check 2 - assign do_eval=1 if do_loop is selected
if do_loop == 1 || do_loop == 2
    do_eval = 1;
end

% Check 3 - do_Covid and dummies for bridge equations
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

disp ('Section 2: Folders and files added');


% -------------------------------------------------------------------------    
%% 3. LOAD DATA & INITIALIZE INPUT FOR ESTIMATION
% -------------------------------------------------------------------------

% Load data
mon_freq = 'Monthly';
quar_freq = 'Quarterly';
blocks_sheet = 'blocks';
[Par,xest,t_m,groups,nameseries,blocks,groups_name,fullnames,datet,Loop] = ...
        common_load_data(excel_datafile,mon_freq,quar_freq,blocks_sheet,Par,m,do_loop,date_today,Loop);

% Subset data if required by user
if do_subset == 1

    % Saving initial number of variables
    nM_init = Par.nM; 
    nQ_init = Par.nQ; 

    % Adjusting data
    xest = xest(:,[var_keep,end]);       
    Par.nM = sum(var_keep<=nM_init);
    Par.nQ = length(var_keep) - Par.nM + 1;

    % Adjusting other parameters
    Par.blocks = Par.blocks([var_keep,end]',:);
    groups = groups([var_keep,end]);
    fullnames = fullnames([var_keep,end]);
    ID_groups = unique(groups(:));
    groups_name = groups_name(ID_groups);

    % Adjust transformations (for bridge equations)
    var_keep_m = var_keep(var_keep<=nM_init); % monthly variables
    Par.trf.transf_m = Par.trf.transf_m(var_keep_m);
    var_keep_q = var_keep(var_keep>nM_init); % quarterly variables
    var_keep_q = [var_keep_q,(nM_init + nQ_init)];
    var_keep_q = var_keep_q - nM_init;
    Par.trf.transf_q = Par.trf.transf_q(var_keep_q);

end

clear p r exceldata block_factors idiosyncratic

disp ('Section 3: Data loaded');

if do_eval == 0

    % ---------------------------------------------------------------------
    %% 4. ESTIMATE DFM AND COMPUTE NEWS DECOMPOSITION
    % ---------------------------------------------------------------------

    % Compute heatmap of input variables
    % NB: for the heatmap we consider  data before Covid correction
    heatmap = common_heatmap(xest,Par,groups,groups_name,fullnames);
    
    disp ('Section 4: Estimation starts');

    % Save initial number of variables
    Par.nM_init = Par.nM;

    % Correct for Covid period and NaN in the estimation sample
    [xest_out,nM_out,blocks_out,r_out,groups_out,groups_name_out,nameseries_out,fullnames_out,transf_m_out,transf_q_out] = ...
        common_NaN_Covid_correct(xest,datet,do_Covid,Par.nM,Par.blocks,Par.r,groups,groups_name,nameseries,fullnames,Par.trf.transf_m,Par.trf.transf_q);

    % Adjust parameters post-correction
    Par.nM = nM_out;
    Par.nQ = size(xest_out,2) - Par.nM;
    Par.blocks = blocks_out;
    Par.r = r_out;
    Par.trf.transf_m = transf_m_out;
    Par.trf.transf_q = transf_q_out;
    groups = groups_out;
    groups_name = groups_name_out;
    nameseries = nameseries_out;
    fullnames = fullnames_out;

    % Prepare the estimation
    warning('off','all');
    datet_in = datet;

    % Run the estimation
    switch country.model
        case 'DFM'
            Res = DFM_estimate(xest_out,Par);
        case 'BEQ'                    
            Res = BEQ_estimate(xest_out,Par,datet_in,nameseries,1,[]);
        case 'BVAR'
            Res = BVAR_estimate(xest_out,Par,datet);
    end

    % Adjust results
    Res.groups = groups; 
    Res.series = nameseries;
    Res.name_descriptor = fullnames;
    GDP_track = [datet Res.X_sm(:,end)];
    
    disp ('Section 4: Estimation completed');
    
    % ---------------------------------------------------------------------
    %% 5. COMPUTE NEWS DECOMPOSITION
    % ---------------------------------------------------------------------
    
    % Computes news relative to nowcast and forecast in 'newsfile'
    switch country.model
        case 'DFM'
            [news_results,news_results_fcst,table_now,table_fcst,prev_news] = ...
                DFM_News_Mainfile(outputfolder,newsfile,groups_name,xest_out,Res,Par,namesave,datet,nameseries,country.model);
        case 'BEQ'
            [news_results,news_results_fcst,prev_news] = ...
                BEQ_News_Mainfile(outputfolder,newsfile,groups_name,xest_out,Res,Par,namesave,datet,nameseries,country.model);
        case 'BVAR'
            [news_results,news_results_fcst,table_now,table_fcst,prev_news] = ...
                BVAR_News_Mainfile(outputfolder,newsfile,groups_name,xest_out,Res,Par,namesave,datet,nameseries,country.model);
    end

    % ---------------------------------------------------------------------
    %% 6. RANGE OF NOWCASTS
    % ---------------------------------------------------------------------  

    if do_range == 1 % Calculates range of nowcast if required by user
        
        disp('Section 6: Range of nowcasts starts');
        [range] = common_range(xest_out,Par,datet,groups_name,groups,country,nameseries);
        disp('Section 6: Range of nowcasts completed');
    
    elseif do_range == 0

        range = [];
        disp('Section 6: No range of nowcasts computed, as do_range is set to 0');

    else

        error('do_range should be either 0 or 1')

    end
    

    % ---------------------------------------------------------------------
    %% 7. ERROR EVALUATION
    % ---------------------------------------------------------------------

    MAE = common_mae(xest,Par,t_m,m,datet,do_Covid,country,MAE,do_mae,nameseries);

    % ---------------------------------------------------------------------
    %% 8. SAVE RESULTS
    % ---------------------------------------------------------------------
    
    common_save_results;

    disp ('Section 8: Results saved');

elseif do_eval == 1

    % Run the evaluation (the script also produces Excel files)
    warning('off','all');
    [Loop,Eval] = common_eval_models(do_loop,Loop,Eval,xest,Par,t_m,m,country,datet,do_Covid,groups);

end

disp ('End nowcast');

toc