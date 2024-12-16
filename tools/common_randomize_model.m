function [rnd_parameters,sel_var] = common_randomize_model(n_iter_mod,Loop,country_model,do_Covid,Par)
% This script provides random model hyper-parameters 
%
% INPUTS
% - n_iter_mod [scalar] = ID of iteration
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
%       o do_random [scalar 0/1] = switch on whether to randomize the models in the loop
% - country_model [string] = model type
% - do_Covid [scalar] = method for Covid correction
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
%
% OUTPUTS
% - rnd_parameters [vector] = vector of (randomized) hyper-parameters
% - sel_var [vector] = ID of variables ot be kept in the sample
% - n_col_results [scalar] = number of columns in results file
%

% Iteration and model
rnd_parameters(1) = n_iter_mod;
rnd_parameters(2) = do_Covid;

% Randomize draws (if applicable)
if Loop.do_random == 1
    rng('shuffle');  % ensures that random will not be the same at each run of the script
end

% Start year (random)
if Loop.min_startyear ~= Loop.max_startyear
    startyear = randsample(Loop.min_startyear:Loop.max_startyear,1); % random start year
else
    startyear = Loop.min_startyear; % in this case, min and max are equal
end
rnd_parameters(3) = startyear;

% Variables (random)
if Loop.min_var ~= Loop.max_var
    rnd_parameters(4) = randsample(Loop.min_var:Loop.max_var,1); % random number of variables
else
    rnd_parameters(4) = Loop.min_var; % in this case, min and max are equal
end
sel_var = sort(randsample(1:Loop.n_var_xest,rnd_parameters(4)),2); % selection of Loop.n_var_iter random variables

switch country_model

    case 'DFM'

        % Number of factors (random)
        if Loop.min_r ~= Loop.max_r
            rnd_parameters(5) = randsample(Loop.min_r:Loop.max_r,1);
        else
            rnd_parameters(5) = Loop.min_r; % in this case, min and max are equal
        end
        
        % Number of lags (random)
        if Loop.min_p ~= Loop.max_p
            rnd_parameters(6) = randsample(Loop.min_p:Loop.max_p,1);
        else
            rnd_parameters(6) = Loop.min_p;
        end

        % Block factors (not random)
        rnd_parameters(7) = Par.block_factors; 

    case 'BEQ'

        % Number of lags for monthly regressors (random)
        if Loop.min_lagM ~= Loop.max_lagM
            rnd_parameters(5) = randsample(Loop.min_lagM:Loop.max_lagM,1);
        else
            rnd_parameters(5) = Loop.min_lagM;
        end
        
        % Number of lags for quarterly regressors (random)
        if Loop.min_lagQ ~= Loop.max_lagQ
            rnd_parameters(6) = randsample(Loop.min_lagQ:Loop.max_lagQ,1);
        else
            rnd_parameters(6) = Loop.min_lagQ;
        end

        % Number of lags for the endogenous variable (random)
        if Loop.min_lagY ~= Loop.max_lagY
            rnd_parameters(7) = randsample(Loop.min_lagY:Loop.max_lagY,1);
        else
            rnd_parameters(7) = Loop.min_lagY;
        end

    case 'BVAR'

        % Number of lags (random)
        if Loop.min_bvar_lags ~= Loop.max_bvar_lags
            rnd_parameters(5) = randsample(Loop.min_bvar_lags:Loop.max_bvar_lags,1);
        else
            rnd_parameters(5) = Loop.min_bvar_lags;
        end
        
end

end