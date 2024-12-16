function[news_results_now,news_results_fcst,prev_news] = BEQ_News_Mainfile(outputfolder,newsfile,groups_name,X_new,Res_new,Par_new,namesave,datet,nameseries,model_new)
% This script computes the news decomposition for bridge equation 
% 
% INPUTS 
% - outputfolder [string] =  name of the output folder (usually the country name)
% - newsfile [string] = name of mat file to compare with for news decomposition
% - groups_name [cell array] = name of the groups (for summing contributions / decomposition)
% - X_new [matrix] = dataset used for current BEQ (with NaN)
% - Res_new [structure] = results of (current) BEQ
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
% - Par_new [structure] = parameters of (current) DFM
%       o startyear [scalar] = starting year of the estimation sample
%       o startmonth [scalar] = starting month of the estimation sample
%       o blocks [vector/matrix] = block identification
%       o r [scalar] = number of estimated factors
%       o block_factors [scalar] = switch on whether to do block factors (=1) or not (=0)
%       o p [scalar] = number of lags
%       o nM [scalar] = number of monthly variables
%       o nQ [scalar] = number of quarterly variables
%       o idio [scalar] = idiosyncratic component specification for DFM: 1 = AR(1), 0 = iid
%       o yearq [scalar] = year for which latest GDP is available
%       o qend [scalar] = quarter for which latest GDP is available
%       o thresh [scalar] = threshold for convergence of EM algorithm in DFM
%       o max_iter [scalar] = number of iterations for the EM algorithm in DFM
%       o lagM [scalar] = number of lags for monthly variables (at quarterly frequency)
%       o lagQ [scalar] = number of lags for quarterly variables
%       o lagY [scalar] = number of lags for endogenous variable
%       o type [scalar] = type of interpolation performed (see BEQ_estimate for details)
% - datet [matrix] = dates (year / month)
% - nameseries [cell vector] = ID of the series in Haver
% - model_new [string] = type of model used
%
% OUTPUTS
% - news_results_now [structure] = news decomposition for the nowcast
%       o dates [cell vector] = dates of previous and current nowcasts
%       o impact_revision [scalar] = effect of data revision
%       o impact_news [vector] = effect of news
%       o impact_reestimation [scalar] = effect of parameter re-estimation
%       o group_sums [vector] = effect of news by groups
%       o now_newdata [scalar] = current nowcast with data changed
%       o now_newpar [scalar] = current nowcast with parameters changed
%       o now_newdata [scalar] = old nowcast
%       o indx_news [vector] = identification of series with news
%       o ctb_indiv [vector] = contribution of individual series to nowcast (last value is mean + carry-over)
%       o ctb_groups [vector] = contribution of series group to nowcast (last value is mean + carry-over)
% - news_results_fcst [structure] = news decomposition for forecast
% - prev_news [logical] = whether the previous comparison is available
%


%% Check if previous newsfile exist
prev_news = exist(strcat(outputfolder,newsfile),'file'); % prev_news is a logical, false = 0; true = 2 for mat. files

%% if previous newsfile exist, store variables needed for the news decomposition

if  prev_news == false

    news_results_now = []; news_results_fcst = [];
    disp('Section 5: News decomposition cannot be computed as a previous nowcast does not exist');

    % Check if the file does not exist because of spelling error
    check_spelling = (exist(strcat(outputfolder,'cur_nowcast'),'file') == 2) && (exist(strcat(outputfolder,newsfile),'file') ~= 2);
    if check_spelling
        warning('The specified newsfile does not exist - but perhaps the .mat file is misspelled as a cur_nowcast.mat file exists.')
    end

    % rename current variables (for contributions)
    Qend_new = Par_new.qend;
    Yearq = Par_new.yearq;

else

    disp('Section 5: Previous nowcast found')
    
    % Rename current variables
    Qend_new = Par_new.qend;
    Yearq = Par_new.yearq;
    newdate = strrep(namesave,'sav_', '');
    nameseries_new = nameseries;
    type_new = Par_new.type;
    
    % Load old vintage of data
    load(strcat(outputfolder,newsfile),'xest_out','namesave','Res','Par','nameseries','country')    
    X_old = xest_out;
    Res_old = Res;
    Par_old = Par;                              
    olddate = strrep(namesave,'sav_', '');      % get the date at which the model was last run
    olddate = strrep(olddate,'.mat', '');
    if ~isfield(Par_old,'type')
        type_old = 0; % if type is not known, then it is likely not a BEQ
    else
        type_old = Par_old.type;
    end
    if ~isfield(country,'model')
        model_old = 'DFM'; % if the information on model is not available, then it is an older version of the code and therefore a DFM
    else
        model_old = country.model;
    end

    % In case on the last date the platform was run several times, namesave
    % is e.g. date_2.mat; get rid of "_2":
    if length(olddate)>11; olddate = olddate(1:11); end  
    clear xest namesave Res Qend Par
    
    % Check 1 - blocks_name input compatibility
    if size(unique(Res_new.groups(:)),1) ~= size(groups_name(:),1)
        error('The number of groups is not identical to the number of group names. Possibly some groups are empty (no variables)')
    end

    % Check 2 - whether input is chronologically correct
    if datenum(olddate) > datenum(newdate)
        warning('The old dataset is associated with a more recent date than the new dataset. Check the inputs to this function.')
    end

    % Check 3 - model types are the same
    if ~strcmp(model_old,model_new)
        disp('Model types do not match.')
        disp('It is not possible to compute the news relative to a different model.')
        news_results_now = []; 
        news_results_fcst = []; 
        prev_news = false;
    end
    
    % Check 4 - specifications of models (lagM, lagQ, lagY, nM) are identical
    if ~(Par_old.lagM == Par_new.lagM && Par_old.lagQ == Par_new.lagQ && Par_old.lagY == Par_new.lagY && Par_old.nM == Par_new.nM && isequal(Par_new.Dum,Par_old.Dum) && type_old == type_new)
        disp('Model specifications (nb. of lags, nb. of variables, type) do not match.')
        disp('It is not possible to compute the news relative to a different specification.')
        news_results_now = []; 
        news_results_fcst = []; 
        prev_news = false;
    end

    % Same number of names for the variables and they are identical
    if length(nameseries_new) == length(nameseries)
        if ~(sum(string(nameseries_new)~=string(nameseries))==0)
            disp('Names across input series do not match.')
            disp('It is not possible to compute the news relative to a different model.')
            news_results_now = []; 
            news_results_fcst = []; 
            prev_news = false;
        end
    else
        disp('Number of names for input series do not match.')
        disp('It is not possible to compute the news relative to a different model.')
        news_results_now = []; 
        news_results_fcst = []; 
        prev_news = false;
    end

    % Check 6 - starts of estimation (startyear/startmonth) are identical
    if ~isfield(Par_old,'startyear') || ~isfield(Par_old,'startmonth') % in case there is none, we suppose they are identical
        Par_old.startyear = Par_new.startyear;
        Par_old.startmonth = Par_new.startmonth;
    end
    if ~(Par_old.startyear == Par_new.startyear && Par_old.startmonth == Par_new.startmonth)
        disp('Starting date of the estimation sample (startyear, startmonth) do not match.')
        disp('It is not possible to compute the news relative to a different specification.')
        news_results_now = []; 
        news_results_fcst = []; 
        prev_news = false;
    end

end

if prev_news ~= false

    %% Prepare inputs for news_function
    
    % Check row size of data matrices
    T_old = size(X_old,1); 
    T_new = size(X_new,1);
    N = size(X_new,2);
    if T_new > T_old
        X_old = [X_old; NaN(T_new-T_old,N)];
    elseif T_new < T_old
        error('The time dimension of the old dataset is larger than in the new dataset.')
    end
    
    % GDP index:
    gdp_index = size(X_new,2);    
    
    % For impact due to revision
    X_rev = X_new;  
    X_rev(isnan(X_old)) = NaN;  
    
    % Find the index of the (month of the) quarter we're interested in nowcasting:
    if Qend_new == 4
        Yearq = Yearq + 1;
    end     
    targetmonth = find(datet(:,1) == Yearq & datet(:,2) == (mod(Qend_new,4)+1)*3 );
    
    %% Start news computation for nowcast and for forecast
    
    for forecast_news = 0:1

        targetmonth = targetmonth + forecast_news*3;    % adjust targetmonth for forecast       
    
        %% Compute the three components of the change in the nowcast: due to news, data revision, or parameter re-estimation
    
        % store old nowcast to compare to nowcast that takes data revisions
        % into account:
        if size(Res_old.X_sm,1)>=targetmonth
            now_old = Res_old.X_sm(targetmonth,gdp_index);
        else
            now_old = nan;
        end
    
        % Calculate impact of news using old parameters
        datet_in = datet;
        re_estimate = 0;
        Res_revdata_oldparam = BEQ_estimate(X_rev,Par_new,datet_in,nameseries_new,re_estimate,Res_old.coeffs);
        now_rev = Res_revdata_oldparam.X_sm(targetmonth,gdp_index);
        impact_revision = now_rev - now_old;
        ctb_revdata = Res_revdata_oldparam.contrib_X_sm;

        % Calculate the impact of news using the old parameters
        Res_newdata_oldparam = BEQ_estimate(X_new,Par_new,datet_in,nameseries_new,re_estimate,Res_old.coeffs);
        now_newdata = Res_newdata_oldparam.X_sm(targetmonth,gdp_index);
        ctb_newdata = Res_newdata_oldparam.contrib_X_sm;

        % Compute the respective impact
        impact_news = ctb_newdata(targetmonth,:) - ctb_revdata(targetmonth,:);
        impact_news = impact_news(1:length(nameseries_new)); % deleting the "news" from constant and dummies (if any) as these are just 0
    
        % Compute impact of re-estimation
        % This is the difference between the nowcast new data / old parameters 
        % and the nowcast new data / new parameters
        now_new = Res_new.X_sm(targetmonth,gdp_index);
        impact_reestimation = now_new - now_newdata;

        % Aggregate contributions by series
        news_results.ctb_indiv = Res_new.contrib_X_sm(targetmonth,:)'; % adding the nowcast with NaN (close to mean) as an additional contribution
        if size(Par_new.Dum,1) > 0
            news_results.ctb_indiv(end-1) = news_results.ctb_indiv(end-1) + news_results.ctb_indiv(end); % adding the contribution of the dummies
            news_results.ctb_indiv(end) = []; % deleting the contribution of the dummies
        end

        % Aggregate contributions by groups
        group_id = unique(Res_new.groups);
        news_results.ctb_groups = zeros(size(group_id));
        for jj = 1:length(group_id)
            news_results.ctb_groups(jj) = sum(news_results.ctb_indiv(Res_new.groups==group_id(jj)),'omitnan');
        end
        if size(Par_new.Dum,1) > 0
            news_results.ctb_groups(jj+1) = news_results.ctb_indiv(end-1) + news_results.ctb_indiv(end); % adding the contribution of the constant and dummies
        else
            news_results.ctb_groups(jj+1) = news_results.ctb_indiv(end); % adding only the contribution of the constant   
        end

        % If there are no news, impact_news is an empty matrix
        if sum(impact_news,'omitnan') ~= 0            
    
            % For legibility, create groups that are plotted. Details to the specific
            % series can be obtained from the table, see below
            group_id = unique(Res_new.groups);
            group_sums = zeros(size(group_id));
            for jj = 1:length(group_id)
                group_sums(jj) = sum(impact_news(Res_new.groups==group_id(jj)),'omitnan');
            end
    
            % Create also the impact for individual variables
            indiv_news = impact_news;

            % Store output
            news_results.dates = {olddate, newdate};
            news_results.impact_revision = impact_revision;
            news_results.impact_news = impact_news;
            news_results.impact_reestimation = impact_reestimation;
            news_results.group_sums = group_sums;
            news_results.indiv_news = indiv_news;
            news_results.now_newdata = now_newdata;
            news_results.now_newpar = now_new; 
            news_results.now_old = now_old;
        
        else
            if ~forecast_news
                fprintf('There are no news (but possibly data revisions).\n')     % Explain a potential change in the Nowcast (only once).
            end
            
            news_results.dates = {olddate, newdate};
            news_results.impact_revision = impact_revision;
            news_results.impact_news = [];
            news_results.impact_reestimation = impact_reestimation;
            news_results.group_sums = zeros(1,size(unique(Res_new.groups(:)),1));
            news_results.indiv_news = zeros(1,length(Res_new.series));
            news_results.now_newdata = now_newdata;
            news_results.now_newpar = now_new; 
            news_results.now_old = now_old;
        end         
        
        if ~forecast_news
            news_results_now = news_results;
        else
            news_results_fcst = news_results;
        end
    
    end

    disp ('Section 5: News decomposition and contributions completed');

else % we still need to compute the contributions 

    % Find the index of the (month of the) quarter we're interested in nowcasting:
    if Qend_new == 4; Yearq = Yearq + 1; end     
    targetmonth = find(datet(:,1) == Yearq & datet(:,2) == (mod(Qend_new,4)+1)*3 );

    for forecast_news = 0:1

        targetmonth = targetmonth + forecast_news*3;    % adjust targetmonth for forecast
    
        % Aggregate contributions by series
        news_results.ctb_indiv = Res_new.contrib_X_sm(targetmonth,:)'; % adding the nowcast with NaN (close to mean) as an additional contribution
        if size(Par_new.Dum,1) > 0
            news_results.ctb_indiv(end-1) = news_results.ctb_indiv(end-1) + news_results.ctb_indiv(end); % adding the contribution of the dummies
            news_results.ctb_indiv(end) = []; % deleting the contribution of the dummies
        end 

        % Aggregate contributions by groups
        group_id = unique(Res_new.groups);
        news_results.ctb_groups = zeros(size(group_id));
        for jj = 1:length(group_id)
            news_results.ctb_groups(jj) = sum(news_results.ctb_indiv(Res_new.groups==group_id(jj)),'omitnan');
        end
        if size(Par_new.Dum,1) > 0
            news_results.ctb_groups(jj+1) = news_results.ctb_indiv(end-1) + news_results.ctb_indiv(end); % adding the contribution of the constant and dummies
        else
            news_results.ctb_groups(jj+1) = news_results.ctb_indiv(end); % adding only the contribution of the constant   
        end     

        if ~forecast_news
            news_results_now = news_results;
        else
            news_results_fcst = news_results;
        end
    
    end

    disp ('Section 5: No news decomposition but contributions completed');

end % end of condition on prev_news

end