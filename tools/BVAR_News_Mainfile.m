function[news_results_now,news_results_fcst,table_now,table_fore,prev_news] = BVAR_News_Mainfile(outputfolder,newsfile,groups_name,X_new,Res_new,Par_new,namesave,datet,nameseries,model_new)
% This script computes the news decomposition. 
% 
% Code adapted from Bańbura, M., and Modugno, M. (2014). "Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%


%% Check if previous newsfile exist
prev_news = exist(strcat(outputfolder,newsfile),'file'); % prev_news is a logical, false = 0; true = 2 for mat. files

%% if previous newsfile exist, store variables needed for the news decomposition

if  prev_news == false

    news_results_now = []; news_results_fcst = []; table_now = []; table_fore = [];
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
    
    % Load old vintage of data
    load(strcat(outputfolder,newsfile),'xest_out','namesave','Res','Par','nameseries','country')
    X_old = xest_out;
    Res_old = Res;
    Par_old = Par;                              
    olddate = strrep(namesave,'sav_', '');      % get the date at which the model was last run
    olddate = strrep(olddate,'.mat', '');
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
        table_now = []; 
        table_fore = [];
        prev_news = false;
    end
    
    % Check 4 - specifications of models (lags) are identical
    if ~(Par_old.bvar_lags == Par_new.bvar_lags)
        disp('Model specifications (lags) do not match.')
        disp('It is not possible to compute the news relative to a different specification.')
        news_results_now = []; 
        news_results_fcst = []; 
        table_now = []; 
        table_fore = [];
        prev_news = false;
    end

    % Same number of names for the variables and they are identical
    if length(nameseries_new) == length(nameseries)
        if ~(sum(string(nameseries_new)~=string(nameseries))==0)
            disp('Number of names for input series do not match.')
            disp('It is not possible to compute the news relative to a different model.')
            news_results_now = []; 
            news_results_fcst = []; 
            table_now = []; 
            table_fore = [];
            prev_news = false;
        end
    else
        disp('Number of names for input series do not match.')
        disp('It is not possible to compute the news relative to a different model.')
        news_results_now = []; 
        news_results_fcst = []; 
        table_now = []; 
        table_fore = [];
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
        table_now = []; 
        table_fore = [];
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
    basemonth = targetmonth;
    
    %% Start news computation for nowcast and for forecast
    
    for forecast_news = 0:1

        targetmonth = targetmonth + forecast_news*3;    % adjust targetmonth for forecast       
    
        %% Compute the three components of the change in the nowcast: Due to News, data revision, or parameter reestimation
    
        % store old nowcast to compare to nowcast that takes data revisions
        % into account:
        if size(Res_old.X_sm,1)>=targetmonth
            now_old = Res_old.X_sm(targetmonth,gdp_index);
        else
            now_old = nan;
        end

        % Calculate impact of news using old parameters
        try
            [now_rev,now_newdata,~,actual,forecast,weight,~,v_miss,~] = ...
                BVAR_compute_news(X_rev,X_new,Res_old,Par_new,targetmonth/3,gdp_index);
        catch
            disp('Not running BVAR news - due to the BVAR not converging in the first place.')
            disp('Please try another specification.')
            news_results_now = []; 
            news_results_fcst = []; 
            table_now = []; 
            table_fore = [];
            prev_news = false;
            news_results.ctb_indiv = zeros(length(Res_new.series),1);
            news_results.ctb_groups = zeros(size(group_id));
            return
        end
    
        % Calculate the impact of news again using the new parameters;
        % the difference to the nowcast 'now_newdata' is then due to reestimation of parameters
        [~,now_newpar,~,~,~,~] = BVAR_compute_news(X_rev,X_new,Res_new,Par_new,targetmonth/3,gdp_index);

    
        % Compute the respective impact
        impact_revision = now_rev - now_old;
        impact_reestimation = now_newpar - now_newdata;
        news = actual - forecast;
        impact_news = weight .* news;         % length equal to number of new datapoints


        % -> now_newpar - now_old = sum(impact_news) + impact_revision + impact_reestimation                   
    
        % Compute contributions
        % NB: this is an approximate solution where contributions are the
        % news between a nowcast based on actual dataset and a nowcast based
        % on hypothetical dataset with NaN
        X_ctb = X_new;
        X_ctb(basemonth-24:end,:) = NaN; % putting 24 months before target to NaN so that forecast (ctb_now_rev) is close to mean

        [ctb_now_rev,~,~,ctb_actual,ctb_forecast,ctb_weight,~,ctb_v_miss,~] = ...
            BVAR_compute_news(X_ctb,X_new,Res_new,Par_new,targetmonth/3,gdp_index);

        impact_ctb = ctb_weight .* (ctb_actual - ctb_forecast);
         
     
        % Aggregate contributions by series
        series_id = 1:length(Res_new.series);
        news_results.ctb_indiv = zeros(length(Res_new.series),1);
        for jj = 1:length(series_id)
            news_results.ctb_indiv(jj) = sum(impact_ctb(logical(ctb_v_miss==series_id(jj))),'omitnan');
        end
        news_results.ctb_indiv(jj+1) = ctb_now_rev; % adding the nowcast with NaN (close to mean) as an additional contribution

        % Aggregate contributions by groups
        group_id = unique(Res_new.groups);
        news_results.ctb_groups = zeros(size(group_id));
        for jj = 1:length(group_id)
            news_results.ctb_groups(jj) = sum(impact_ctb(logical(Res_new.groups(ctb_v_miss)==group_id(jj))),'omitnan');
        end
        news_results.ctb_groups(jj+1) = ctb_now_rev; % adding the nowcast with NaN (close to mean) as an additional contribution    

        

        % If there are no news, impact_news is an empty matrix
        if ~isempty(impact_news)            
    
            % For legibility, create groups that are plotted. Details to the specific
            % series can be obtained from the table, see below
            group_id = unique(Res_new.groups);
            group_sums = zeros(size(group_id));
            for jj = 1:length(group_id)
                group_sums(jj) = sum(impact_news(logical(Res_new.groups(v_miss)==group_id(jj))),'omitnan');
            end
    
            % Create also the impact for individual variables
            indiv_news = zeros(1,length(Res_new.series));
            for jj = 1:length(indiv_news)
                indiv_news(jj) = sum(impact_news(logical(v_miss==jj)),'omitnan');
            end

            % Store output
            news_results.dates = {olddate, newdate};
            news_results.impact_revision = impact_revision;
            news_results.impact_news = impact_news;
            news_results.impact_reestimation = impact_reestimation;
            news_results.group_sums = group_sums;
            news_results.indiv_news = indiv_news;
            news_results.now_newdata = now_newdata;
            news_results.now_newpar = now_newpar; 
            news_results.now_old = now_old;
            news_results.actual = actual;
            news_results.forecasts = forecast;
            news_results.weight = weight;
            news_results.indx_news = v_miss;
    
            % Print output to the Command Window            
            % names of the series as a variable in the table is necessary if a
            % series has two or more new data points, because RowNames need to be
            % unique.
            tablevariables = Res_new.name_descriptor(v_miss)';
            tablevariables{end+1} = 'Revision';
            tablevariables{end+1} = 'Reestimation';
            T1 = table(tablevariables(:), 'VariableNames', {'Series'});
    
            % Note that the forecast values displayed here were computed using the
            % revised data, therefore they may not be identical to the ones in
            % Res_old
            T2 = table([actual; zeros(2,1)], [forecast;zeros(2,1)], [weight;zeros(2,1)] ,[impact_news; impact_revision; impact_reestimation],...
                'VariableNames', {'Actual', 'Predicted_value', 'Weight', 'Impact'});
            
            if forecast_news
                table_fore = [T1,T2];
            else
                table_now = [T1,T2];
            end
            
            
            if ~forecast_news             % only print for nowcast
                fprintf('\n\n')
                disp(repmat('=',1,80))
                fprintf('News impact calculation: Data Vintage %s compared to vintage %s \n',datestr(newdate, 'mmmm dd, yyyy'),datestr(olddate, 'mmmm dd, yyyy'))
    
                % Check difference:
                discpr = abs(now_newpar - now_old) - abs(sum(impact_revision + sum(group_sums) + impact_reestimation));
                if (discpr >1e-10)
                    warning('Difference to last Nowcast does not sum to the decomposition. Discrepancy: %5.10f',discpr)
                end
                fprintf('\n')
                fprintf('Updated Nowcast based on the %s vintage: %f \n',datestr(newdate, 'mmmm dd, yyyy'),now_newpar)
                fprintf('Old Nowcast based on the %s vintage: %f \n',datestr(olddate, 'mmmm dd, yyyy'),now_old)
                fprintf('\n')
    
                disp('Details on the series'' with new data points:') 
                disp(table_now)
    
                % Display table with group sums
                T3 = table([group_sums'; impact_revision + impact_reestimation],'VariableNames',...
                    {'Groups_impact'},'Rownames', [groups_name, {'Revision and re-estimation'}]');
                disp(T3)
                disp(repmat('=',1,80))
            end
    
        else
            if ~forecast_news
                fprintf('There are no news (but possibly data revisions).\n')     % Explain a potential change in the Nowcast (only once).
            end
            table_now = []; table_fore = [];
            news_results.dates = {olddate, newdate};
            news_results.impact_revision = impact_revision;
            news_results.impact_news = [];
            news_results.impact_reestimation = impact_reestimation;
            news_results.group_sums = zeros(1,size(unique(Res_new.groups(:)),1));
            news_results.indiv_news = zeros(1,length(Res_new.series));
            news_results.now_newdata = now_newdata;
            news_results.now_newpar = now_newpar; 
            news_results.now_old = now_old;
            news_results.indx_news = v_miss;
        end         
        
        if ~forecast_news
            news_results_now = news_results;
        else
            news_results_fcst = news_results;
        end
    
    end

    disp ('Section 5: News decomposition and contributions completed');

else

    % GDP index:
    gdp_index = size(X_new,2);    


    % Find the index of the (month of the) quarter we're interested in nowcasting:
    if Qend_new == 4; Yearq = Yearq + 1; end     
    targetmonth = find(datet(:,1) == Yearq & datet(:,2) == (mod(Qend_new,4)+1)*3 );
    basemonth = targetmonth;

    for forecast_news = 0:1

        targetmonth = targetmonth + forecast_news*3;    % adjust targetmonth for forecast
    
        % Compute contributions
        % NB: this is an approximate solution where contributions are the
        % news between a nowcast based on actual dataset and a nowcast based
        % on hypothetical dataset with NaN

        X_ctb = X_new;
        X_ctb(basemonth-24:end,:) = NaN; % putting 24 months before target to NaN so that forecast (ctb_now_rev) is close to mean

        try
            [ctb_now_rev,~,~,ctb_actual,ctb_forecast,ctb_weight,~,ctb_v_miss,~] = ...
                BVAR_compute_news(X_ctb,X_new,Res_new,Par_new,targetmonth/3,gdp_index);
        catch
            disp('Not running BVAR contributions - due to the BVAR not converging in the first place.')
            disp('Please try another specification.')
            news_results_now.ctb_indiv = zeros(length(Res_new.series)+1,1); % +1 for the mean
            news_results_now.ctb_groups = zeros(size(unique(Res_new.groups)));
            news_results_fcst.ctb_indiv = zeros(length(Res_new.series)+1,1);
            news_results_fcst.ctb_groups = zeros(size(unique(Res_new.groups))); % +1 for the mean
            return
        end

        impact_ctb = ctb_weight .* (ctb_actual - ctb_forecast);


        % Aggregate contributions by series
        series_id = 1:length(Res_new.series);
        news_results.ctb_indiv = zeros(length(Res_new.series),1);
        for jj = 1:length(series_id)
            news_results.ctb_indiv(jj) = sum(impact_ctb(logical(ctb_v_miss==series_id(jj))),'omitnan');
        end
        news_results.ctb_indiv(jj+1) = ctb_now_rev; % adding the nowcast with NaN (close to mean) as an additional contribution

        % Aggregate contributions by groups
        group_id = unique(Res_new.groups);
        news_results.ctb_groups = zeros(size(group_id));
        for jj = 1:length(group_id)
            news_results.ctb_groups(jj) = sum(impact_ctb(logical(Res_new.groups(ctb_v_miss)==group_id(jj))),'omitnan');
        end
        news_results.ctb_groups(jj+1) = ctb_now_rev; % adding the nowcast with NaN (close to mean) as an additional contribution       

        if ~forecast_news
            news_results_now = news_results;
        else
            news_results_fcst = news_results;
        end
    
    end

    disp ('Section 5: No news decomposition but contributions completed');

end % end of condition on prev_news





end