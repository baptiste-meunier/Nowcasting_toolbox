% This scripts saves the entire workspace and copies the results in the
% excel file specified in "excel_outputfile". 
% The subfolder output stores the results of each run by
% date in matlab (.mat) format.
%
% Adapted from an initial version by S. Delle Chiaie and F. Kurcz
%

%% Extract nowcast and percentage of available data

% Nowcast_position finds the third month of a quarter
indxGDP = find(~isnan(xest_out(:,end)),1,'last');
lastnow = GDP_track(indxGDP + 3,3); % quarter after last available
forecast = GDP_track(indxGDP + 6,3); % two quarters after last available

% Add percentage of variables available
input_now = xest_out((indxGDP+1):(indxGDP+3),1:(end-1)); % inputs for quarter after last available
share_now = 100 - 100*sum(isnan(input_now),'all')/numel(input_now);
input_fore = xest_out((indxGDP+4):(indxGDP+6),1:(end-1)); % inputs for two quarters after last available
share_fore = 100 - 100*sum(isnan(input_fore),'all')/numel(input_fore);


%% Prepare the MAE 

% If real-time nowcast quarter exceeds model nowcast quarter, then we are backcasting
yr_q_mdl = [datet(indxGDP + 3,1), quarter(datetime(datet(indxGDP + 3,1),datet(indxGDP + 3,2),1))]; % model year and quarter of nowcast
yr_q_real = [year(datetime("today")),quarter(datetime("today"))]; % real-time year and quarter of nowcast
date_comp = sum(yr_q_real > yr_q_mdl);

% Get the MAE
% NB: these are MAE adjusted for outliers (as in ECB, 2009)
if date_comp >= 1 % first prediction is a backcast

    mae_now = MAE.Bac.mae;
    mae_for = MAE.Now.mae;
    fda_now = MAE.Bac.fda;
    fda_for = MAE.Now.fda;

else % first prediction is a nowcast

    mae_now = MAE.Now.mae;
    mae_for = MAE.For.mae;
    fda_now = MAE.Now.fda;
    fda_for = MAE.For.fda;

end 


%% Prepare table for output

% use a cell to store both dates and numeric values in one array
if prev_news == false

    chartdata = {strrep(namesave,'sav_', ''), lastnow, strcat(num2str(round(share_now)),'%'), strcat(num2str(round(100*fda_now)),'%'), ... % current nowcast
                 mae_now,lastnow - mae_now,lastnow + mae_now, ... % MAE
                 [], ...
                 [], [], ... % old nowcast ([] since no prev_news)
                 [], ...
                 nan(1,length(groups_name)),[], ... % last [] is for revisions to data and model
                 nan(1,1), ... % blank space between group and individual news
                 nan(1,length(nameseries))};

    quarterahead_chartdata = {strrep(namesave,'sav_', ''), forecast, strcat(num2str(round(share_fore)),'%'), strcat(num2str(round(100*fda_for)),'%'), ... % current forecast
                              mae_for,forecast - mae_for,forecast + mae_for, ... % MAE
                              [], ...
                              [], [], ... % old forecast ([] since no prev_news)
                              [], ...
                              nan(1,length(groups_name)),[], ... % last [] is for revisions to data and model
                              nan(1,1), ... % blank space between group and individual news
                              nan(1,length(nameseries))};

else

    chartdata = {datestr(news_results.dates(2)), lastnow, strcat(num2str(round(share_now)),'%'), strcat(num2str(round(100*fda_now)),'%'), ... % current nowcast
                 mae_now,lastnow - mae_now,lastnow + mae_now, ... % MAE
                 [], ...
                 datestr(news_results.dates(1)), news_results.now_old, ... % old nowcast
                 [], ...
                 news_results.group_sums, news_results.impact_revision + news_results.impact_reestimation, ... % news decomposition
                 nan(1,1), ... % blank space between group and individual news
                 news_results.indiv_news}; % news decomposition for individual variables

    quarterahead_chartdata = {datestr(news_results_fcst.dates(2)), forecast, strcat(num2str(round(share_fore)),'%'), strcat(num2str(round(100*fda_for)),'%'), ... % current forecast
                              mae_for,forecast - mae_for,forecast + mae_for, ... % MAE
                              [], ...
                              datestr(news_results_fcst.dates(1)), news_results_fcst.now_old, ... % old forecast
                              [], ...
                              news_results_fcst.group_sums, news_results_fcst.impact_revision + news_results_fcst.impact_reestimation, ...
                              nan(1,1), ... % blank space between group and individual news
                              news_results_fcst.indiv_news}; % news decomposition

end


%% Prepare table for (approximate) contributions

% Sort the individual contributions
% NB: we always sort only for the nowcast (first missing quarter)
[ctb_indiv_sorted,idx_sort] = sort(news_results.ctb_indiv(1:end-1)); % not sorting the last because this is the mean
ctb_indiv_fcst_sorted = news_results_fcst.ctb_indiv(idx_sort);

% Create the matrix for contributions
contrib = [sum(news_results.ctb_groups),sum(news_results_fcst.ctb_groups); ...
           NaN,NaN; ...
           news_results.ctb_groups',news_results_fcst.ctb_groups'; ...
           NaN,NaN; ...
           NaN,NaN; ...
           ctb_indiv_sorted,ctb_indiv_fcst_sorted; ... % sort without the mean
           news_results.ctb_indiv(end),news_results_fcst.ctb_indiv(end)]; % adding the mean on the bottom

% Create the names of the contributions
nameseries_sorted = nameseries(idx_sort);
names_contrib = ['Nowcast'; ...
                 NaN; ...
                 string(groups_name)';'Mean'; 
                 NaN; ...
                 NaN; ...
                 string(nameseries_sorted)'; ... % no mean in the names of series
                 'Mean'];

% Sort the full names and groups of the individual series
fullnames_sorted = fullnames(idx_sort);
groups_sorted = groups(idx_sort);
gnames_sorted = cell(size(groups_sorted));
for i = 1:length(groups_sorted)
    gnames_sorted(i) = groups_name(groups_sorted(i));
end

% Create information on individual contributions
info_indiv_contrib = cell(size(contrib,1)-1,1); % -1 because we have no information for the mean (last row of contrib / names_contrib)
num_cell = length(contrib) - length(fullnames); % +1 because of NaN in contrib / names_contrib
info_indiv_contrib(num_cell:end,1) = fullnames_sorted';
info_indiv_contrib(num_cell:end,2) = num2cell(groups_sorted');
info_indiv_contrib(num_cell:end,3) = gnames_sorted';
info_indiv_contrib(num_cell-1,:) = {'Full names','Group ID','Group name'};

% Create a matrix with individual contributions by groups
% NB: this is for use in the Excel charts
matrix_indiv_contrib = cell(size(contrib,1)-1,max(groups_sorted));
for vv = 1:length(ctb_indiv_sorted)
    gg = groups_sorted(vv);
    matrix_indiv_contrib(num_cell+vv-1,gg) = num2cell(ctb_indiv_sorted(vv));
end

% Adding the headers
% NB: also include a warning that individual contributions are for the first prediction (backcast or nowcast)
matrix_indiv_contrib(num_cell-1,:) = groups_name;
matrix_indiv_contrib(num_cell-2,1) = {'Contributions below relate to first prediction (back-cast or nowcast)'};


%% Save the entire workspace in excel and store in the output folder
common_save_excel(chartdata,quarterahead_chartdata,Par,excel_outputfile,newsfile,prev_news,groups_name,contrib,names_contrib,info_indiv_contrib,matrix_indiv_contrib,nameseries,fullnames,heatmap,datet,country,Res,groups,do_range,range,indxGDP);


%% Save the entire workspace in matlab
% ...results in the summary file 'country.name'_GDP_tracking.xlsx

% Rename variables to save them
Qend = Par.qend;
%namesave = namesave_orig;

% Clean workspace
clear Qend_* namesave_* flag* store* Nowcast_position pos Yearq check_spelling

% Autorename the workspace if the name already exists
[namesave] = common_autorename(namesave,outputfolder,rootfolder);
save(strcat(outputfolder,namesave)); 

if strcmp('cur_nowcast.mat',newsfile)
    % Save old nowcast as old_nowcast
    if exist(strcat(outputfolder,'cur_nowcast.mat'),'file') == 2
        movefile(strcat(outputfolder,'cur_nowcast.mat'),strcat(outputfolder,'old_nowcast.mat'));
    end
    
    % Save current output as cur_nowcast
    save(strcat(outputfolder,'cur_nowcast')); 
end

% Create a vintage the output Excel for the nowcast date
vintage_excel = [strrep(namesave,'.mat',''),'_',country.name,'_tracking.xlsx'];
filename_vintage = fullfile(outputfolder, vintage_excel);
copyfile(excel_outputfile,filename_vintage);
