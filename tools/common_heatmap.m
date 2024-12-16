function [heatmap] = common_heatmap(xest,Par,groups,groups_name,fullnames)
% This scripts computes the heatmap = z-scores for input data entering in the model 
% Variables are re-organized by groups
%
% INPUT 
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
%
% OUTPUT
% - heatmap [stucture]
%       o zscores [matrix] = z-scores for the actual data (re-ordered by group)
%       o names [cell array] = tables with name of groups and variables (re-ordered by group)
%       o zscores_agg [matrix] = z-scores for the actual data averaged (unweighted) by group
%       o names_agg [cell array] = tables with name of groups
%

%% STEP 1. CHECK CONSISTENCY

% Check 1 - dimensions of groups and fullnames match
if length(groups) ~= length(fullnames)
    error('Dimensions of groups and fullnames do not match')
end

% Check 2 - dimensions of groups and xest match
if length(groups) ~= size(xest,2)
    error('Dimensions of groups and xest (number of columns) do not match')
end

% Check 3 - number of groups and groups_name match
if max(groups) > length(groups_name)
    error('Some values in groups (ID of groups) are higher than the number of entries in groups_name (names of groups)')
end


%% STEP 2. COMPUTE Z-SCORES

% Mean and standard_deviations
z_m = mean(xest,'omitnan');
z_s = std(xest,'omitnan');

% Fill the data for quarterly values
heatmap_init = fillmissing(xest,'next','MaxGap',3);

% Smooth the m-o-m variations as in the DFM
% This is based on the formula used in Mariano, R., and Murasawa, Y. (2003), "A new coincident index of business cycles based on monthly and quarterly series," Journal of Applied Econometrics, 18, 427–443
smoothed_heatmap = 1/9*(heatmap_init(1:end-4,1:Par.nM) + 2*heatmap_init(2:end-3,1:Par.nM) + ...
                   3*heatmap_init(3:end-2,1:Par.nM) + 2*heatmap_init(4:end-1,1:Par.nM) + ...
                   heatmap_init(5:end,1:Par.nM));
smoothed_heatmap = [nan(4,Par.nM);smoothed_heatmap];
heatmap_init(:,1:Par.nM) = smoothed_heatmap;

% Standardised values - for the actual data
zscores = (heatmap_init - repmat(z_m,size(heatmap_init,1),1))./repmat(z_s,size(heatmap_init,1),1);


%% STEP 3. RE-ORGANIZE THE HEATMAP BY GROUPS

% Prepare containers
heatmap.names = cell(size(xest,2),2);
heatmap.zscores = nan(size(zscores));
counter = 1;

% Add the target variable on the bottom
heatmap.names(end,1) = groups_name(groups(end));
heatmap.names(end,2) = fullnames(end);
heatmap.zscores(:,end) = zscores(:,end);

% Re-organize values and variable names
groups_ID = unique(groups);
for vv = 1:length(groups_ID)

    % Get indices of variables belonging to the group
    ID_vv = groups_ID(vv);
    idx_vv = find(groups(1:end-1)==ID_vv);

    % Add the variables to the matrix
    for ww = 1:length(idx_vv)

        num_ww = idx_vv(ww); % column number of the variable
        heatmap.zscores(:,counter) = zscores(:,num_ww);
        heatmap.names(counter,1) = groups_name(ID_vv);
        heatmap.names(counter,2) = fullnames(num_ww);
        counter = counter + 1;

    end

end


%% STEP 4. COMPUTE THE AVERAGE HEATMAPS

% Prepare containers
% +1 because target variable is isolated in the bottom line
heatmap.names_agg = cell(length(groups_name)+1,1);
heatmap.zscores_agg = nan(size(zscores,1),length(groups_name)+1);
counter = 1;

% Add the target variable on the bottom
heatmap.names_agg(end,1) = fullnames(end);
heatmap.zscores_agg(:,end) = zscores(:,end);

% Re-organize values and variable names
groups_ID = unique(groups);
for vv = 1:length(groups_ID)

    % Get indices of variables belonging to the group
    ID_vv = groups_ID(vv);
    idx_vv = groups(1:end-1)==ID_vv; % -1 to avoid the target variable (in last position)

    % Compute the average 
    temp_zscores = zscores(:,idx_vv);
    temp_mean_zscores = mean(temp_zscores,2,"omitnan");

    % Add to the matrix
    heatmap.zscores_agg(:,counter) = temp_mean_zscores;
    heatmap.names_agg(counter) = groups_name(ID_vv);
    counter = counter + 1;

end

% Check 1 - Dimension of heatmap.zscores
if size(heatmap.zscores,1) ~= size(xest,1)
    error('Number of observation do not match between data and zscores')
end

% Check 2 - Dimension of heatmap.zscores_agg
if size(heatmap.zscores_agg,1) ~= size(xest,1)
    error('Number of observation do not match between data and zscores_agg')
end



end
