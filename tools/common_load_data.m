function[Par,xest,t_m,groups,nameseries,blocks,blocks_name,full_names,datet,Loop] = common_load_data(excel_datafile,mon_freq,quar_freq,blocks_sheet,Par_in,m,do_loop,date_today,Loop)
% This script loads data
% fetching from Haver does not work currently
% fetching should be from an excelsheet so exceldata should be 1 
% 
% Adapted from an initial version by S. Delle Chiaie and F. Kurcz
%
% INPUTS
% - excel_datafile [string] = name of the excel with data
% - mon_freq [string] = name of the sheet (in excel_datafile) with monthly indicators
% - quar_freq [string] = name of the sheet (in excel_datafile) with quarterly indicators
% - blocks_sheet [string] = name of the sheet (in excel_datafile) with blocks
% - Par_in [structure] = initial parameters of DFM
%       o startyear [scalar] = starting year of the estimation sample
%       o startmonth [scalar] = starting month of the estimation sample
%       o block_factors [scalar] = switch on whether to do block factors (=1) or not (=0)
%       o r [scalar] = number of estimated factors
%       o p [scalar] = number of lags
%       o idio [scalar] = idiosyncratic component specification: 1 = AR(1), 0 = iid
%       o thresh [scalar] = threshold for convergence of EM algorithm
%       o max_iter [scalar] = number of iterations for the EM algorithm
% - m [scalar] = number of months ahead
% - do_loop [scalar] = switch on whether to do a loop (1/2) or not (0)
% - date_today [vector] = date of the nowcast (year / month)
%
% OUTPUTS
% Data starting after startyear/startmonth is trimmed
% - Par [structure] = parameters of DFM
%       o block_factors [scalar] = switch on whether to do block factors (=1) or not (=0)
%       o r [scalar] = number of estimated factors
%       o p [scalar] = number of lags
%       o idio [scalar] = idiosyncratic component specification: 1 = AR(1), 0 = iid
%       o thresh [scalar] = threshold for convergence of EM algorithm
%       o max_iter [scalar] = number of iterations for the EM algorithm
%       o blocks [vector/matrix] = block identification
%       o nM [scalar] = number of monthly variables
%       o nQ [scalar] = number of quarterly variables
%       o yearq [scalar] = year for which latest GDP is available
%       o qend [scalar] = quarter for which latest GDP is available
%       o trf [structure] = transformation of the variable
%       o blocks_full [vector/matrix] = block identification (for do_loop=2)
%       o size_data_in [scalar] = lengths of raw data set
% - xest [matrix] = values
% - t_m [matrix] = dates (year / month)
% - groups [vector] = identification of groups
% - nameseries [cell vector] = ID of the series in Haver
% - blocks [matrix] = block indentification (raw, can differ from blocks in Par)
% - blocks_name [cell vector] = name of the blocks
% - full_names [cell vector] = complete name of the series
% - datet [matrix] = dates (year / month). Difference with t_m is that datet relates to the data when adding NaN for months ahead
%

% Set parameters
Par = Par_in;


%% Load data from Excel file

[A,B] = xlsread(excel_datafile,mon_freq);
[C,D] = xlsread(excel_datafile,quar_freq);
[blocks,~] = xlsread(excel_datafile,blocks_sheet);
[~,G] = xlsread(excel_datafile,'Groups');
transf_m = A(1,:);                  % load the transformation used for monthly series
transf_q = C(1,:);                  % load the transformation used for quarterly series
groups = [A(2,:) C(2,:)];           % load group definition
nameseries = [B(3,2:end) D(3,2:end)];
full_names = [B(4,2:end) D(4,2:end)];
blocks_name = G;

% Check dimensions
if size(A,1) ~= size(B,1)
    error("Dimensions of dates and series do not match for monthly data. Please check the input file")
end
if size(C,1) ~= size(D,1)
    error("Dimensions of dates and series do not match for quarterly data. Please check the input file")
end

% Data matrices untransformed
seriesm = A(5:end,:);                  % monthly
seriesq = C(5:end,:);                  % quarterly

% Check for discontinued and suspended series
common_CheckDataAvailability(full_names);              


%% Convert dates

[Year_m, Month_m] = datevec(datetime(B(5:end,1), 'InputFormat', 'dd/MM/yyyy'));
t_m = [Year_m, Month_m];
[Year_q, Month_q] = datevec(datetime(D(5:end,1), 'InputFormat', 'dd/MM/yyyy'));
t_q = [Year_q, Month_q];
        

%% Transformations

% Transforming data
data_m = common_transform_data(seriesm,transf_m,1);
data_q = common_transform_data(seriesq,transf_q,1);

% Computing dates
t_m = t_m(2:end,:); % adjust for growth rates
t_q = t_q(2:end,:);

% Ensure the data / dates goes until date_today
% Notably useful for out-of-sample exercises because the publication delays are inferred from last data point
% So last data point needs to match the intended date (otherwise it will wrongly infer publication delays) 
end_t_m = find(t_m(:,1)==date_today(1) & t_m(:,2)==date_today(2), 1);
end_t_q = find(t_q(:,1)==date_today(1) & t_q(:,2)==date_today(2), 1);
if isempty(end_t_m) && isempty(end_t_q)

    % Number of missing months
    nmth_m = 12*t_m(end,1) + t_m(end,2);
    nmth_target = 12*date_today(1) + date_today(2);
    nfill_m = nmth_target - nmth_m;

    % Fill data
    t_m = [t_m; nan(nfill_m,2)];
    t_m = common_complete_dates_fct(t_m,nfill_m);
    data_m = [data_m; nan(nfill_m,size(data_m,2))];

else
    
    % Check that we don't have "excess" data. For example, observations for
    % April 2024 while we are in March 2024 - otherwise it might indicate
    % an issue
    if ~isempty(end_t_q)
        if end_t_q > size(t_q,1) || end_t_m > size(t_m,1)
            error("Presence of (monthly or quarterly) observations after the intended today's date. Please check the inputs and, if running an evaluation, parameters Eval.data_update_lastyear and Eval.data_update_lastmonth.")
        end
    else
        if end_t_m > size(t_m,1)
            error("Presence of monthly observations after the intended today's date. Please check the inputs and, if running an evaluation, parameters Eval.data_update_lastyear and Eval.data_update_lastmonth.")
        end
    end

end


%% Merge monthly and quarterly data

% Ensure that t_m goes further than t_q (otherwise create additional NaN)
nmth_m = 12*t_m(end,1)+t_m(end,2);
nmth_q = 12*t_q(end,1)+t_q(end,2);
nfill_m = nmth_q - nmth_m;
if nfill_m > 0
    t_m = [t_m; nan(nfill_m,2)];
    t_m = common_complete_dates_fct(t_m,nfill_m);
    data_m = [data_m; nan(nfill_m,size(data_m,2))];
end

% Get dimensions
[mT,nM] = size(data_m);
[~,nQ] = size(data_q);

% Fill in missing values for quarterly variables
temp = NaN(mT,nQ);
for i=1:size(data_q,1)           
    inx = find(t_m(:,1)==t_q(i,1) & t_m(:,2)==t_q(i,2));
    temp(inx,:) = data_q(i,:);
end

data = [data_m temp];   % Final transformed dataset


%% Define starting date of dataset

start = find(t_m(:,1)==Par.startyear & t_m(:,2)==Par.startmonth);
if isempty(start)
    error(['Starting date (Par.startyear and Par.startmonth) not found in the dataset. Please set it no earlier than ',num2str(t_m(1,1)),'M',num2str(t_m(1,2)),'. This constraint is due to the Excel dataset, if feasible you can also consider adding more observations in it.'])
end
data = data(start:end,:);
t_m = t_m(start:end,:);


%% Finalise settings
if do_loop == 0 || do_loop == 1
    if Par_in.block_factors == 1
        % if block factors are included use:
        Par.blocks = blocks;
        Par.r = ones(size(Par.blocks,2),1); % Number of block factors
        Par.r(1) = Par_in.r;    
    else
        Par.blocks = ones(size(data,2),1);
        Par.r = ones(1,1); 
        Par.r(1) = Par_in.r;
    end
end

if do_loop == 2
    Par.blocks_full = blocks;
    Par.size_data_in = size(data,2);
end


xest = [data; NaN(m,nM + nQ)];                  % data file, forecast horizon appended
datet = [t_m; NaN(m,2)];                        % include dates
datet = common_complete_dates_fct(datet,m);     % complete dates series
Par.nM = nM;
Par.nQ = nQ;    

% Find index of quarter of most recent GDP data
gdp_index = find(~isnan(seriesq(:,end)),1,'last') - 1; % subtract 1 bc. of growth rates
Par.yearq = t_q(gdp_index,1);
Par.qend = t_q(gdp_index,2)/3;        % Par.qend gives the quarter for which latest GDP is available

% Save transformations (for bridge equations)
Par.trf.transf_q = transf_q;
Par.trf.transf_m = transf_m;

% Load list of models to loop over if required by user
if do_loop == 2

    % Read list of models
    opts = detectImportOptions(Loop.list_name, 'Sheet', 'Model list');
    opts = setvartype(opts, "batch", 'char');
    Loop.mdl_list = readtable(Loop.list_name,opts,'Sheet', 'Model list'); % read list of models
    disp('Section 3: List of models to check loaded')

    % Append list with models using alternative covid corrections
    if Loop.alter_covid == 1 
        h = height(Loop.mdl_list);
        for i=1:h
            if Loop.mdl_list{i,"do_Covid"} == 0
                cov_cor = [1 2 3];
            elseif Loop.mdl_list{i,"do_Covid"} == 1
                cov_cor = [0 2 3];
            elseif Loop.mdl_list{i,"do_Covid"} == 2
                cov_cor = [0 1 3];
            elseif Loop.mdl_list{i,"do_Covid"} == 3
                cov_cor = [0 1 2];
            end
            
            % Append list of models
            for j=1:length(cov_cor)                         
                new_T = Loop.mdl_list(i,:);
                new_T{1,"do_Covid"} = cov_cor(j);
                new_T(1,"batch") = {strcat(Loop.mdl_list{i,"batch"},'_','covcor_',num2str(cov_cor(j)))};
                Loop.mdl_list = [Loop.mdl_list;new_T];
            end
        end
    end % end of condition on Loop.alter_covid == 1
end % end of condition on do_loop==2

end