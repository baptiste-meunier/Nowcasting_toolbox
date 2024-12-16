function [] = common_save_excel(chartdata,quarterahead_chartdata,Par,excel_outputfile,newsfile,prev_news,groups_name,contrib,names_contrib,info_indiv_contrib,matrix_indiv_contrib,nameseries,fullnames,heatmap,datet,country,Res,groups,do_range,range,indxGDP)
% This scripts saves the output in the excel file specified in "excel_outputfile". 
%
% Adapted from an initial version by S. Delle Chiaie and F. Kurcz
%
% INPUTS 
% - chartdata [cell array] = data to be pasted on Excel for current quarter 
% - quarterahead_chartdata [x] = data to be pasted on Excel for next quarter
% - Par [structure] = parameters of DFM
%       o blocks [vector/matrix] = block identification
%       o r [scalar] = number of estimated factors
%       o block_factors [scalar] = switch on whether to do block factors (=1) or not (=0)
%       o p [scalar] = number of lags
%       o nM [scalar] = number of monthly variables
%       o nQ [scalar] = number of quarterly variables
%       o idio [scalar] = idiosyncratic component specification: 1 = AR(1), 0 = iid
%       o yearq [scalar] = year for which latest GDP is available
%       o qend [scalar] = quarter for which latest GDP is available
%       o thresh [scalar] = threshold for convergence of EM algorithm
%       o max_iter [scalar] = number of iterations for the EM algorithm
% - filename [string] = Excel file containing tracking and news decomposition
% - newsfile [string] = name of mat file to compare with for news decomposition
% - prev_news [logical] = whether the previous comparison is available
% - groups_name [cell array] = names of groups
% - contrib [matrix] = contributions of groups and individual variables to nowcasts
% - names_contrib [string array] = names of group and individual contributions
% - info_indiv_contrib [cell array] = full names and groups ID / names for the individual series
% - matrix_indiv_contrib [cell array] = matrix with contributions allocated to groups (for the Excel charts)
% - nameseries [cell array] = names of the series
% - fullnames [cell array] = full names of the series
% - heatmap [stucture]
%       o zscores [matrix] = z-scores for the actual data (re-ordered by group)
%       o names [cell array] = tables with name of groups and variables (re-ordered by group)
% - datet [matrix] = dates (year/month)
% - country [structure] = country parameters
%       o name [string] = name of the country
%       o stat [string] = transformation of target data (qoq or qoqAR)
%       o model [string] = model type
% - Res [structure] = results of DFM
%       o X_sm [matrix] = Kalman-smoothed dataset (NaN replaced with expectations)
%       o F [matrix] = smoothed states 
%                      Rows give time
%                      Columns are organized according to Res.C
%       o C [matrix] = observation matrix
%                      Rows correspond to each series
%                      Columns give factor loadings
%                      For example, col. 1-5 give loadings for the first block and are organized in reverse-chronological order (f^G_t, f^G_t-1, f^G_t-2, f^G_t-3, f^G_t-4)
%       o R [matrix] = covariance for observation matrix residuals
%       o A [matrix] = transition matrix
%                      Square matrix that follows the same organization scheme as Res.C's columns
%                      Identity matrices are used to account for matching terms on the left and right-hand side.
%                      For example, we place an I4 matrix to account for matching (f_t-1; f_t-2; f_t-3; f_t-4) terms.
%       o Q [matrix] = covariance for transition equation residuals
%       o Mx [vector] = series mean
%       o Wx [vector] = series standard deviation
%       o Z_0 [vector] = initial value of state
%       o V_0 [matrix] = initial value of covariance matrix
%       o r [scalar] = number of estimated factors
%       o p [scalar] = number of lags in transition equation
%       o L [scalar] = log-likelihood
%       o groups [vector] = indicator for group belonging
%       o series [cell vector] = mnemotechnics of series in X_sm (generally their Haver handle)
%       o name_descriptor [cell vector] = full names of series in X_sm
% - indxGDP [scalar] = index of the last available value for the target
%

% Checking that Excel output file exists
if ~(exist(excel_outputfile,'file') == 2)
    error('The specified excel_outputfile does not exist. Check the name and run the Mainfile again.')
end


%% Prepare parameters

% Get dates for the nowcast / forecast
Qend = Par.qend;        % quarter for which we have latest GDP data
Yearq = Par.yearq;      % year for which we have latest GDP data

% Locate nowcast and forecast quarter and year
QuartersQ = {'Q1','Q2','Q3','Q4','Q1','Q2'};
Qsheet = QuartersQ{Qend+1};         % quarter we are nowcasting
Qsheet_fore = QuartersQ{Qend + 2};  % quarter we are forecasting
if Qend == 4 % advancing year if needed
    Yearq = Yearq + 1;
    Yearq_fore = Yearq; 
elseif Qend == 3
    Yearq_fore = Yearq + 1; 
else
    Yearq_fore = Yearq;
end

% Names of headers (generic)
header_names = ['Date','Current nowcast','Variable coverage', 'FDA (last 10y)', ...
                'MAE','Prediction - MAE','Prediction + MAE', ...
                ' ', ...
                'Date (old nowcast)','Old nowcast',' ', ...
                string(groups_name),'Revisions (data and model)' ...
                ' ', ...
                string(nameseries)];


%% Write for current quarter
if strcmp('cur_nowcast.mat',newsfile)

    SheetName = strcat(num2str(Yearq),Qsheet);  % sheet name for quarter we are nowcasting
    SheetName_fore = strcat(num2str(Yearq_fore),Qsheet_fore);   % sheet name for quarter we are forecasting

    % Try opening the sheet for current quarter
    % Otherwise we would need a new
    try
        [num,~] = xlsread(excel_outputfile,SheetName);
        is_sheetname = 1;
    catch
        num = [];
        is_sheetname = 0;
    end

else

    SheetName = 'SpecialRun';
    [num,~] = xlsread(excel_outputfile,SheetName);
    is_sheetname = 1; % this one should always exist

end

Row = num2str(size(num,1) + 2); % +2 to add one line and take into account the fact that the first line (text) is not read here

% Export data:
Ranger1 = ['A',Row,':K',Row];                   
Ranger2 = ['L',Row];

% Specify in output that this is a new model if that's the case
if (prev_news == false) && (str2double(Row) > 2) % meaning there are already some data in the sheet but you try to put a new model or file
    chartdata(end,8) = {'Model/file changed'};
    new_model = 1;
else
    new_model = 0;
end

try

    % Write predictions and date of predictions
    xlswrite(excel_outputfile, chartdata(1,1:11), ...
             SheetName, Ranger1)

    % Write news decomposition
    if new_model == 1 % if new model, write new columns names

        writecell([groups_name,'Revisions (data and model)',' ',nameseries],excel_outputfile, ...
                  'Sheet',SheetName, ...
                  'Range',Ranger2) % write new headers
        xlswrite(excel_outputfile,[' ',' ']', ...
                SheetName,['I',Row,':J',Row]) % replacing the 'old' nowcasts which do not make sense

    elseif ~is_sheetname % if first sheet, write column headers

        xlswrite(excel_outputfile,header_names, ...
                SheetName,'A1')
        xlswrite(excel_outputfile,[' ',' ']', ...
                SheetName,'I2:J2') % replacing the 'old' nowcasts which do not make sense

    else % otherwise, write contributions

        xlswrite(excel_outputfile,cell2mat(chartdata(1,12:end)), ...
                 SheetName,Ranger2)

    end

catch exception % catch so the program doesn't break down because of a locked Excel file

    if strcmp(exception.identifier,'MATLAB:xlswrite:LockedFile') % generally the issue is that the file is locked by the user
        fprintf(['The excel sheet ''', excel_outputfile, ''' is locked. Close the file, return to the Command Window by clicking into it and press any key to continue. \n'])

        pause

        % Do the same as above
        xlswrite(excel_outputfile,chartdata(1,1:11), ...
                 SheetName,Ranger1) 

        % Write news decomposition
        if new_model == 1

            writecell([groups_name,'Revisions (data and model)',' ',nameseries],excel_outputfile, ...
                  'Sheet',SheetName, ...
                  'Range',Ranger2) % write new headers 
            xlswrite(excel_outputfile,[' ',' ']', ...
                SheetName,['I',Row,':J',Row]) % replacing the 'old' nowcasts which do not make sense

        elseif ~is_sheetname % if first sheet, write column headers

            xlswrite(excel_outputfile,header_names, ...
                    SheetName,'A1')
            xlswrite(excel_outputfile,[' ',' ']', ...
                    SheetName,'I2:J2') % replacing the 'old' nowcasts which do not make sense

        else

            xlswrite(excel_outputfile,cell2mat(chartdata(1,12:end)), ...
                     SheetName,Ranger2) 
        end

        clear exception  

    end 
end


%% Write for following quarter
% NB: it does not write for the first (because there is no news yet)
if strcmp('cur_nowcast.mat',newsfile)
    try

        % Try to read in row (will work if the sheet exists)
        [numnew,~] = xlsread(excel_outputfile,SheetName_fore);
        Rownew = num2str(size(numnew,1) + 2);
        Rangernew1 = ['A',Rownew,':K',Rownew]; 
        Rangernew2 = ['L',Rownew]; 

        % Specify in output that this is a new model if that's the case
        if prev_news == false && str2double(Rownew) > 2 % meaning there are already some data in the sheet but you try to put a new model
            quarterahead_chartdata(end,9) = {'Model/file changed'};
            new_model = 1;
        else
            new_model = 0;
        end
        
        % Write the last line
        xlswrite(excel_outputfile,quarterahead_chartdata(1,1:11), ...
                 SheetName_fore,Rangernew1);
        if new_model == 1
            writecell([groups_name,'Revisions (data and model)',' ',nameseries],excel_outputfile, ...
              'Sheet',SheetName_fore, ...
              'Range',Rangernew2) % write new headers
        else
            xlswrite(excel_outputfile,cell2mat(quarterahead_chartdata(1,12:end)), ...
                     SheetName_fore,Rangernew2);
        end
    
    catch % otherwise it means this is the first we are predicting this quarter
      
        % Create a new sheet with headers ...
        xlswrite(excel_outputfile,header_names, ...
            SheetName_fore,'A1')
        
        % ... and write in second row
        % NB1: not the contributions because this is a new quarter
        % NB2: and not the 'old' nowcast because this is the first time it is (officially) done
        xlswrite(excel_outputfile,quarterahead_chartdata(end,1:7), ...
                 SheetName_fore,'A2:G2')

    end

end 


%% Write contributions 

% Value of contributions
xlswrite(excel_outputfile,contrib, ...
         'Contributions','B2');

% Headers of contributions
header_contrib = {'',strcat(num2str(Yearq),Qsheet),strcat(num2str(Yearq_fore),Qsheet_fore)};
xlswrite(excel_outputfile,header_contrib, ...
         'Contributions','A1');

% Names of contributions
xlswrite(excel_outputfile,names_contrib, ...
         'Contributions','A2');

% Write information on individual series
xlswrite(excel_outputfile,info_indiv_contrib,...
    'Contributions','D2');

% Write matrix of individual contributions
xlswrite(excel_outputfile,matrix_indiv_contrib,...
    'Contributions','G2');


%% Write range of nowcasts

if do_range == 1

    % Values of range
    range_min = [min(range.out(:,1)),min(range.out(:,2))];
    range_max = [max(range.out(:,1)),max(range.out(:,2))];
    values_range = vertcat(range.out,[NaN,NaN],range_min,range_max,contrib(1,:));
    xlswrite(excel_outputfile,values_range, ...
             'Range','B2');
    
    % Header of range
    header_range = {'',strcat(num2str(Yearq),Qsheet),strcat(num2str(Yearq_fore),Qsheet_fore)};
    xlswrite(excel_outputfile,header_range, ...
             'Range','A1');
    
    % Names
    names_range = vertcat(string(range.list(:,1)),NaN,'Min','Max','Full model');
    xlswrite(excel_outputfile,names_range, ...
             'Range','A2');

end


%% Write the heatmaps

% Groups and full names of individual series
xlswrite(excel_outputfile,heatmap.names, ...
         'Heatmap','A3');

% Groups (for aggregate heatmap)
num_cell = 3 + size(heatmap.names,1) + 4;
str_cell = strcat('A',num2str(num_cell));
xlswrite(excel_outputfile,heatmap.names_agg, ...
         'Heatmap',str_cell);

% Dates (last 18 months)
datet_18 = datet(indxGDP+6-17:indxGDP+6,:);
start_18 = datetime(datet_18(1,1),datet_18(1,2),15);
end_18 = datetime(datet_18(end,1),datet_18(end,2),15);
datetime_18 = start_18:calmonths(1):end_18;

% Dates for the heatmap of individual series
xlswrite(excel_outputfile,exceltime(datetime_18), ...
         'Heatmap','C2');

% Dates for the aggregate heatmap
num_cell = 3 + size(heatmap.names,1) + 3; % +3 only (intead of +4) because these are the headers
str_cell = strcat('C',num2str(num_cell));
xlswrite(excel_outputfile,exceltime(datetime_18), ...
         'Heatmap',str_cell);

% Write heatmap of individual series
heatmap_18_init = heatmap.zscores(indxGDP+6-17:indxGDP+6,:)';
heatmap_18 = num2cell(heatmap_18_init);
heatmap_18(isnan(heatmap_18_init)) ={'NaN'};
xlswrite(excel_outputfile,heatmap_18, ...
         'Heatmap','C3');

% Write aggregate heatmap
heatmap_agg_18_init = heatmap.zscores_agg(indxGDP+6-17:indxGDP+6,:)';
heatmap_agg_18 = num2cell(heatmap_agg_18_init);
heatmap_agg_18(isnan(heatmap_agg_18_init)) ={'NaN'};
num_cell = 3 + size(heatmap.names,1) + 4;
str_cell = strcat('C',num2str(num_cell));
xlswrite(excel_outputfile,heatmap_agg_18, ...
         'Heatmap',str_cell);


end % End of function