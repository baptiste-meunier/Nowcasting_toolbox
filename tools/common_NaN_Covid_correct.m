function[x_out,nM_out,blocks_out,r_out,groups_out,groups_name_out,nameseries_out,fullnames_out,transf_m_out,transf_q_out] = common_NaN_Covid_correct(x_in,datet_in,do_Covid,nM_init,blocks_init,r_init,groups_init,groups_name_init,nameseries_init,fullnames_init,transf_m_init,transf_q_init)
% This script corrects the dataset for Covid observations and for NaN
% Only series with at least 3 observations are kept
% Also adjust blocks (because they can't be launched only on quarterly variables)
% Subsequent model parameters (nM,nQ,r,blocks,groups) are changed
% 
% INPUTS 
% - x_in [matrix] = values of the input series
% - datet_in [matrix] = dates (year / month)
% - do_Covid [scalar] = switch on correction method
%                   0 = do nothing
%                   1 = add dummies (one for Q2 2020 and one for Q3 2020)
%                   2 = put to NaN (from Feb. 2020 to Sep. 2020 included)
%                   3 = outlier-correction (with outliers replaces by NaN) = parameters of DFM
%                   4 = add dummies (one for Q1 2020 and one for Q2 2020)
% - nM_init [scalar] = number of monthly variables (in Par structure)
% - blocks_init [matrix] = blocks belonging (in Par structure)
% - r_init [scalar / vector] = number of factors
% - groups_init [vector] = identification of groups
% - groups_name_init [cell vector] = name of the groups
% - namseseries_init [cell array] = name of the series
% - fullnames_init [cell array] = full names of the series
% - transf_m_init [vector] = transformations for monthly variables
% - transf_q_init [vector] = transformations for quarterly variables
%
% OUTPUTS
% - x_out [matrix] = corrected values of the input series
% - nM_out [scalar] = number of corrected monthly variables
% - blocks_out [matrix] = adjusted blocks belonging
% - r_out [scalar / vector] = adjusted number of factors
% - groups_out [vector] = adjusted identification of groups
% - groups_name_out [cell vector] = adjusted name of the groups
% - nameseries_out [cell array] = adjusted name of the series
% - fullnames_out [cell array] = adjusted full names of the series
% - transf_m_out [vector] = adjusted transformations for monthly variables
% - transf_q_out [vector] = adjusted transformations for quarterly variables
%


%% Step 1: correct for NaN and Inf

% Transform 'Inf' into NaN
x_in(isinf(x_in)) = NaN;

% Take only variables with at least two observations
idx_keep = find(sum(~isnan(x_in)) > 2);

% Check that input variables are not there twice
% If that's the case, remove the (first) index from idx_keep
for ii = 1:size(x_in,2)
    col_ii = x_in(:,ii);
    for jj = ii:size(x_in,2)
        if ii ~= jj
            col_jj = x_in(:,jj);
            if isequaln(col_ii,col_jj)
                idx_keep = idx_keep(idx_keep ~= ii);
            end
        end
    end
end

% Adjust input data
x_temp = x_in(:,idx_keep);

% Adjust names
nameseries_temp = nameseries_init(idx_keep);
fullnames_temp = fullnames_init(idx_keep);

% Adjust transformations
idx_keep_m = idx_keep(idx_keep<=nM_init); % monthly variables
transf_m_out = transf_m_init(idx_keep_m);
idx_keep_q = idx_keep(idx_keep>nM_init); % quarterly variables
idx_keep_q = idx_keep_q - nM_init;
transf_q_out = transf_q_init(idx_keep_q);

% Check that the target series has enough data points
if max(idx_keep) < size(x_in,2)
    error('The target variable (at the right of the dataset) has not enough observations (less than 3). Please consider taking another target series or starting later the out-of-sampl evaluation')
end

% Adjust number of monthly series
nM_temp = sum(idx_keep <= nM_init);

% Adjust groups
groups_temp = groups_init(idx_keep); % changing groups
ID_groups = unique(groups_temp(:)); % ID of the groups in x_temp
groups_name_temp = groups_name_init(ID_groups); % subsetting also the names of the groups

% Adjust blocks
blocks_temp = blocks_init(idx_keep',:);
idx_keep = sum(blocks_temp,1) > 0; % to delete blocks that would now be empty
blocks_temp = blocks_temp(:,idx_keep);
r_temp = r_init(idx_keep);


%% Step 2: adjust for Covid observations

switch do_Covid
    
    case 0

        x_out = x_temp;
        nM_out = nM_temp;
        blocks_out = blocks_temp;
        r_out = r_temp;
        groups_out = groups_temp;
        groups_name_out = groups_name_temp;
        nameseries_out = nameseries_temp;
        fullnames_out = fullnames_temp;

    case 1

        % Create dummies
        ind_2020q2 = find(datet_in(:,1)==2020 & datet_in(:,2)==6);
        d2020_q2 = zeros(size(x_temp,1),1);
        d2020_q2(ind_2020q2,1) = 1;
        ind_2020q3 = find(datet_in(:,1)==2020 & datet_in(:,2)==9);
        d2020_q3 = zeros(size(x_temp,1),1);
        d2020_q3(ind_2020q3,1) = 1;

        % Check that there are at least 5 other observations at the date that the dummy takes value 1
        % Otherwise this would cause the DFM to not converge
        idx_NotAllNan = find(sum(~isnan(x_temp),2)>4);

        % Add dummies to dataset (unless the dummy is all 0)
        nM_out = nM_temp;
        x_out = x_temp;
        if sum(d2020_q2) > 0 && ismember(ind_2020q2,idx_NotAllNan)
            x_out = horzcat(x_out(:,1:nM_out),d2020_q2,x_out(:,(nM_out+1):end));
            nameseries_temp = horzcat(nameseries_temp(1:nM_out),{'Covid_20q2'},nameseries_temp(nM_out+1:end));
            fullnames_temp = horzcat(fullnames_temp(1:nM_out),{'Covid_20q2'},fullnames_temp(nM_out+1:end));
            nM_out = nM_out + 1;
        end
        if sum(d2020_q3) > 0 && ismember(ind_2020q3,idx_NotAllNan)
            x_out = horzcat(x_out(:,1:nM_out),d2020_q3,x_out(:,(nM_out+1):end));
            nameseries_temp = horzcat(nameseries_temp(1:nM_out),{'Covid_20q3'},nameseries_temp(nM_out+1:end));
            fullnames_temp = horzcat(fullnames_temp(1:nM_out),{'Covid_20q3'},fullnames_temp(nM_out+1:end));
            nM_out = nM_out + 1;
        end
        nM_add = nM_out - nM_temp;

        % Adjust the rest of the parameters if dummies have been added
        if nM_add > 0    

            % Adjust r and blocks (if block_factors == 1)
            if size(r_temp,1) > 1
                blocks_out = vertcat(blocks_temp(1:nM_temp,:), ... % all (monthly) variables
                                 zeros(nM_add,size(blocks_temp,2)), ... % inserting lines for new dummies
                                 blocks_temp((nM_temp+1):end,:));
                blocks_out(nM_temp+1:nM_temp+nM_add,1) = 1; % putting 1 to add dummies to the global factor
                temp_col = vertcat(zeros(nM_temp,1), ...
                                   ones(nM_add,1), ...
                                   zeros(size(blocks_temp,1)-nM_temp,1)); % additional column with the "Covid" factor
                blocks_out = horzcat(blocks_out,temp_col);
                r_out = [r_temp;1];
            else
                blocks_out = vertcat(blocks_temp,ones(nM_add,1)); 
                r_out = r_temp;
            end
    
            % Adjust groups and transformations
            id_max = max(groups_temp);
            groups_out = horzcat(groups_temp(:,1:nM_temp), ... % all (monthly) variables before the dummies
                                 repmat(id_max+1,1,nM_add), ... % inserting new group for dummies
                                 groups_temp(:,(nM_temp+1):end));
            groups_name_out = horzcat(groups_name_temp,{'Covid'});
            transf_m_out = [transf_m_out,zeros(1,nM_add)];

        else

            blocks_out = blocks_temp;
            r_out = r_temp;
            groups_out = groups_temp;
            groups_name_out = groups_name_temp;

        end

        nameseries_out = nameseries_temp;
        fullnames_out = fullnames_temp;

    case 2

        ind_start = find(datet_in(:,1)==2020 & datet_in(:,2)==2);
        ind_end = find(datet_in(:,1)==2020 & datet_in(:,2)==9);
        if ~isempty(ind_start) && isempty(ind_end) % case when Feb. 2020 is in x_temp, but not Sept. 2020
            ind_end = size(x_temp,1); % then NaN should go until the end (which is before Sep. 2020)
        end
        x_out = x_temp;
        x_out(ind_start:ind_end,:) = NaN;
        nM_out = nM_temp;
        blocks_out = blocks_temp;
        r_out = r_temp;
        groups_out = groups_temp;
        groups_name_out = groups_name_temp;
        nameseries_out = nameseries_temp;
        fullnames_out = fullnames_temp;

    case 3

        [x_out,~] = common_outliers(x_temp,0);
        nM_out = nM_temp;
        blocks_out = blocks_temp;
        r_out = r_temp;
        groups_out = groups_temp;
        groups_name_out = groups_name_temp;
        nameseries_out = nameseries_temp;
        fullnames_out = fullnames_temp;

    case 4

        % Create dummies
        ind_2020q1 = find(datet_in(:,1)==2020 & datet_in(:,2)==3);
        d2020_q1 = zeros(size(x_temp,1),1);
        d2020_q1(ind_2020q1,1) = 1;
        ind_2020q2 = find(datet_in(:,1)==2020 & datet_in(:,2)==6);
        d2020_q2 = zeros(size(x_temp,1),1);
        d2020_q2(ind_2020q2,1) = 1;

        % Check that there are at least 5 other observations at 1 of dummy
        % Otherwise this would cause the DFM to not converge
        idx_NotAllNan = find(sum(~isnan(x_temp),2)>4);

        % Add dummies to dataset (unless the dummy is all 0)
        nM_out = nM_temp;
        x_out = x_temp;
        if sum(d2020_q1) > 0 && ismember(ind_2020q1,idx_NotAllNan)
            x_out = horzcat(x_out(:,1:nM_out),d2020_q1,x_out(:,(nM_out+1):end));
            nameseries_temp = horzcat(nameseries_temp(1:nM_out),{'Covid_20q1'},nameseries_temp(nM_out+1:end));
            fullnames_temp = horzcat(fullnames_temp(1:nM_out),{'Covid_20q1'},fullnames_temp(nM_out+1:end));
            nM_out = nM_out + 1;
        end
        if sum(d2020_q2) > 0 && ismember(ind_2020q2,idx_NotAllNan)
            x_out = horzcat(x_out(:,1:nM_out),d2020_q2,x_out(:,(nM_out+1):end));
            nameseries_temp = horzcat(nameseries_temp(1:nM_out),{'Covid_20q2'},nameseries_temp(nM_out+1:end));
            fullnames_temp = horzcat(fullnames_temp(1:nM_out),{'Covid_20q2'},fullnames_temp(nM_out+1:end));
            nM_out = nM_out + 1;
        end
        nM_add = nM_out - nM_temp;

        % Adjust the rest of the parameters if dummies have been added
        if nM_add > 0    

            % Adjust r and blocks (if block_factors == 1)
            if size(r_temp,1) > 1
                blocks_out = vertcat(blocks_temp(1:nM_temp,:), ... % all (monthly) variables
                                 zeros(nM_add,size(blocks_temp,2)), ... % inserting lines for new dummies
                                 blocks_temp((nM_temp+1):end,:));
                blocks_out(nM_temp+1:nM_temp+nM_add,1) = 1; % putting 1 to add dummies to the global factor
                temp_col = vertcat(zeros(nM_temp,1), ...
                                   ones(nM_add,1), ...
                                   zeros(size(blocks_temp,1)-nM_temp,1)); % additional column with the "Covid" factor
                blocks_out = horzcat(blocks_out,temp_col);
                r_out = [r_temp;1];
            else
                blocks_out = vertcat(blocks_temp,ones(nM_add,1)); 
                r_out = r_temp;
            end
    
            % Adjust groups and transformations
            id_max = max(groups_temp);
            groups_out = horzcat(groups_temp(:,1:nM_temp), ... % all (monthly) variables before the dummies
                                 repmat(id_max+1,1,nM_add), ... % inserting new group for dummies
                                 groups_temp(:,(nM_temp+1):end));
            groups_name_out = horzcat(groups_name_temp,{'Covid'});
            transf_m_out = [transf_m_out,zeros(1,nM_add)];

        else

            blocks_out = blocks_temp;
            r_out = r_temp;
            groups_out = groups_temp;
            groups_name_out = groups_name_temp;

        end

        nameseries_out = nameseries_temp;
        fullnames_out = fullnames_temp;

end


%% Step 3: Check that blocks do not contain only quarterly variables
% Otherwise blocks are discarded and variables launch only on the global factor

% Control blocks
idx_keep = sum(blocks_out(1:nM_out,:),1) > 0;

% Adjust blocks and factors
blocks_out = blocks_out(:,idx_keep);
r_out = r_out(idx_keep);

end % end of main function