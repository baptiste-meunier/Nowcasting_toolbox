function [range] = common_range(xest_out,Par,datet,groups_name,groups,country,nameseries)
% This function creates a range of nowcasts by restimating the model,
% removing one group at a time, then removing two groups at a time. Thus,
% the number of estimated models depends on the number of groups. 
% Function is activated by user via do_range switch.
%
% INPUTS 
% - xest_out [matrix] = dataset
% - Par [structure] = parameters of the model
% - datet [matrix] = dates (year / month)
% - groups_name [cell vector] = name of the groups
% - groups [vector] = identification of groups
% - country [structure] = information on model and country used
% - nameseries [cell vector] = ID of the series in Haver
%
% OUTPUTS
% - range
%


%% 1 Create list to loop over

range.list = {};
c = 0;
groups_input = groups(1:(end-1));
for i = 1:length(groups_name) % Create list to loop over
    if ismember(i,groups_input)
       range.list{c+1,1} = strcat("No ",groups_name(i)); % name
       range.list{c+1,2} = [~(groups_input==i),1]; % selection of vars
       range.list{c+1,3} = sum(range.list{c+1,2}(1:Par.nM)); % number of m vars
       range.list{c+1,4} = sum(range.list{c+1,2}((Par.nM+1):end)); % number of q vars
       c = c + 1;
    end
end
combs = nchoosek(1:length(groups_name),2);
for j = 1:size(combs,1)
    if and(ismember(combs(j,1),groups_input), ismember(combs(j,2),groups_input))
        range.list{c+1,1} = strcat("No ",groups_name(combs(j,1))," and ",groups_name(combs(j,2))); % name
        range.list{c+1,2} = [~or(groups_input==combs(j,1),groups_input==combs(j,2)),1]; % selection of vars
        range.list{c+1,3} = sum(range.list{c+1,2}(1:Par.nM)); % number of m vars
        range.list{c+1,4} = sum(range.list{c+1,2}((Par.nM+1):end)); % number of q vars
        c = c + 1;
    end
end

%% 2 Loop over RON_list

xest_main = xest_out; % backup parameters and data of main nowcast
Par_main = Par;
range.out = zeros(height(range.list),2);

for ii = 1:height(range.list)

    disp(['Estimating range ',num2str(ii),'/',num2str(height(range.list))]);

    %% 2.1 Cut data set

    range.data = xest_main(:,find(range.list{ii,2}));
    xest_temp = range.data;

    %% 2.2 Adjust Par structure

    Par.nM = range.list{ii,3};
    Par.nQ = range.list{ii,4};
    Par.blocks = Par_main.blocks(find(range.list{ii,2}),:);
    Par.r = Par_main.r;

    sel_var = find(range.list{ii,2});
    sel_var_m = sel_var(sel_var<=Par_main.nM);
    Par.trf.transf_m = Par_main.trf.transf_m(sel_var_m);
    sel_var_q = sel_var(sel_var>Par_main.nM);
    sel_var_q = sel_var_q - Par_main.nM;
    Par.trf.transf_q = Par_main.trf.transf_q(sel_var_q);

    % Check blocks (must contain at least one monthly variable)
    idx_keep = sum(Par.blocks(1:Par.nM,:),1) > 0;
    Par.blocks = Par.blocks(:,idx_keep);
    Par.r = Par.r(idx_keep);

    %% 2.3 Pass to estimation
    
    switch country.model
            case 'DFM'
                Res_range = DFM_estimate(xest_temp,Par);
            case 'BEQ'
                re_estimate = 1;
                coeffs_in = [];
                Res_range = BEQ_estimate(xest_temp,Par,datet,nameseries,re_estimate,coeffs_in);
            case 'BVAR'
                Res_range = BVAR_estimate(xest_temp,Par,datet);
    end

    %% 3 Write results
    
    % Extract now- and forecast
    track_range = [datet Res_range.X_sm(:,end)];
    indx_range = find(~isnan(xest_main(:,end)),1,'last');
    now_range = track_range(indx_range + 3,3);
    fore_range = track_range(indx_range + 6,3);
        
    % Write 
    range.out(ii,1) = now_range;
    range.out(ii,2) = fore_range;

end

end