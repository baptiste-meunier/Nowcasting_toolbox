function [X_out,date_out] = BVAR_filldata(X_in,date_in)
% This script extends the data (X) and dates (date) to ensure they start at
% 1st month of a quarter and end at last month of a quarter. Data is filled
% with NaN
%
% INPUTS
% - X_in [matrix] = data
% - date_in [matrix] = monthly dates (year/month). Quarters 1, 2, 3, 4 are indicated by  3, 6, 9, 12 
%
% OUTPUTS
% - X_out [matrix] = extended data
% - date_out [matrix] = extended monthly dates (year/month).
%

% Check start date
if mod(date_in(1,2),3)~=1

    if mod(date_in(1,2),3)==0 % case when first date is 3rd month of quarter
        n_nan = 2;
        date_out = [date_in(1,1),date_in(1,2)-2; ...
                    date_in(1,1),date_in(1,2)-1; ...
                    date_in];
    elseif mod(date_in(1,2),3)==2 % case when first date is 2nd month of quarter
        n_nan = 1;
        date_out = [date_in(1,1),date_in(1,2)-1; ...
                    date_in];
    end
    X_out = [NaN(n_nan,size(X_in,2)); ...
             X_in];
    
else

    X_out = X_in;
    date_out = date_in;

end

% Check end date
if mod(date_out(end,2),3)~=0

    if mod(date_out(end,2),3)==1 % case when last date is 1st month of quarter
        n_nan = 2;
        date_out = [date_out; ...
                    date_out(end,1),date_out(end,2)+1; ...
                    date_out(end,1),date_out(end,2)+2];
    elseif mod(date_out(end,2),3)==2 % case when first date is 2nd month of quarter
        n_nan = 1;
        date_out = [date_out; ...
                    date_out(end,1),date_out(end,2)+1];
    end
    X_out = [X_out; ...
             NaN(n_nan,size(X_out,2));];
    
end

end