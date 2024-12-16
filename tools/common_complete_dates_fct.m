function[datet] = common_complete_dates_fct(datet,m)
% This script completes date
% NB: datet should have NaN for dates to complete
%
% INPUTS
% - datet [matrix] = dates (year / month) with NaN for dates to fill
% - m [scalar] = number of months ahead (should correspond to the NaN)
% 
% OUTPUT
% - datet [matrix] = dates (year / month) filled
%

for i=1:m
    if datet(end-m+i-1,2)==12
        datet(end-m+i,2) = 1;
        datet(end-m+i,1) = datet(end-m+i-1,1)+1;
    else
        datet(end-m+i,2) = datet(end-m+i-1,2)+1;
        datet(end-m+i,1) = datet(end-m+i-1,1);
    end
end