function[] = common_CheckDataAvailability(fullnames)
% [] = CheckDataAvailability(fullnames), for a cell array of strings
% 'fullnames', returns a warning if a series name contains either DISC or
% SUSP. This is the usual Haver procedure to indicate that a series is
% discontinued or suspended.
%
% Code by S. Delle Chiaie and F. Kurcz
%

if iscell(fullnames)
    % matlab recommends use of 'contains' but it is only available for
    % MATLAB R2018
    indx1 = ~cellfun('isempty',strfind(fullnames, 'DISC'));
    indx2 = ~cellfun('isempty',strfind(fullnames, 'SUSP'));
    indx = indx1 | indx2;
    indx_abs = find(indx);
    for j = 1:sum(indx)
        warning('%s is suspended / discontinued.',char(fullnames(indx_abs(j))))
    end
end