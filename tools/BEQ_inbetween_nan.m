function X_out = BEQ_inbetween_nan(X_in)
% This scripts runs a linear interpolation of in-between NaN
%
% INPUT
% - X_in [matrix] = input data with observations in row and series in columns
%
% OUTPUT
% - X_out [matrix] = data with in-between NaN interpolated
%

X_out = fillmissing(X_in,'linear', ... % states linear interpolation
                    1, ... % fills by columns (= series)
                    'EndValues','none'); % ensures trailing NaN are not filled

end