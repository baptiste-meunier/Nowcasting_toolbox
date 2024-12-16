function X_new = BEQ_Run_extrapolation_BVAR(X_in,Par_BVAR)
% This script runs the extrapolation of data using BVAR
%
% Code adapted from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
% Changes relative to Banbura et al. (2023) are
%     - no transformation and trimming
%     - i2 (index for last non-NaN) is series-specific
%     - give average of nDraw as X_new
%
% INPUTS
% - X_in [matrix] = input data
% - Par_BVAR [structure] = parameters of the BVAR
%       o nDraw [scalar] = number of draws (0 = point forecast)
%       o KK [scalar] = xxx
%       o iRW [scalar] = xxx
%       o lags [scalar] = number of lags
%       o lambda [scalar] = overall shrinkage parameters
%
% OUTPUTS
% - X_new [matrix] = interpolated data
%

% Get position of last non-NaN value for each column
i2 = NaN(1,size(X_in,2));
for nn = 1:size(X_in,2)
    i2(nn) = find(~isnan(X_in(:,nn)), 1, 'last');
end

% Run forecast
X_new = BEQ_BVAR_fcst(X_in,Par_BVAR);

% Replacing values with actual data
fcstsize3 = max(1,Par_BVAR.nDraw);
for nn = 1:size(X_in,2)
    X_new(1:i2(nn),nn,:) = repmat(X_in(1:i2(nn),nn),[1 1 fcstsize3]);
end

% Averaging over nDraw
if fcstsize3 > 1
    X_new = mean(X_new,3,"omitnan");
end

end