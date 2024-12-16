function [parameters,parameters_names,results_names] = common_initiate_loop(n_iter,country_model)
% This script creates containers for evaluation of models wihtin a loop
%
% INPUTS
% - n_iter [scalar] = number of random models to be tested
% - country_model [string] = model type
%
% OUTPUTS
% - parameters [matrix] = container for hyper-parameters
% - parameters_name [cell array] = names of parameters
% - results_names [cell array] = names of columns in results file
% - n_col_results [scalar] = number of columns in results file
%

switch country_model
    case 'DFM'
        parameters = nan(n_iter,7);
        parameters_names = {'n iteration','do_Covid','start year (est.)','nb variables','nb factors','nb lags','blocks',' variables','groups'};
    case 'BEQ'
        parameters = nan(n_iter,7);
        parameters_names = {'n iteration','do_Covid','start year (est.)','nb variables','nb lagM','nb lagQ','nb lagY','variables','groups'};
    case 'BVAR'
        parameters = nan(n_iter,5);
        parameters_names = {'n iteration','do_Covid','start year (est.)','nb variables','nb lags','variables','groups'};
end

other_names = {'RMSE all','RMSE pre-Covid (<2020)','RMSE Covid (2020)','RMSE post-Covid (>2020)','RMSE no Covid', ...
               'FDA all','FDA pre-Covid (<2020)','FDA Covid (2020)','FDA post-Covid (>2020)','FDA no Covid', ...
               '1m - RMSE all','1m - RMSE pre-Covid (<2020)','1m - RMSE Covid (2020)','1m - RMSE post-Covid (>2020)','1m - RMSE no Covid', ...
               '1m - FDA all','1m - FDA pre-Covid (<2020)','1m - FDA Covid (2020)','1m - FDA post-Covid (>2020)','1m - FDA no Covid', ...
               '2m - RMSE all','2m - RMSE pre-Covid (<2020)','2m - RMSE Covid (2020)','2m - RMSE post-Covid (>2020)','2m - RMSE no Covid', ...
               '2m - FDA all','2m - FDA pre-Covid (<2020)','2m - FDA Covid (2020)','2m - FDA post-Covid (>2020)','2m - FDA no Covid', ...
               '3m - RMSE all','3m - RMSE pre-Covid (<2020)','3m - RMSE Covid (2020)','3m - RMSE post-Covid (>2020)','3m - RMSE no Covid', ...
               '3m - FDA all','3m - FDA pre-Covid (<2020)','3m - FDA Covid (2020)','3m - FDA post-Covid (>2020)','3m - FDA no Covid', ...
               'variables','groups'};
results_names = horzcat({'batch'},parameters_names(1:end-2),other_names);

end