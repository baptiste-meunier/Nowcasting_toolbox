function [converged, decrease] = DFM_em_converged(loglik, previous_loglik, threshold, check_increased)
% This function checks the convergence of the EM algorithm
%
% Code from Bańbura, M., and Modugno, M. (2014). “Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%
% NB: Convergence means the slope of the log-likelihood function falls below 'threshold',
%     i.e., |f(t) - f(t-1)| / avg < threshold,
%     where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
%     'threshold' defaults to 1e-4.
%     This stopping criterion is from Numerical Recipes in C p423
%

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        %fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end

end