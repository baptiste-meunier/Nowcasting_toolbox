function [b0, Sig] = BEQ_OLS(Y,Xlags,cflag,sflag)
% This script estimates the model Y = Xlags*beta + u from Ordinary Least Squares
% Y might be nobs x k matrix: in this case k regressions are run
% Any missing data are skipped.
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - Y [vector/matrix] = series to be forecasted
% - Xlags [vector/matrix] = regressors
% - cflag [switch 0/1] = if 1, add constant to Z 
% - sflag [switch 0/1] = if 1, standardise data before regression
%
% OUTPUTS
% - b [vector/matrix] = coefficients
% - Sig [matrix] = covariance matrix of residuals (Choleski decomposition)
%

% Add constant & standardise 
  if  cflag
      Xlags = [ones(size(Xlags(:,1))) Xlags];
  end

% Clear all obs with missing data
  ix    = ~isnan(sum(Y,2)) & ~isnan(sum(Xlags,2));
  Yc    =  Y(ix,:);
  Zc    =  Xlags(ix,:);
  
  if  sflag
      stdY = std(Yc);
      if stdY== 0
          stdY = 1;
      end
      Yc = Yc./repmat(stdY,size(Yc,1),1);
      stdZ = std(Zc);
      stdZ(stdZ==0) = 1;
      Zc = Zc./repmat(stdZ,size(Zc,1),1);
  else
      stdY = ones(1,size(Yc,2));
      stdZ = ones(1,size(Zc,2));
  end
  
  idx0  = ~all(Zc==0);
  b0 = zeros(size(Zc,2),1);
%   sig_b0 = zeros(size(Zc,2),1);
  
  Zc = Zc(:,idx0);
 
% Estimate beta & resids
  iZc   = (Zc'*Zc)\eye(size(Zc,2));
  b     = iZc * Zc' * Yc;
  b     = b.*repmat(stdY,size(b,1),1);
  b     = b./repmat(stdZ',1,size(b,2));
  U     = Y - Xlags(:,idx0)*b;
  Sig   = cov(U(ix,:));
  b0(idx0) = b;

% Std dev and t-stat of b
%   sig_b = sqrt(diag(iZc)*diag(Sig)');
%   t_val = b ./ sig_b;
%   
%   sig_b0(idx0) = sig_b; 
%   t_val0(idx0) = t_val;

  
% Choleski decomp of Sig  
  Sig   = chol(Sig)';              % ! Adjust this later on !

end
 