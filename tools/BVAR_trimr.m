function z = BVAR_trimr(x,n1,n2)
% This function return a matrix (or vector) x stripped of the specified rows.
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% -----------------------------------------------------
% USAGE: z = trimr(x,n1,n2)
% where: x = input matrix (or vector) (n x k)
%       n1 = first n1 rows to strip
%       n2 = last  n2 rows to strip
% NOTE: modeled after Gauss trimr function
% -----------------------------------------------------
% RETURNS: z = x(n1+1:n-n2,:)
% -----------------------------------------------------
%
% written by:
% James P. LeSage, Dept of Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

  [n junk] = size(x);
  if (n1+n2) >= n; 
     error('Attempting to trim too much in trimr');
  end;
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
  
