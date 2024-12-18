function z = BVAR_lag(x,n,v)
% This function creates a matrix or vector of lagged values
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% USAGE: z = lag(x,n,v)
% where: x = input matrix or vector, (nobs x k)
%        n = order of lag
%        v = (optional) initial values (default=0)
% e.g.
%     z = lag(x) creates a matrix (or vector) of x, lagged 1 observations
%     z = lag(x,n) creates a matrix (or vector) of x, lagged n observations
%     z = lag(x,n,v) creates a matrix (or vector) of x, lagged n observations,
%         with initial values taking a value v.
% ------------------------------------------------------
% RETURNS: z = matrix (or vector) of lags (nobs x k)
% ------------------------------------------------------
% NOTES: if n <= 0, z = [] is returned. While you may find this
%        preverse, it is sometimes useful.
%-------------------------------------------------------
% SEE ALSO: mlag() 
%-------------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

switch(nargin)

case 1
   n = 1; v = 0;
   zt = ones(n,BVAR_cols(x))*v;
   z = [ zt; BVAR_trimr(x,0,n)];

case 2
   v = 0;
   if n < 1
   z = [];
   return;
   end;
   zt = ones(n,BVAR_cols(x))*v;
   z = [ zt; BVAR_trimr(x,0,n)];

case 3
   if n < 1
   z = [];
   return;
   end;
   zt = ones(n,BVAR_cols(x))*v;
   z = [ zt; BVAR_trimr(x,0,n)];

otherwise
error('lag: wrong # of input arguments');
end;

  
