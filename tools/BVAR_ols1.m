function r = BVAR_ols1(y,x) 
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

if (nargin ~= 2); error('Wrong # of arguments to ols');
else
 [nobs nvar] = size(x); [nobsy junk] = size(y);
 if (nobs ~= nobsy); error('x and y must have same # obs in ols');
 end;
end;

r.nobs=nobs;
r.nvar=nvar;
r.bhatols=(x'*x)\(x'*y);
r.yhatols=x*r.bhatols;
r.resols=y-r.yhatols;
r.sig2hatols=(r.resols'*r.resols)/(nobs-nvar);
r.sigbhatols=r.sig2hatols*(inv(x'*x));
r.XX=(x'*x);
r.R2=var(r.yhatols)/var(y);


