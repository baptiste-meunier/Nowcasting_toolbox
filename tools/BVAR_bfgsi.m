function H = BVAR_bfgsi(H0,dg,dx)
% This code runs the optimization
% 
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% H = bfgsi(H0,dg,dx)
% dg is previous change in gradient; dx is previous change in x;
% 6/8/93 version that updates inverse hessian instead of hessian
% itself.
% Copyright by Christopher Sims 1996.  This material may be freely
% reproduced and modified.
%

if size(dg,2)>1
   dg=dg';
end
if size(dx,2)>1
   dx=dx';
end

dgp = dg';
dxp = dx';

Hdg = H0*dg;
dgdx = dgp*dx;
if (abs(dgdx) >1e-12)
   H = H0 + (1+(dgp*Hdg)/dgdx)*(dx*dxp)/dgdx - (dx*Hdg'+Hdg*dxp)/dgdx;
else
   disp('bfgs update failed.')
   disp(['|dg| = ' num2str(sqrt(dgp*dg)) '|dx| = ' num2str(sqrt(dxp*dx))]);
   disp(['dg''*dx = ' num2str(dgdx)])
   disp(['|H*dg| = ' num2str(Hdg'*Hdg)])
   H=H0;
end
%save H.dat H

end