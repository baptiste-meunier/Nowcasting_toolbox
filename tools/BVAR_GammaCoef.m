function r = BVAR_GammaCoef(mode,sd,plotit)
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

r.k=(2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2;
r.theta=sqrt(sd^2/r.k);

if plotit==1
    xxx=[0:.000001:mode+5*sd];
    plot(xxx,[xxx.^(r.k-1).*exp(-xxx./r.theta)*r.theta^-r.k/gamma(r.k)],'k--','LineWidth',2)
end

end