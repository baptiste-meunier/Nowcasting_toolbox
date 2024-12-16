function mInv = BVAR_inv2(m)
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%

    mInv = m\eye(size(m));
    
%     mInv(abs(mInv)<1e-16) = 0;

end