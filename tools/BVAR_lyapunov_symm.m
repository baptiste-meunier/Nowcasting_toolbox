function [x,u] = BVAR_lyapunov_symm(a,b,qz_criterium,lyapunov_complex_threshold,method)
% Solves the Lyapunov equation x-a*x*a' = b, for b and x symmetric matrices.
% If a has some unit roots, the function computes only the solution of the stable subsystem.
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%  
% INPUTS:
%   a                           [double]    n*n matrix.
%   b                           [double]    n*n matrix.
%   qz_criterium                [double]    scalar, unit root threshold for eigenvalues in a.
%   lyapunov_complex_threshold  [double]    scalar, complex block threshold for the upper triangular matrix T.
%   method                      [integer]   Scalar, if method=0 [default] then U, T, n and k are not persistent.  
%                                                      method=1 then U, T, n and k are declared as persistent 
%                                                               variables and the schur decomposition is triggered.    
%                                                      method=2 then U, T, n and k are declared as persistent 
%                                                               variables and the schur decomposition is not performed.
% OUTPUTS
%   x:      [double]    m*m solution matrix of the lyapunov equation, where m is the dimension of the stable subsystem.
%   u:      [double]    Schur vectors associated with unit roots  
%
% ALGORITHM
%   Uses reordered Schur decomposition
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if nargin<5
    method = 0;
end
if method
    persistent U T k n
else
    if exist('U','var')
        clear('U','T','k','n')
    end
end

u = [];

if size(a,1) == 1
    x=b/(1-a*a);
    return
end

if method<2
    [U,T] = schur(a);
    e1 = abs(ordeig(T)) > 2-qz_criterium;
    k = sum(e1);       % Number of unit roots. 
    n = length(e1)-k;  % Number of stationary variables.
    if k > 0
        % Selects stable roots
        [U,T] = ordschur(U,T,e1);
        T = T(k+1:end,k+1:end);
    end
end

B = U(:,k+1:end)'*b*U(:,k+1:end);

x = zeros(n,n);
i = n;


while i >= 2

    if abs(T(i,i-1))<lyapunov_complex_threshold
        if i == n
            c = zeros(n,1);
        else
            c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
                T(i,i)*T(1:i,i+1:end)*x(i+1:end,i);
        end
        q = eye(i)-T(1:i,1:i)*T(i,i);
        x(1:i,i) = q\(B(1:i,i)+c);
        x(i,1:i-1) = x(1:i-1,i)';

        if (i == 19)
            temp = (q\(B(1:i,i)+c));
        end 
        i = i - 1;
    else
        if i == n
            c = zeros(n,1);
            c1 = zeros(n,1);
        else
            c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
                T(i,i)*T(1:i,i+1:end)*x(i+1:end,i) + ...
                T(i,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1);
            c1 = T(1:i,:)*(x(:,i+1:end)*T(i-1,i+1:end)') + ...
                 T(i-1,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1) + ...
                 T(i-1,i)*T(1:i,i+1:end)*x(i+1:end,i);

        end
        q = [  eye(i)-T(1:i,1:i)*T(i,i) ,  -T(1:i,1:i)*T(i,i-1) ; ...
               -T(1:i,1:i)*T(i-1,i)     ,   eye(i)-T(1:i,1:i)*T(i-1,i-1) ];
        z =  q\[ B(1:i,i)+c ; B(1:i,i-1) + c1 ];
        x(1:i,i) = z(1:i);
        x(1:i,i-1) = z(i+1:end);
        x(i,1:i-1) = x(1:i-1,i)';
        x(i-1,1:i-2) = x(1:i-2,i-1)';

        i = i - 2;
    end

end


if i == 1
    c = T(1,:)*(x(:,2:end)*T(1,2:end)') + T(1,1)*T(1,2:end)*x(2:end,1);
    x(1,1) = (B(1,1)+c)/(1-T(1,1)*T(1,1));
end

x = U(:,k+1:end)*x*U(:,k+1:end)';
u = U(:,1:k);
