function [logML,betadraw,drawSIGMA] = BVAR_logMLVAR_formcmc(par,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,draw,hyperpriors,priorcoef);
% This function computes the log-posterior (or the logML if hyperpriors=0), 
% and draws from the posterior distribution of the coefficients and of the 
% covariance matrix of the residuals of the BVAR 
% 
% BVAR from of Giannone, D., Lenza, M., and Primiceri, G. E. (2015). 
% “Prior Selection for Vector Autoregressions”, 
% The Review of Economics and Statistics, 97(2), 436–451
%
% Code from Cimadomo, J., Giannone, G., Lenza, M., Monti, F., & Sokol, A. 
% (2022). "Nowcasting with large Bayesian vector autoregressions", 
% Journal of Econometrics, 231, 500-519
%
% Last modified: 07/01/2013


%% hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=par(1);
d=n+2;
theta=MIN.theta;
miu=MIN.miu;

if mn.psi==0;
    psi=SS*(d-n-1);
    if sur==1;
        theta=par(2);
        if noc==1;
            miu=par(3);
        end
    elseif sur==0;
        theta = 1;
        if noc==1;
            miu=par(2);
        end
    end
elseif mn.psi==1;
    psi=par(2:n+1);
    if sur==1;
        theta=par(n+2);
        if noc==1;
            miu=par(n+3);
        end
    elseif sur==0;
        theta = 1;
        if noc==1;
            miu=par(n+2);
        end
    end
end

if mn.alpha==0;
    alpha=2;
elseif mn.alpha==1;
    alpha=par(end);
end


%% return a very low value of the posterior if the parameters are outside the bounds
if sum([lambda;psi;theta;miu;alpha]<[MIN.lambda;MIN.psi;MIN.theta;MIN.miu;MIN.alpha])>0 | ...
        sum([lambda;psi;theta;miu]>[MAX.lambda;MAX.psi;MAX.theta;MAX.miu])>0
    logML=-10e15;
    betadraw=[];
    drawSIGMA=[];
    return
else
    
    %% priors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=1+n*lags;
    omega=zeros(k,1);
    omega(1)=Vc;
    for i=1:lags
        omega(1+(i-1)*n+1:1+i*n)=(d-n-1)*(lambda^2)*(1/(i^alpha))./psi;
    end
    
    % prior scale matrix for the covariance of the shocks
    PSI=diag(psi);
    
    Td=0;
    xdsur=[];
    ydsur=[];
    xdnoc=[];
    ydnoc=[];
    % dummy observations if sur and/or noc = 1
    if sur==1;
                
        xdsur=[1/theta (1/theta)*repmat(y0,1,lags)];
        ydsur=(1/theta)*y0;
        
        y=[y;ydsur];
        x=[x;xdsur];
        Td=1;
    end
    
    if noc==1;
        
        ydnoc=(1/miu)*diag(y0); ydnoc(pos,pos)=0;
        xdnoc=[zeros(n,1) (1/miu)*repmat(diag(y0),1,lags)];
        
        y=[y;ydnoc];
        x=[x;xdnoc];
        Td=Td+n;
    end
    
    
    %% output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % posterior mode of the VAR coefficients
    xT = x';
    xTx = xT*x;
%     betahat=(xTx+diag(1./omega))\(xT*y+diag(1./omega)*b);
        betahat=(xTx+diag(1./omega))\(xT*y+b./omega);
    % VAR residuals
    epshat=y-x*betahat;
    epspesT = epshat'*epshat;
betaTerm = (betahat-b)'*((betahat-b)./omega);
    
    % logML
    T=T+Td;
    
%     aaa=diag(sqrt(omega))*(x'*x)*diag(sqrt(omega));
    aaa=sqrt(omega)'.*xTx.*sqrt(omega);

%     bbb=diag(1./sqrt(psi))*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*diag(1./sqrt(psi));
    bbb= sqrt(psi)'.\(epspesT + betaTerm)./sqrt(psi);

    eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
    eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
    
    logML = - n*T*log(pi)/2 + sum(gammaln((T+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
        - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
    
    if sur==1 | noc==1;
        yd=[ydsur;ydnoc];
        xd=[xdsur;xdnoc];
        
        % prior mode of the VAR coefficients
        % betahatd=(xd'*xd+diag(1./omega))\(xd'*yd+diag(1./omega)*b);
        betahatd=b;     % this is the case for our priors (the line above delivers the same but is numerically not very stable)

        
        % VAR residuals at the prior mode
        epshatd=yd-xd*betahatd;
        
%         aaa=diag(sqrt(omega))*(xd'*xd)*diag(sqrt(omega));
    aaa=sqrt(omega)'.*(xd'*xd).*sqrt(omega);
    
%         bbb=diag(1./sqrt(psi))*(epshatd'*epshatd + (betahatd-b)'*diag(1./omega)*(betahatd-b))*diag(1./sqrt(psi));
    bbb= sqrt(psi)'.\(epshatd'*epshatd + (betahatd-b)'*((betahatd-b)./omega))./sqrt(psi);

        
        eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
        eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
        
        % normalizing constant
        norm = - n*Td*log(pi)/2 + sum(gammaln((Td+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
            - Td*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
        
        logML=logML-norm;
        
        if hyperpriors==1;
            logML=logML+logGammapdf(lambda,priorcoef.lambda.k,priorcoef.lambda.theta);
            if sur==1;
                logML=logML+logGammapdf(theta,priorcoef.theta.k,priorcoef.theta.theta);
            end
            if noc==1;
                logML=logML+logGammapdf(miu,priorcoef.miu.k,priorcoef.miu.theta);
            end
            if mn.psi==1;
                logML=logML+sum(logIG2pdf(psi/(d-n-1),priorcoef.alpha.PSI,priorcoef.beta.PSI));
            end
        end
    end
    
    % takes a draw from the posterior of SIGMA and beta, if draw is on
    if draw==0
        betadraw=[];
        drawSIGMA=[];
    elseif draw==1
        S=(PSI + epspesT + betaTerm);       
        [V E]=eig(S);
        Sinv=real(V*(V'./abs(diag(E)))); %added real
        eta=mvnrnd(zeros(1,n),Sinv,T+d);
        drawSIGMA=(eta'*eta)\eye(n);
        drawSIGMA = (drawSIGMA+drawSIGMA')/2;
        %[cholSIGMA,junk]=chol(drawSIGMA);
        %betadraw=betahat+mvnrnd(zeros(k,1),(x'*x+diag(1./omega))\eye(k),n)'*cholSIGMA;
        cholSIGMA=cholred(drawSIGMA);
        cholZZinv = cholred((xTx+diag(1./omega))\eye(k));
        betadraw=betahat + cholZZinv'*randn(size(betahat))*cholSIGMA;
    end
end

function r=logGammapdf(x,k,theta);
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logIG2pdf(x,alpha,beta);
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);

function C = cholred(S)
[v,d] = eig(S);
d = diag(real(d));
v = real(v);
warning off
scale = mean(diag(S))*1e-12;
J = (d>scale);
C = zeros(size(S));
C(J,:) = v(:,J)'.*(d(J).^(1/2));