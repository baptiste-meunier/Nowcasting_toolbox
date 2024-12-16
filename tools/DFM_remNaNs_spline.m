function [X,indNaN] = DFM_remNaNs_spline(X,options)
% This function fills missing values (NaN) in the input dataset of the DFM
%
% Code from Bańbura, M., and Modugno, M. (2014). “Maximum likelihood 
% estimation of factor models on datasets with arbitrary pattern of 
% missing data”, Journal of Applied Econometrics, 29(11), 133–160
%

% dbstop('3')
[T,N]=size(X);
k=options.k;
indNaN=isnan(X);
switch options.method
    case 1  % replace all the missing values
        for i = 1:N  
            x = X(:,i);
            x(indNaN(:,i))=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(indNaN(:,i))=x_MA(indNaN(:,i));
            X(:,i)=x;
        end
    case 2 %replace missing values after removing leading and closing zeros
        rem1=(sum(indNaN,2)>N*0.8);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            t1 = min(find(~isnanx));
            t2 = max(find(~isnanx));
            x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
            isnanx = isnan(x);
            x(isnanx)=median(x(isnanx==0,:));
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(isnanx)=x_MA(isnanx);
            X(:,i)=x;
        end
    case 3 %only remove rows with leading and closing zeros
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN = isnan(X);
    case 4%remove rows with leading and closing zeros & replace missing values
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            t1 = min(find(~isnanx));
            t2 = max(find(~isnanx));
            x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
            isnanx = isnan(x);
            x(isnanx)=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(isnanx)=x_MA(isnanx);
            X(:,i)=x;
        end
    case 5 %replace missing values  
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            t1 = min(find(~isnanx));
            t2 = max(find(~isnanx));
            x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
            isnanx = isnan(x);
            x(isnanx)=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(isnanx)=x_MA(isnanx);
            X(:,i)=x;
        end
end
