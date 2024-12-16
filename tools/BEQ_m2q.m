 function [X_q,date_q] = BEQ_m2q(X,date,aggtype)
% This script convertss monthly dates and data into quarterly format
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
%  INPUTS
%  - X [matrix] = monthly data
%  - date [matrix] = monthly dates (year/month)
%                    quarters 1, 2, 3, 4 are indicated by  3, 6, 9, 12 
%  - aggtypes [string] = 'S'  Sum of monthly data
%                      = 'A'  Average of monthly data
%                      = 'M'  Aggregation of month-on-month rates
%                      = '3'  3rd month of each quarter
%
%  OUTPUTS
%  - date_q [matrix] = quarterly dates (number of observations divided by 3)
%  - X_q [matrix] = quarterly data (number of observations divided by 3)
%

if size(X,1) ~= size(date,1)
 error('Data and date vectors must be of equal length')
end   

% Expand date to start in 1st and end in 3rd months of the quarter 
d1 = BEQ_Add2Date(date(1,:), ...
                  1-BEQ_monOfQ(date(1,2))); % meaning we want to extend (backward) to first month of quarter 
d2 = BEQ_Add2Date(date(end,:), ...
                  3-BEQ_monOfQ(date(end,2))); % meaning we want to extend (forward) to last month of quarter

% Trim to expanded dates 
[X,date] = BEQ_TrimData(X,date,d1,d2,'M'); % trims or extend data to dates d1 and d2
  
% Generate date vector and data
[~, k] = size(X);
 
for j = 1:k
  if       aggtype{j} == 'S' ;  X(:,j) = filter([1 1 1]       ,1,X(:,j));
  elseif   aggtype{j} == 'A' ;  X(:,j) = filter([1 1 1]    ./3,1,X(:,j));
  elseif   aggtype{j} == 'M' ;  X(:,j) = filter([1 2 3 2 1]./3,1,X(:,j));
  end
end

X_q = X(3:3:end,:);
X_q(1,ismember(aggtype,'M')) = NaN;
date_q    = date(3:3:end,:);

end 