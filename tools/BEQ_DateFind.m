function date_ix = BEQ_DateFind(Dates,date)
% This script finds the indices of a vector date in Date vector Dates
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - Dates [matrix] = vector of dates (year/month or quarter)
% - date  [vector] = date to find (can be several of them)
%
% OUTPUT
% - i [scalar] = index of date in Date
%                if date(i,:) is not found in Dates then NaN is returned
%
%

  k       = size(date,1);
  date_ix = nan * zeros(k,1);
  
  for i = 1:k
      j = find(Dates(:,1) == date(i,1) & Dates(:,2) == date(i,2));
     
     if ~isempty(j)
         date_ix(i) = j;
     end    
  end

end