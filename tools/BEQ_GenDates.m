function Dates = BEQ_GenDates(date1,date2,F)
% This script generates a Dates vector from date1 to date2 with frequency F
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - trim0 [vector] = start date (year/month or quarter)
% - trim1 [vector] = end date (year/month or quarter)
% - F [string] = frequency ('Q' or 'M')
%
% OUTPUT
% - Dates [vector] = date vector (year/month or quarter)
%

% Checks
  if ~((F == 'Q') || (F == 'M')) 
     error('Frequency must be either "Q" or "M"');
  end
 
% Years  
  Dates = (date1(1,1):date2(1,1))';

% Quarterly  
  if F == 'Q'
      Dates = [kron(Dates,ones(4,1)) , kron(ones(size(Dates,1),1),[3;6;9;12])];
      
      a     = date1(1,2) / 3;
      b     = date2(1,2) / 3;
      Dates = Dates(a:(end-(4-b)),:);
  end
  
% Monthly  
  if F == 'M'
      Dates = [kron(Dates,ones(12,1)) , kron(ones(size(Dates,1),1),(1:12)')];
     
      a = date1(1,2);
      b = date2(1,2);
      Dates = Dates(a:(end-(12-b)),:);
  end

end