function m = BEQ_monOfQ(m)
% This script delivers the month within the quarter 
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUT
% - m [scalar] = month of the year
%
% OUTPUT
% - m [scalar] = month in quarter (=1, 2, or 3)
%

  m = mod(m,3);
  m = m + 3*(m==0);

end