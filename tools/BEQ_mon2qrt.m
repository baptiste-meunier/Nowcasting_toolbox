 function datq = BEQ_mon2qrt(datm)
% This script creates the quarterly date in which month datm lies
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUT
% - datm [vector] = monthly date (year/month)
%
% OUTPUT
% - datq [vector] = quarterly date (year/quarter)
%

  if sum(size(datm) ~= [1 2]) > 0
      datq = NaN*zeros(1,2);
      return
  end 
  
  datq(1,1) = datm(1,1);
  datq(1,2) = ceil(datm(1,2)/3)*3;

 end