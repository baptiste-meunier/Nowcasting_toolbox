function Date = BEQ_Add2Date(Date,n)
% This script calculates the new date which is obtained from shifting Date
% forward by n months. n must be an integer number and might be negative.
%
% Code from Bańbura, M., Belousova, I., Bodnár, K., & Tóth, M. B.
% (2023). "Nowcasting employment in the euro area", Working Paper Series,
% No 2815, European Central Bank
%
% INPUTS
% - Date [vector] = date (year/month)
% - n [scalar] = shift
%
% OUTPUT
% - Date [vector] = shifted date (year/month)
%

  years  = floor(n/12);
  mons   =   mod(n,12);
        
  if Date(2) + mons > 12
     years  = years  +  1;
     mons   = mons  - 12;
  end

  Date = Date + [years, mons];
  
end
