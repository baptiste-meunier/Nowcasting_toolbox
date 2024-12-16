function xlcol_addr = common_num2xlcol(col_num)
% This script transforms xlcol_addr into the corresponding column in Excel
% 
% INPUT
% - col_num [scalar] = column number
% 
% OUTPUT
% - xlcol_addr [char] = column ID in Excel (upper case)
% 

    n=1;

    while col_num>26*(26^n-1)/25
        n=n+1;
    end

    base_26=zeros(1,n);
    tmp_var=-1+col_num-26*(26^(n-1)-1)/25;
    
    for k=1:n
        divisor=26^(n-k);
        remainder=mod(tmp_var,divisor);
        base_26(k)=65+(tmp_var-remainder)/divisor;
        tmp_var=remainder;
    end
    xlcol_addr=char(base_26); % Character vector of xlcol address
end