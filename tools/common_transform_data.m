function [x_out] = common_transform_data(x_in,transf_in,revert)
% This script transform the data
% It can also revert the data back to levels (for bridge equations)
% 
% INPUTS 
% - x_in [matrix] = original series
% - transf_in [vector] = indicator for transformations
%                    0 = levels
%                    1 = one-period dlog change
%                    2 = first difference
%                    3 = one-period dlog change times 4 (= annualized qoq for quarterly data)
% - revert [scalar -1/1] = switch on whether to do transform data in stationary data (1) or back to levels (-1)
%
% OUTPUTS
% - x_out [matrix] = transformed dataset
%

if revert == 1

    x_out = real(log(x_in(2:end,:))-log(x_in(1:end-1,:)))*100;                                          % One-period growth rates
    x_out(:,transf_in==0) = x_in(2:end,transf_in==0);                                                   % No transformation
    x_out(:,transf_in==2) = x_in(2:end,transf_in==2)-x_in(1:end-1,transf_in==2);                        % First difference
    x_out(:,transf_in==3) = real(log(x_in(2:end,transf_in==3))-log(x_in(1:end-1,transf_in==3)))*400;    % One-period growth rates times 4 (= annualized for quarterly data)

elseif revert == -1

    x_out = ones(size(x_in,1)+1,size(x_in,2));
    
    for cc = 1:size(x_out,2)
        for rr = 2:size(x_out,1)
            if ~isnan(x_in(rr-1,cc))
                switch transf_in(cc)
                    case 1
                        x_out(rr,cc) = exp(log(x_out(rr-1,cc))+x_in(rr-1,cc)/100); % one-period growth rates
                    case 2
                        x_out(rr,cc) = x_out(rr-1,cc)+x_in(rr-1,cc); % first difference
                    case 3
                        x_out(rr,cc) = exp(log(x_out(rr-1,cc))+x_in(rr-1,cc)/400); % one-period growth rates times 4 (annualized)
                end
            else
                if rr < size(x_out,1)
                    if isnan(x_in(rr,cc))
                        x_out(rr,cc) = NaN;
                    end
                else
                    x_out(rr,cc) = NaN;
                end
            end
        end
    end 

    x_out = x_out(2:end,:); % deleting first row (which was for initializing)
    x_out(:,transf_in==0) = x_in(:,transf_in==0); % putting levels the same

end % end of condition on revert

end