function [ M_matrix ] = estimateMforDifferentRegions( values, value_low, value_high, del2P )
% Average local M over given regions in the M grid,
% and stores the esimates in a matrix
%
% Arg:
%   values (matrix): values that you want to limit with respect to 
%                    (r matrix or pO2 matrix)
%   value_low (int or double): 
%   value_high (int or double): 
%   del2P (matrix): M grid, which you average M from
%
% Returns:
%   M_matrix (matrix): estimated values of M, averaged over a region parted
%                      by value_low and value_high 

s = value_low(2)-value_low(1);
M_matrix = NaN(length(value_high), length(value_low));
i = 1;
str = size(M_matrix);
for low = value_low
    if low < value_high(1)
        highStart = value_high(1);
        j = 1;
    elseif low == value_high(1)
        highStart = low+s;
        j = 2;
        g = 2;
    else
        highStart = low+s;
        j = g+1;
        g = g+1;
    end
    if j > str(1)
        continue
    end
    for high = highStart:s:value_high(end)
        M_list = del2P(values > low & values < high);
        M_est = mean(M_list);
        M_matrix(j, i) = M_est;       
        j = j+1;
    end
    i = i+1;
end
end

