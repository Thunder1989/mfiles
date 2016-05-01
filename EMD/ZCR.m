function [count, rate] = ZCR(x)
%   zero crossing rate, input is a vector or a matrix
%   if x is a vector returns the zero crossing rate of the vector
%   if x is a matrix returns a row vector with the zero crossing rate of
%   the columns values

count = sum(abs(diff(x>0))); %count of ZC
rate = count/length(x);

end
