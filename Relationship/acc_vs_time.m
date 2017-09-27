% clear
% clc
bldg = [642];
% bldg = [596];
for i = 1:length(bldg)
    bid = bldg(i);
    for week = 1:8
        correlation_onepair(bid,week);
    end
end