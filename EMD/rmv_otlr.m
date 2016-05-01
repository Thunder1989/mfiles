%remove outlier
function output = rmv_otlr(ts)

tmp=ts(1);
for k=2:length(ts)
    if abs(abs(ts(k))-abs(tmp))>3*min(abs(tmp),abs(ts(k))) && tmp~=0
        ts(k) = tmp;
        continue;
    end
    tmp = ts(k);
end
output = EWMA(ts', 50);%too much noise in the data, so do a exp. moving average

end