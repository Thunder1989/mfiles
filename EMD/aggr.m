% simply do EMD + aggr on a given trace
% ts: the trace given to process
% T: length of the data in unit of day

function imf_aggr = aggr(ts,T)
th1 = 48*T/6;%6h
th2 = 48*T/(1/3);%0.5h

len = length(ts);
imf_aggr = zeros(1, len);

%load data => apply emd => aggregate imfs
% window = 10;
% data = EWMA(data', window)';%too much noise in the data, so do a exp. moving average
% fprintf('===========computing the imfs of %s/%s===========\n', folder(fn).name, file(n).name)
imfs = emd(ts);%each row is an IMF
freq = ZCR(imfs');
imf_aggr = sum(imfs(freq>th1 & freq<th2,:),1);

end