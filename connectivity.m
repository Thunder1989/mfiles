%% calculate the connectivity score
seq = 1:length(sub);
list = [];
% sub(abs(sub)>=0.08) = 1;
% sub(abs(sub)<0.08) = 0;

while numel(seq)>3
    first = seq(1);
    seq = seq(2:end);
    comb = nchoosek(seq,3);
    w = 0;
    index = [];
    for i=1:size(comb,1)
        idx = first;
        idx = [idx comb(i,:)];
        tmp = sub(idx, idx);
        tmp_weight = sum(sum(abs(tmp)));
        if tmp_weight > w
            w = tmp_weight;
            index = idx;
        end
    end
    list = [list; index];
    for i=1:length(index)
        seq(seq==index(i))=[];
    end
end