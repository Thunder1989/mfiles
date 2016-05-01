% ma bi de, zhe shi wo xie guo zui e xin de code
clear
clc
path = 'C:\Users\dh5gm\Dropbox\SDB\KETI_tmp\';
folder = dir(path);
first = 1371024000;
last = 1371628800;
period = last - first;

num = size(folder,1);
files = cell((num-2)*4,1);
scores = cell((num-2)*4,(num-2)*4);
for i=1:length(scores)
    scores{i,i} = 0;
end

%get all the file names
count=1;
for n=3:num
    path1 = strcat(path, folder(n).name,'\');%path to the csv files
    file = dir(strcat(path1,'*.csv'));
    len=size(file,1);
    for i=1:len
        filename = [path1, file(i).name];
        files{count} = filename;
        count = count+1;
    end
end

for i=1:length(files)
    i
    ts1 = csvread(files{i});
    for j=1:i-1
        j
        ts = first:last;
        t1 = ts1(:,1);
        d1 = ts1(:,2);
        ts2 = csvread(files{j});
        t2 = ts2(:,1);
        d2 = ts2(:,2);
        if ceil(j/4)==ceil(i/4) %within the same room, directly compute the score
            d1 = rmv_otlr(d1);
            d1 = aggr(d1,7);
            d2 = rmv_otlr(d2);
            d2 = aggr(d2,7);
            
            tmp = 1;
            len = min(length(t1),length(t2));
            corr = zeros(period/(15*60),1);
            count = 1;
            for k = 2:len
                if t1(k)>t1(tmp)+15*60 && k-tmp>40 %15-min window
                    w1 = d1(tmp:k);
                    w2 = d2(tmp:k);
                    res = corrcoef(w1,w2);
                    corr(count) = res(1,2);
%                     fprintf('=======computing a score=======\n');
                    count = count+1;
                    tmp = k;
                end
            end
            corr = corr(corr~=0);
            scores{i,j} = corr;
            
        else %two traces from diff rooms
            % find common time period ---start
            start = max(t1(1),t2(1));
            if start>first
                ts = ts(ts>start);
            end
            stop = min(t1(end),t2(end));
            if stop<last
                ts = ts(ts<stop);
            end

            delta = diff(t1);
            idx = (delta>=60);
            p_end = t1([false idx']);
            dt = delta(idx);
            p_start = p_end - dt;
            for k=1:length(p_start)
                ts = ts(ts<=p_start(k) | ts>=p_end(k));
            end
            
            delta = diff(t2);
            idx = (delta>=60);
            p_end = t2([false idx']);
            dt = delta(idx);
            p_start = p_end - dt;
            for k=1:length(p_start)
                ts = ts(ts<=p_start(k) | ts>=p_end(k));
            end
            
            delta = diff(ts);
            idx = (delta>1);
            p_end = ts([false idx]);
            p_start = p_end - delta(idx);
            t = [p_start; p_end];
            t = reshape(t,1,numel(t));
            limit = [ts(1) t ts(end)];
            % find common time period ---end
                         
            %get the min length of two traces
            for k=2:2:length(limit)-1
                bd1 = limit(k);
                bd2 = limit(k+1);
                t1 = t1(t1<=bd1 | t1>=bd2);
            end
            len1 = length(t1);
            for k=2:2:length(limit)-1
                bd1 = limit(k);
                bd2 = limit(k+1);
                t2 = t2(t2<=bd1 | t2>=bd2);
            end
            len2 = length(t2);           
            min_length = min(len1,len2);
            
            %resample traces
            %on trace1
            for k=2:2:length(limit)-1
                bd1 = limit(k);
                bd2 = limit(k+1);
                d1 = d1(t1<=bd1 | t1>=bd2);
            end
            len = length(d1);
            if len > min_length+500
                [p,q] = rat(min_length/len);
                d1 = resample(d1, p, q);
                d1 = d1(1:min_length);
            else
                d1 = d1(1:min_length);
            end
            %on trace2
            for k=2:2:length(limit)-1
                bd1 = limit(k);
                bd2 = limit(k+1);
                d2 = d2(t2<=bd1 | t2>=bd2);
            end
            len = length(d2);
            if len > min_length+500
                [p,q] = rat(min_length/len);
                d2 = resample(d2, p, q);
                d2 = d2(1:min_length);
            else
                d2 = d2(1:min_length);
            end

            %do the same thing as traces in one room
            ts = [];
            if length(t1) == min_length
                ts = t1;
            else
                ts = t2;
            end
            d1 = rmv_otlr(d1);
            d1 = aggr(d1,7);
            d2 = rmv_otlr(d2);
            d2 = aggr(d2,7);
  
            tmp = 1;
            corr = zeros(period/(15*60),1);
            count = 1;
            for k = 2:min_length
                if ts(k)>ts(tmp)+15*60 && k-tmp>40 %15-min window
                    w1 = d1(tmp:k);
                    w2 = d2(tmp:k);
                    res = corrcoef(w1,w2);
                    corr(count) = res(1,2);
%                     fprintf('=======computing a score=======\n');
                    count = count + 1;
                    tmp = k;
                end
            end
            corr = corr(corr~=0);
            scores{i,j} = corr;
        end
    end
end