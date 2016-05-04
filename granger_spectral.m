%checking holes in data
%remove hole period from ts
%511 has only less than 3 days of data
clear
clc
path = '/Users/hdz/Downloads/code/Data/KETI_oneweek/';
folder = dir(path);
first = 1377241200; %8/23/2013 0:0:0am PST
last = 1378018799; %8/31/2013 11:59:59pm PST
ts = first:last;
roomNum = size(folder,1);
data = {};
sensorName = {};
gt_type = [];
gt_room = [];
ctr = 0; % # of streams
for n = 4:roomNum
    path1 = strcat(path, folder(n).name,'/');   %path to each room
    file = dir(strcat(path1,'*.csv'));
    sensorNum = size(file,1);
    for i=1:sensorNum
        if i~=4 %skip pir
            ctr = ctr+1;
            gt_type = [gt_type i];
            gt_room = [gt_room n-3];

            filename = [path1, file(i).name];
            sensorName{ctr} = strcat(folder(n).name, '/', file(i).name);
            input = csvread(filename);
            input = input(:,1);
            ts(ts<input(1)) = [];
            ts(ts>input(end)) = [];

            delta = diff(input); %input is a column vector
            idx = (delta>=60*5);
            t_end = input([false idx']); %find the ending ts for each hole period
            dt = delta(idx);
            t_start = t_end - dt;
            for j=1:length(t_start)
                ts(ts>t_start(j) & ts<t_end(j)) = [];
            end
        end
    end
end

%% get valid timestamp intervals without holes for streams
delta = diff(ts);
idx = (delta>1);
t_end = ts([false idx]); %end of the hole
t_start = t_end - delta(idx);
t = [t_start; t_end];
t = reshape(t,1,numel(t));
boundary = [ts(1) t ts(end)];
fprintf('valid timestamp intervals done!\n');

%% get the shortest data length
% align the readings and resample
data_length = intmax;
for n = 4:roomNum
%     figure
    path1 = strcat(path, folder(n).name,'/');%path to the csv files
    file = dir(strcat(path1,'*.csv'));
    len = size(file,1);
    for i = 1:len
        filename = [path1, file(i).name];
        input = csvread(filename);
        ts = input(:,1);
        input = input(:,2);
        
        input = input(ts>=boundary(1) & ts<=boundary(end));
        ts = ts(ts>=boundary(1) & ts<=boundary(end));
        
%         for j=2:2:length(boundary)-1
%             hole_start = boundary(j);
%             hole_end = boundary(j+1);
%             input(ts>hole_start & ts<hole_end) = [];
%             ts(ts>hole_start & ts<hole_end) = [];
%         end
        
        if length(input) < data_length
            data_length = length(input);
        end
    end
end
fprintf('minimum length computed!\n');

%% align the readings and resample
%data_length = 1263; %with 511
data_length = 6485; %w/o 511
input_data = zeros(ctr, data_length);
ctr = 0;
for n = 4:roomNum
    path1 = strcat(path, folder(n).name,'/');   %path to each room
    file = dir(strcat(path1,'*.csv'));
    sensorNum = size(file,1);
    for i = 1:sensorNum
        filename = [path1, file(i).name];
        data = csvread(filename);
        ts = data(:,1);
        data = data(:,2);
        
        data = data(ts>=boundary(1) & ts<=boundary(end));
        ts = ts(ts>=boundary(1) & ts<=boundary(end));

%         for j=2:2:length(boundary)-1
%             hole_start = boundary(j);
%             hole_end = boundary(j+1);
%             data(ts>hole_start & ts<hole_end) = [];
%         end

        len = length(data);
        if len > data_length+300
            [p,q] = rat(data_length/len);
            data = resample(data, p, q, 5, 20);
            data = data(1:data_length);
        else
            data = data(1:data_length);
        end
        
        ctr = ctr + 1;
        input_data(ctr,:) = data;
        
%         %remove outlier
%         if i~=3
%             tmp = data(1);
%             for k=2:length(data)
%                 if abs(abs(data(k))-abs(tmp))>3*min(abs(tmp),abs(data(k))) && tmp~=0
%                     data(k) = tmp;
%                     continue;
%                 end
%                 tmp = data(k);
%             end
%         end

    end
end
fprintf('data processing done!\n');

%% normalization
data = input_data;
data(4:5:end,:) = []; %remove pir since mostly are zeros
sensorNum = size(data,1); % N stream

for i = 1:sensorNum
    tmp = data(i,:);
%     t_min = min(tmp);
%     t_max = max(tmp);
%     data(i,:) = (tmp-t_min) / (t_max-t_min);
    u = mean(tmp);
    d = std(tmp);
    data(i,:) = (tmp-u) / d;
end
fprintf('normalization done!\n');

%% spectral clustering method
clc
w = zeros(sensorNum, sensorNum);
fprintf('lasso computing started...\n');
res = [];
for b = 0.02:0.02:0.02
%     b = 0.14;
    for i = 1:sensorNum
        cur = data(i,:); %1 by d
        src = data;
        src(i,:) = []; %N-1 by d
        % 0.015~0.03 for max-min normalization, 0.06 gives no zero rows
        % 0.14~0.2 for u-std normalization, 0.14 gives no zero rows - 0.02-0.14 all resonable
        coef = lasso(src', cur', 'Lambda', b); %N-1 by 1
        idx = 1;
        for j=1:length(coef)
            if(idx==i) 
                idx = idx+1;
            end
            w(i,idx) = coef(j);
            idx = idx+1;
        end
    %     fprintf('lasso computing itr %d done!\n', i);
    end
% fprintf('lasso computing done!\n');

    [g_, idx] = sort(gt_type); %ground truth indexing
    w_ = w(idx,idx); %re-ordered by true type ID
    w_ = max(w_, w_'); %symmetrize N by N
    d = diag(sum(w_,2));
%     find(diag(d)==0)
%     fprintf('all zero rows found!\n');
%     pause
    l = d - w_; %unormalized Laplacian
    [evc, evl] = eig(l); %N by N, each column in evc is an eigenvector
    idx = find(diag(evl)>0);
    input = evc(:,idx(1:4));
    c_idx = kmeans(input,4,'Distance','cosine');
%     adjrand(c_idx, g_)
%     figure
%     spy(w_)
    res = [res adjrand(c_idx, g_)];
end
fprintf('lasso-based spectral clustering done...\n');
% figure
% plot(0.015:0.005:0.015, res)

%% baseline w
data = input_data;
data(4:5:end,:) = []; %remove pir since mostly are zeros
sensorNum = size(data,1); % N stream, each row is a stream
% normalization
for i = 1:sensorNum
    tmp = data(i,:);
%     t_min = min(tmp);
%     t_max = max(tmp);
%     data(i,:) = (tmp-t_min) / (t_max-t_min);
    u = mean(tmp);
    d = std(tmp);
    data(i,:) = (tmp-u) / d;
end
options = [];
options.NeighborMode = 'KNN';
options.WeightMode = 'Cosine';

% options.t = 1;
[g_, idx] = sort(gt_type);
% [g_, idx] = sort(gt_room);
res_ = [];
% for i = 5:2:5
    options.k = 13;
    w_c = constructW(data, options);
    w_ = w_c(idx,idx);
    figure
    spy(w_)
    w_ = max(w_, w_'); %symmetric N by N
    d = diag(sum(w_,2));
%     find(diag(d)==0)
%     fprintf('all zero rows found!\n');
%     pause
    l = d - w_
    [evc, evl] = eigs(l,4,'sa'); %N by N, each column is a principal component
%     idx = find(diag(evl)>0);
%     input = evc(:,idx(1:4));
    c_idx = kmeans(evc,4,'Distance','cosine');
    adjrand(c_idx, g_)
%     res_ = [res_ adjrand(c_idx, g_)];
% end
% figure
% plot(5:2:21, res_)

%% plot each for inspection
sensorNum = size(data,1); % N stream
figure
hold on
for i = 1:sensorNum
    cur = data(i,:); %1 by d
    plot(cur)
    pause
end