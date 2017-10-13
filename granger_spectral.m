%checking holes in data
%remove hole period from ts
%511 has only less than 3 days of data
% clear
% clc
path = 'D:\TraneData\KETI_oneweek\';
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
for n = 3:roomNum
    path1 = strcat(path, folder(n).name,'/');   %path to each room
    file = dir(strcat(path1,'*.csv'));
    sensorNum = size(file,1);
    for i=1:sensorNum
        if i~=4 %skip pir
            ctr = ctr+1;
            gt_type = [gt_type i];
            gt_room = [gt_room n-3];

%             filename = [path1, file(i).name];
%             sensorName{ctr} = strcat(folder(n).name, '/', file(i).name);
%             input = csvread(filename);
%             input = input(:,1);
%             ts(ts<input(1)) = [];
%             ts(ts>input(end)) = [];
% 
%             delta = diff(input); %input is a column vector
%             idx = (delta>=60*5);
%             t_end = input([false idx']); %find the ending ts for each hole period
%             dt = delta(idx);
%             t_start = t_end - dt;
%             for j=1:length(t_start)
%                 ts(ts>t_start(j) & ts<t_end(j)) = [];
%             end
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
load('keti_processed_ssc.mat');
data = input_data;
data(4:5:end,:) = []; %remove pir since mostly are zeros
sensorNum = size(data,1); % N stream

for i = 1:sensorNum
    tmp = data(i,:);
%     t_min = min(tmp);
%     t_max = max(tmp);
%     data(i,:) = (tmp-t_min) / (t_max-t_min);
    u = mean(tmp);
    D = std(tmp);
    data(i,:) = (tmp-u) / D;
end
fprintf('normalization done!\n');

%% spectral clustering on Dantzig
clc
W = zeros(sensorNum, sensorNum);
fprintf('lasso computing started...\n');
res = [];
residual = zeros(size(data));
for b = 0.08:0.02:0.08
%     b = 0.14;
    for i = 1:sensorNum
        cur = data(i,:); %1 by T
        src = data;
        cur_type = gt_type(i);
        mask_idx = find(gt_type==cur_type);
        for j = 1:length(mask_idx)
            id = mask_idx(j);
            src(id,:) = zeros(size(cur)); %extra constraint - set each of the same type to zero
        end
        src(i,:) = zeros(size(cur)); %set self to 0
        % 0.015~0.03 for max-min normalization, 0.06 gives no zero rows
        % 0.14~0.2 for u-std normalization, 0.14 gives no zero rows - 0.02-0.14 all resonable
        coef = lasso(src', cur', 'Lambda', b); %d by 1
        W(i,:) = coef;
        residual(i,:) = abs(cur' - src' * coef);
    %     fprintf('lasso computing itr %d done!\n', i);
    end
% fprintf('lasso computing done!\n');

    [g_, idx] = sort(gt_type); %ground truth indexing
    W_ = W(idx,idx); %re-ordered by true type ID
    W_ = max(W_, W_'); %symmetrize N by N
    D = diag(sum(W_,2));
%     find(diag(d)==0)
%     fprintf('all zero rows found!\n');
%     pause
    L = D - W_; %unormalized Laplacian
    [evc, evl] = eig(L); %N by N, each column in evc is an eigenvector
    idx = find(diag(evl)>0);
    input = evc(:,idx(1:4));
    c_idx = kmeans(input,4,'Distance','cosine');
    adjrand(c_idx, g_)
    res = [res adjrand(c_idx, g_)];
end
fprintf('lasso-based spectral clustering done...\n');

%% CC
clc
res = [];
k = 4;
for i=1:10
    corr = corrcoef(data');
    corr = corr - eye(size(corr));
    W_ = corr;
    [g_, idx] = sort(gt_type); %ground truth indexing
    W_ = W_(idx,idx); %re-ordered by true type ID
    W_ = max(W_, W_'); %symmetrize w_, N by N
    D = diag(sum(W_,2));
    L = D - W_; %unormalized Laplacian
    [evc, evl] = eig(L); %each column of evc is an eigenvector
    idx = find(diag(evl)>=0);
    input = evc(:,idx(1:k));
    % input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
    c_idx = kmeans(input,k,'Distance','cosine');
    ari = adjrand(c_idx, g_);
    res = [res ari];
end
mean(res)
std(res)

%% Cosine
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
    D = std(tmp);
    data(i,:) = (tmp-u) / D;
end
options = [];
options.NeighborMode = 'KNN';
options.WeightMode = 'Cosine';
options.bTrueKNN = 1;
% options.t = 1;
[g_, idx] = sort(gt_type);
% [g_, idx] = sort(gt_room);

res = [];
for run=1:10
    for i = 5:2:5
    options.k = i;
    W_c = constructW(data, options);
    W_ = W_c(idx,idx);
%     figure
%     spy(W_)
    W_ = max(W_, W_'); %symmetric N by N
    D = diag(sum(W_,2));
%     find(diag(d)==0)
%     fprintf('all zero rows found!\n');
%     pause
    L = D - W_;
    [evc, evl] = eigs(L,4,'sa'); %N by N, each column is a principal component
%     [evc, evl] = eig(L); %N by N, each column in evc is an eigenvector
%     idx = find(diag(evl)>0);
%     input = evc(:,idx(1:4));
    c_idx = kmeans(evc,4,'Distance','cosine');
%     adjrand(c_idx, g_)
    res = [res adjrand(c_idx, g_)];
    end
end
mean(res)
std(res)
% figure
% plot(5:2:21, res_)

%% DTW
clc
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
    D = std(tmp);
    data(i,:) = (tmp-u) / D;
end

[g_, idx] = sort(gt_type);
% [g_, idx] = sort(gt_room);
res = [];
W_ = zeros(sensorNum,sensorNum);
combo = nchoosek(1:sensorNum,2);
for i = 1:length(combo)
    pair = combo(i,:);
    p1 = pair(1);
    p2 = pair(2);
    W_(p1,p2) = dtw(data(p1,:),data(p2,:));
end
save('dtw.mat','W_');
W_ = W_ - min(W_(:));
W_ = W_ ./ max(W_(:)); %normalized to (0,1)
W_ = W_(idx,idx);
W_ = 1 - W_; %convert distance to correlation-like
W_ = max(W_, W_'); %symmetric N by N
for run=1:10
    D = diag(sum(W_,2));
%     find(diag(d)==0)
%     fprintf('all zero rows found!\n');
%     pause
    L = D - W_; %unormalized Laplacian
    [evc, evl] = eig(L); %N by N, each column in evc is an eigenvector
    idx = find(diag(evl)>0);
    input = evc(:,idx(1:4));
    c_idx = kmeans(input,4,'Distance','cosine');
    adjrand(c_idx, g_);
    res = [res adjrand(c_idx, g_)];
end
mean(res)
std(res)

%% PCA
clc
load('keti_aligned_nopir.mat')
k = 4;
res = [];
corr = pca(data);
input = data*corr(:,1:k);
for run=1:10
    c_idx = kmeans(input, k);
    ari = adjrand(c_idx, gt_type);
    res = [res ari];
end
mean(res)
std(res)

%% plot each for inspection
sensorNum = size(data,1); % N stream
figure
hold on
for i = 1:sensorNum
    cur = data(i,:); %1 by d
    plot(cur)
    pause
end