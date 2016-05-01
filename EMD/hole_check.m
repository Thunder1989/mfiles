%checking holes in data
%% remove hole period from ts
clear
clc
% path = '/Users/hdz_1989/Downloads/SDB/KETI/tmp/';
% path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI/';
path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI_int/';
folder = dir(path);
first = 1371024000;
last = 1371628800;
ts = first:last;
num = size(folder,1);
% figure
data = {};
for n=4:num
%     figure
    path1 = strcat(path, folder(n).name,'/');%path to the csv files
    file = dir(strcat(path1,'*.csv'));
    len=size(file,1);
    for i=1:len
        if i~=4%skip pir
            filename = [path1, file(i).name];
            input = csvread(filename);
            input = input(:,1);
            if input(1)>first
                ts(ts<input(1))=[];
            end
            if input(end)<last
                ts(ts>input(end))=[];
            end

            delta = diff(input);
            idx = (delta>=60*5);
            t_end = input([false idx']);
            dt = delta(idx);
            t_start = t_end - dt;
            for j=1:length(t_start)
                ts(ts>t_start(j) & ts<t_end(j))=[];
            end
        end
%         data{n-2,i} = delta;
%         subplot(len,1,i)
%         grid on
%         plot(delta,'k')
    end
end

%% get the common intervals of timestamp
delta = diff(ts);
idx = (delta>1);
t_end = ts([false idx]);
t_start = t_end - delta(idx);
t = [t_start; t_end];
t = reshape(t,1,numel(t));
limit = [ts(1) t ts(end)];

%% get the shortest data length
%align the readings and resample
% data_length = 999999;
% for n=3:num
% %     figure
%     path1 = strcat(path, folder(n).name,'/');%path to the csv files
%     file = dir(strcat(path1,'*.csv'));
%     len=size(file,1);
%     for i=1:len
%         filename = [path1, file(i).name];
%         input = csvread(filename);
%         ts = input(:,1);
%         input = input(:,2);
%         for j=2:2:length(limit)-1
%             bd1 = limit(j);
%             bd2 = limit(j+1);
%             input(ts>bd1 & ts<bd2) = [];
%         end
%         if length(input)<data_length
%             data_length = length(input);
%         end
%     end
% end

%% align the readings and resample
% clear
% clc
% path = '/Users/hdz_1989/Downloads/SDB/KETI_tmp/';
% folder = dir(path);
data_length = 101243;
len = size(folder,1);
T = 7;
th1 = 48*T/6;%6h
th2 = 48*T/(1/3);%20min
imf_aggr = zeros(4*(size(folder,1)-3), data_length);
FN = '7day_corr.txt';
FID = fopen(FN,'w');

%load data => apply emd => aggregate imfs
for fn=4:len
    path1 = strcat(path, folder(fn).name,'/');%path to the csv files
    file = dir(strcat(path1,'*.csv'));
    num = size(file,1);
    for n = 1:num
        filename = [path1, file(n).name];
        data = csvread(filename);
        ts = data(:,1);
        data = data(:,2);
        for j=2:2:length(limit)-1
            bd1 = limit(j);
            bd2 = limit(j+1);
            data(ts>bd1 & ts<bd2) = [];
        end

        len = length(data);
        if len > data_length+200
            [p,q] = rat(data_length/len);
            data = resample(data, p, q);
            data = data(1:data_length);
        else
            data = data(1:data_length);
        end
        %outlier removal
        if n~=3
            tmp=data(1);
            for k=2:length(data)
                if abs(abs(data(k))-abs(tmp))>3*min(abs(tmp),abs(data(k))) && tmp~=0
                    data(k) = tmp;
                    continue;
                end
                tmp = data(k);
            end
        end
        window = ceil(data_length/(T*24*60)*10);%2min moving window, might be removed later
        data = EWMA(data', window)';%too much noise in the data, so do a exp. moving average
        fprintf('===========computing the imfs of %s/%s===========\n', folder(fn).name, file(n).name)
        imfs = emd(data);%each row is an IMF
        freq = ZCR(imfs');
        imf_aggr(((fn-4)*4+n),:) = sum(imfs(freq>th1 & freq<th2,:),1);
    end
end

%compute xcorr
num = size(imf_aggr,1);
xcor_mat = eye(num);
for i=1:num
  for j=1:(i-1)
    coef = corrcoef(imf_aggr(i,:),imf_aggr(j,:));
    xcor_mat(i, j) = coef(1,2);
    xcor_mat(j, i) = xcor_mat(i, j);
  end
end
%     fprintf(FID, '==========xcormat obtained with %d days data==========\n', T);
for i=1:num
    for j=1:num
        fprintf(FID, '%.2f\t', xcor_mat(i, j));
    end
    fprintf(FID,'\n');
end 
fclose(FID);

%% compute hole ratio to whole trace
% clear
% clc
% path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI_tmp/';
% folder = dir(path);
% first = 1371024000;
% last = 1373673540;
% period = last-first;
% num = size(folder,1);
% ratio = zeros(num-4,5);
% for n=5:num
%     path1 = strcat(path, folder(n).name,'/');%path to the csv files
%     file = dir(strcat(path1,'*.csv'));
%     len=size(file,1);
%     for i=1:len
%         filename = [path1, file(i).name];
%         input = csvread(filename);
%         input = input(:,1);
%     %     input = importdata(filename);
%     %     input = input.data;
%     % %     outlier removal
%     %     if n~=3
%     %         tmp=input(1);
%     %         for k=2:length(input)
%     %             if abs(abs(input(k))-abs(tmp))>3*min(abs(tmp),abs(input(k))) && tmp~=0
%     %                 input(k) = tmp;
%     %                 continue;
%     %             end
%     %             tmp = input(k);
%     %         end
%     %         input = EWMA(input', 50);%too much noise in the data, so do a exp. moving average
%     %     end
%         delta = diff(input);
%         hole_ts = delta(delta>=60*5);
%         ratio(n-4,i) = sum(hole_ts)/period;
%     end
% end