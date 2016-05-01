%% identify broken sensors
data = csvread('C:\Users\dh5gm\Dropbox\SDB\KETI_7\413\co2.csv');
% figure
% hold on
% grid on
% data = input(:,2);
% % subplot(2,1,1)
% plot(data,'r--')
tmp = EWMA(data(:,2)',100);
% % subplot(2,1,2)
% plot(data,'b')

tmp = detrend(tmp);
plot(tmp,'b')
L = length(tmp);
N = 2^nextpow2(L);
D = fft(tmp, N)/L;
figure
plot(2*abs(D(1:N/2+1)))
%% Rice data
clear
clc
path1 = '/Users/hdz/Downloads/code/Data/Rice_new/';
fd1 = dir(path1);
num = size(fd1,1);
% vector = zeros((num-2)*4,4);
% vector = [];
vector = cell(607,14);
point = cell(607,1);
class = zeros(607,1);
ctr =1;
pause on
for n=3:num
    path2 = strcat(path1, fd1(n).name,'/');%path to each room
    files = dir(strcat(path2,'*.csv'));%path 
    k = size(files,1);
    if strcmp(fd1(n).name,'co2')
        id=1;
        continue
%     elseif strcmp(fd1(n).name,'hum')
%         id=2;
%     elseif strcmp(fd1(n).name,'pressure')
%         id=3;
%     elseif strcmp(fd1(n).name,'rmt')
%         id=4;
%     elseif strcmp(fd1(n).name,'status')
%         id=5;
%     elseif strcmp(fd1(n).name,'stpt')
%         id=6;
    elseif strcmp(fd1(n).name,'flow')
        id=7;
%     elseif strcmp(fd1(n).name,'temp_other')
%         id=8;
%     elseif strcmp(fd1(n).name,'occupancy')
%         id=21;
    else
        continue
    end

    for i=1:k
        filename = [path2, files(i).name];
        tmp = csvread(filename);
        tmp = sortrows(tmp,1);
        ts = tmp(:,1); %timestamp
        tmp = tmp(:,2); %data
        figure
        plot(ts,tmp)
%         imfs = emd(tmp');
%         n = size(imfs,1);
%         figure
%         for j=1:n
%             subplot(n,1,j)
%             plot(imfs(j,:))
%         end
        pause
        continue
        if mean(tmp)==0 %|| length(tmp)<700
%             delete(filename)
%             fprintf('length: %d\n', length(tmp));
            continue
        end
%         tmp = EWMA(tmp',20);
        
%new feature extraction
%get from every 45-minute window with 1/3 or none overlap
%and get the stat on these vectors as feature
        mini = [];
        maxi = [];
        medi = [];
        rmsq = [];
        vrn = [];
        skw = [];
        kur = [];
        q1 = [];
        q3 = [];
        qr = [];
        slope = [];
        p_num = [];
        p_dis = [];
        p_amp = [];

%         start = ts(1);
        start = 1370077200;
        %   different time scales:
        %   1. hourly
        %   2. daily
        %   3. diurnal
        while start<ts(end)
            stop = start + 60*60*8;
            if sum(ts>=start&ts<stop)==0
                start = start + 60*60*24;
                continue;
            else
                d = tmp(ts>=start&ts<stop);
                ts_cur = ts(ts>=start&ts<stop);
                start = start + 60*60*24;
            end
            
            mini = [mini min(d)];
            maxi = [maxi max(d)];
            medi = [medi median(d)];
            rmsq = [rmsq rms(d)];
            vrn = [vrn var(d)];
            skw = [skw skewness(d)];
            kur = [kur kurtosis(d)];
            q1 = [q1 quantile(d, 0.25)];
            q3 = [q3 quantile(d, 0.75)];
            qr = [qr iqr(d)];
            p = polyfit(1:length(d),d',1); 
            slope = [slope p(1)];
            try
                [y, idx] = findpeaks(d);
                n = numel(y);
                amp = mean(y);
                dis = mean(diff(ts_cur(idx)));
            catch ME
                n = 0;
                amp = 0;
                dis = -999;
            end
            p_num = [p_num n];
            p_amp = [p_amp amp];
            p_dis = [p_dis dis];

%             pks = [pks numel(findpeaks(d))];
        end
%         new = [];
% %         new = [min(pk) max(pk) median(pk) var(pk)];
%         new = [new min(mdn) max(mdn) median(mdn) var(mdn)];
%         new = [new min(vrn) max(vrn) median(vrn) var(vrn)];
%         new = [new min(range) max(range) median(range) var(range)];
%         new = [new min(pks) max(pks) median(pks) var(pks)];
%         new = [new id];
%         vector = [vector; new];
        vector{ctr,1} = mini;
        vector{ctr,2} = maxi;
        vector{ctr,3} = medi;
        vector{ctr,4} = rmsq;
        vector{ctr,5} = vrn;
        vector{ctr,6} = skw;
        vector{ctr,7} = kur;
        vector{ctr,8} = q1;
        vector{ctr,9} = q3;
        vector{ctr,10} = qr;
        vector{ctr,11} = slope;
        vector{ctr,12} = p_num;
        vector{ctr,13} = p_amp;
        vector{ctr,14} = p_dis;
        point{ctr} = filename;
        class(ctr) = id;
        ctr = ctr+1;
    end
end

minlen = length(vector{1,1});
for i=2:size(vector,1)
    if length(vector{i,1})<minlen && ~isempty(vector{i,1})
        minlen = length(vector{i,1});
    end
end

f = zeros(length(vector),4*14+1);
for i=1:size(f,1)
%     tmp1 = vector{i,1};
%     tmp2 = vector{i,2};
%     tmp3 = vector{i,3};
%     tmp4 = vector{i,4};
%     tmp5 = vector{i,5};
%     tmp6 = vector{i,6};
%     tmp7 = vector{i,7};
%     tmp8 = vector{i,8};
%     tmp9 = vector{i,9};
%     tmp10 = vector{i,10};
%     tmp11 = vector{i,11};  
%     f(i,:) = [tmp1(1:minlen) tmp2(1:minlen) tmp3(1:minlen) tmp4(1:minlen) ...
%         tmp5(1:minlen) tmp6(1:minlen) tmp7(1:minlen) tmp8(1:minlen) tmp9(1:minlen) ...
%         tmp10(1:minlen) tmp11(1:minlen) class(i)];
    tmp = [];    
    for j=1:14
        src = vector{i,j};
        tmp = [tmp min(src) max(src) median(src) var(src)];
    end
    tmp = [tmp class(i)];
    f(i,:) = tmp;
end
f(isnan(f)) = -999;

%% KETI data
clear
clc
path1 = '/Users/hdz/Downloads/code/Data/KETI_oneweek/';
fd1 = dir(path1);
num = size(fd1,1);
% vector = zeros((num-2)*4,4);
vector = [];
point = cell(300,1);
ctr =1;
pause on
for n=3:num
    path2 = strcat(path1, fd1(n).name,'/');%path to each room
    files = dir(strcat(path2,'*.csv'));%path 
    k = size(files,1);
    for i=1:k
        if strcmp(files(i).name,'pir.csv')
            continue
        end
        filename = [path2, files(i).name];
        tmp = csvread(filename);
        ts = tmp(:,1); %timestamp
        tmp = tmp(:,2); %data reading
        if strcmp(files(i).name,'co2.csv')
            id=1;
        elseif strcmp(files(i).name,'humidity.csv')
            id=2;
            continue
        elseif strcmp(files(i).name,'light.csv')
%             id=3;
            continue;
        elseif strcmp(files(i).name,'temperature.csv')
            id=4;
            continue
        end
%         tmp = EWMA(tmp',20);
        
        figure
        plot(ts,tmp)
        pause

%old feature tuple, for the entire stream
%         vector(ctr,4) = id;
%         vector(ctr,3) = std(tmp);
%         vector(ctr,2) = median(tmp);
%         [N,X] = hist(tmp);
%         vector(ctr,1) = X(N==max(N));
%         X(N==max(N))=[];
%         N(N==max(N))=[];
% %         vector(ctr,2) = X(N==max(N));
% %         X(N==max(N))=[];
% %         N(N==max(N))=[];
%         ctr = ctr+1;

%new feature extraction
%get from every 15-minute window with 50% overlap
%and get the stat on these vectors as feature
        %peak1, mean, std
%         pk = [];
        mdn = [];
        vrn = [];
        start = ts(1);
        while start<ts(end)
            stop = start + 45*60;
            d = tmp(ts>=start&ts<stop);
            start = start + 45*60;
            if isempty(d)
                continue
            end
            vrn = [vrn var(d)];
            mdn = [mdn median(d)];
%             [N,X] = hist(d,5);
%             pk = [pk X(N==max(N))];
        end
        new = [];
%         new = [min(pk) max(pk) median(pk) var(pk)];
        new = [min(mdn) max(mdn) median(mdn) var(mdn)];
        new = [new min(vrn) max(vrn) median(vrn) var(vrn)];
        new = [new id];
        vector = [vector; new];
        point{ctr} = filename;
        ctr = ctr+1;
    end
end

% [a, C] = kmeans(vector,4);
%% SDH BACnet data
clear
clc
path1 = '/Users/hdz/Downloads/code/Data/SDH_new/';
fd1 = dir(path1);
num = size(fd1,1);
% vector = zeros(1474,4);
vector = cell(1512,11);
point = cell(1512,1);
class = zeros(1512,1);
ctr = 1;
pause on
for n=3:num
    path2 = strcat(path1, fd1(n).name, '/');%path to each type
    fd2 = dir(path2);
    
%     if strcmp(fd1(n).name,'pos')
%         id=5;
%     elseif strcmp(fd1(n).name,'stpt')
%         id=6;
%     elseif strcmp(fd1(n).name,'rmt')
%         id=4;
    if strcmp(fd1(n).name,'flow')
        id=7;
%     elseif strcmp(fd1(n).name,'other_t')
%         id=8;
%     elseif strcmp(fd1(n).name,'pwr')
%         id=9;
%     elseif strcmp(fd1(n).name,'ctrl')
%         id=10;
%     elseif strcmp(fd1(n).name,'occu')
%         id=11;
%     elseif strcmp(fd1(n).name,'spd')
%         id=12;
%     elseif strcmp(fd1(n).name,'status')
%         id=13;
    else
        continue
    end

    for j=3:size(fd2,1)
        filename = [path2, fd2(j).name];
%         s = dir(filename);
%         if s.bytes == 0
%             fprintf('%s is empty', filename);
%             delete(filename)
%         else
        tmp = csvread(filename); %only read the second col
        ts = tmp(:,1);
        tmp = tmp(:,2);
        %data mean == 0, skip the file
        if mean(tmp)==0 %|| length(tmp)<700
%             delete(filename)
%             fprintf('length: %d\n', length(tmp));
            continue
        end
        %temp convert from F to C
        if strcmp(fd1(n).name,'rmt')
            tmp = (tmp-32)*5/9;
        elseif strcmp(fd1(n).name,'other_t')
            tmp = (tmp-32)*5/9;
        end
           
        figure
        plot(ts,tmp)
%         imfs = emd(tmp');
%         n = size(imfs,1);
%         figure
%         for k=1:n
%             subplot(n,1,k)
%             plot(imfs(k,:))
%         end
        pause

%         if mean(data)==0 && median(data)==0
%             filename
%             delete(filename)
%         end

        mini = [];
        maxi = [];
        medi = [];
        rmsq = [];
        vrn = [];
        skw = [];
        kur = [];
        q1 = [];
        q3 = [];
        qr = [];
        slope = [];
        start = ts(1);
        while start<ts(end)
            stop = start + 60*60;
            if sum(ts>=start&ts<stop)==0
                start = start + 30*60;
                continue;
            else
                d = tmp(ts>=start&ts<stop);
                start = start + 30*60;
            end
            
            mini = [mini min(d)];
            maxi = [maxi max(d)];
            medi = [medi median(d)];
            rmsq = [rmsq rms(d)];
            vrn = [vrn var(d)];
            skw = [skw skewness(d)];
            kur = [kur kurtosis(d)];
            q1 = [q1 quantile(d, 0.25)];
            q3 = [q3 quantile(d, 0.75)];
            qr = [qr iqr(d)];
            p = polyfit(1:length(d),d',1); 
            slope = [slope p(1)];
%             pks = [pks numel(findpeaks(d))];
        end
%         new = [];
% %         new = [min(pk) max(pk) median(pk) var(pk)];
%         new = [new min(mdn) max(mdn) median(mdn) var(mdn)];
%         new = [new min(vrn) max(vrn) median(vrn) var(vrn)];
%         new = [new min(range) max(range) median(range) var(range)];
%         new = [new min(pks) max(pks) median(pks) var(pks)];
%         new = [new id];
%         vector = [vector; new];
        vector{ctr,1} = mini;
        vector{ctr,2} = maxi;
        vector{ctr,3} = medi;
        vector{ctr,4} = rmsq;
        vector{ctr,5} = vrn;
        vector{ctr,6} = skw;
        vector{ctr,7} = kur;
        vector{ctr,8} = q1;
        vector{ctr,9} = q3;
        vector{ctr,10} = qr;
        vector{ctr,11} = slope;
        point{ctr} = filename;
        class(ctr) = id;
        ctr = ctr+1;
    end
end

minlen = length(vector{1,1});
for i=2:size(vector,1)
    if length(vector{i,1})<minlen && ~isempty(vector{i,1})
        minlen = length(vector{i,1});
    end
end
f = zeros(length(vector),minlen*11+1);
for i=1:size(f,1)
    tmp1 = vector{i,1};
    tmp2 = vector{i,2};
    tmp3 = vector{i,3};
    tmp4 = vector{i,4};
    tmp5 = vector{i,5};
    tmp6 = vector{i,6};
    tmp7 = vector{i,7};
    tmp8 = vector{i,8};
    tmp9 = vector{i,9};
    tmp10 = vector{i,10};
    tmp11 = vector{i,11};  
    f(i,:) = [tmp1(1:minlen) tmp2(1:minlen) tmp3(1:minlen) tmp4(1:minlen) ...
        tmp5(1:minlen) tmp6(1:minlen) tmp7(1:minlen) tmp8(1:minlen) tmp9(1:minlen) ...
        tmp10(1:minlen) tmp11(1:minlen) class(i)];
end
f(isnan(f)) = -999;

% color = 'rgbk';
% figure
% hold on
% for i=1:size(vector,1)
%     id = vector(i,4);
%     plot(vector(i,3), vector(i,1), 'color', color(id), 'Line', 'o');
% end

%% Soda data feature
clear
clc
path1 = 'C:\Users\dh5gm\Desktop\data\New\Soda\';
fd1 = dir(path1);
% folder = dir(strcat(path,'*Active*.csv'));
num = size(fd1,1);
vector = [];
point = cell(1400,1);
ctr = 1;
for n=3:num
    
    if strcmp(fd1(n).name,'pos')
        id=5;
    elseif strcmp(fd1(n).name,'stpt')
        id=6;
    elseif strcmp(fd1(n).name,'rmt')
        id=4;
    elseif strcmp(fd1(n).name,'flow')
        id=7;
    elseif strcmp(fd1(n).name,'other_t')
        id=8;
    elseif strcmp(fd1(n).name,'pwr')
        id=9;
    elseif strcmp(fd1(n).name,'ctrl')
        id=10;
    elseif strcmp(fd1(n).name,'occu')
        id=11;
    elseif strcmp(fd1(n).name,'spd')
        id=12;
    elseif strcmp(fd1(n).name,'status')
        id=13;
    elseif strcmp(fd1(n).name,'pressure')
        id=14;
    elseif strcmp(fd1(n).name,'tmr')
        id=15;
    end
    
    path2 = strcat(path1, fd1(n).name, '\');%path to each room
    fd2 = dir(path2);
    for j=3:size(fd2,1)
        filename = [path2, fd2(j).name];
%         s = dir(filename);
%         if s.bytes == 0
%             fprintf('%s is empty', filename);
%             delete(filename)
%         else
        tmp = csvread(filename); %only read the second col
        ts = tmp(:,1);
        tmp = tmp(:,2);
        %data mean == 0, skip the file
        if mean(tmp)==0 %|| length(tmp)<700
%             delete(filename)
%             fprintf('length: %d\n', length(tmp));
            continue
        end
        %temp convert from F to C
        if strcmp(fd1(n).name,'rmt')
            tmp = (tmp-32)*5/9;
        elseif strcmp(fd1(n).name,'other_t')
            tmp = (tmp-32)*5/9;
        end

%         if mean(data)==0 && median(data)==0
%             filename
%             delete(filename)
%         end

%             vector(ctr,4) = id;
%             vector(ctr,3) = std(tmp);
%             vector(ctr,2) = median(tmp);
%             [N,X] = hist(tmp);
%             vector(ctr,1) = X(find(N==max(N),1));
% %             X(find(N==max(N),1))=[];
% %             N(find(N==max(N),1))=[];
%     %         vector(ctr,2) = X(find(N==max(N),1));
%     %         X(find(N==max(N),1))=[];
%     %         N(find(N==max(N),1))=[];
%             ctr = ctr+1;
%     %         end

%new feature extraction
%get from every 15-minute window with 50% overlap
%and get the stat on these vectors as feature
        %peak1, mean, std
%             pk = [];
        mdn = [];
        vrn = [];
        start = ts(1);
        while start<ts(end)
            stop = start + 45*60;
            d = tmp(ts>=start&ts<stop);
            start = start + 45*60;
            if isempty(d)
                continue;
            end
            vrn = [vrn var(d)];
            mdn = [mdn median(d)];
%                 [N,X] = hist(d,5);
%                 pk = [pk X(N==max(N))];
        end

        if isempty(mdn) || isempty(vrn)
            continue
        end
        new = [];
%             new = [min(pk) max(pk) median(pk) var(pk)];
        new = [min(mdn) max(mdn) median(mdn) var(mdn)];
        new = [new min(vrn) max(vrn) median(vrn) var(vrn)];
        new = [new id];
        vector = [vector; new];
        point{ctr} = fd2(j).name;
        ctr = ctr+1;
    end
    
%     file = dir(strcat(path2, '\', '*200906M.dat'));%path 
%     if size(file,1) == 0
%         continue
%     end
% 
%     filename = strcat(path2, '\', file.name);
% %     tmp = load(filename); %only read the second col
%     tmp = tmp(:,end);
%     tmp = tmp(~isnan(tmp));
%     vector(ctr,4) = id;
%     vector(ctr,3) = std(tmp);
%     vector(ctr,2) = (median(tmp)-32)*5/9;
% %     if std(data) > 1
% %         fprintf('stream: %s\n', fd1(n).name);
% %     end
% %         data = EWMA(data',20);
%     [N,X] = hist(tmp);
%     vector(ctr,1) = (X(find(N==max(N),1))-32)*5/9;
%     X(find(N==max(N),1))=[];
%     N(find(N==max(N),1))=[];
% %         vector(ctr,2) = X(find(N==max(N),1));
% %         X(find(N==max(N),1))=[];
% %         N(find(N==max(N),1))=[];
%     ctr = ctr+1;
end

%% plot 3D
figure
hold on
x = vector(:,1)';
y = vector(:,2)';
z = vector(:,3)';
color = 'rgbk';
for i=1:length(x)
    plot3(x(i),y(i),z(i), 'color', color(mod(i-1,4)+1), 'Line', 'o');
%     text(x(i),y(i),z(i), num2str(i),'VerticalAlignment','bottom',...
 %                                                                  'HorizontalAlignment','right');
end
grid on

%% uncomfortable zones
% T-[73,81] in F [22.8, 27.2], H-[19,40]
path1 = 'C:\Users\dh5gm\Dropbox\SDB\KETI_oneweek\';
fd1 = dir(path1);
num = size(fd1,1);
for n=3:3
    path2 = strcat(path1, fd1(n).name,'\');%path to each room
    files = dir(strcat(path2,'*.csv'));%path 
    k = size(files,1);
    i=4;
    filename = [path2, files(i).name];
    tmp = csvread(filename); %only read the second col
    pir_ts = tmp(:,1);
    pir_data = tmp(:,2);
    
    
    i=2;
    filename = [path2, files(i).name];
    tmp = csvread(filename); %only read the second col
    hum_ts = tmp(:,1);
    hum = tmp(:,2);
    i=5;
    filename = [path2, files(i).name];
    tmp = csvread(filename); %only read the second col
    temp_ts = tmp(:,1);
    tmp = tmp(:,2);
    k = min(length(hum),length(tmp));
    ts = temp_ts(1:k);
    flag = zeros(1,k);
    for i=1:k
%         if hum(i)>=19 && hum(i)<=40 && temp(i)>=22.8 && temp(i)<=27.2
        if tmp(i)>=22.8 && tmp(i)<=27.2
            flag(i)=0;
        else
            flag(i)=30;
        end
    end
    figure
    hold on
    plot(ts, flag, 'ko');
    [hAx,hLine1,hLine2] = plotyy(ts, hum(1:k), ts, tmp(1:k));
%     set(hLine2,'LineStyle','--', 'color', 'r');
    grid on
end

%% uncomfortable distribution SDH_BACnet_Old (winter)
% comfort T range [67, 76]F = [19.4, 24.4]C
clear
clc
path = 'C:\Users\dh5gm\Desktop\data\SDH_temp\';
dis_t_stpt = [];
dis_t_cz = [];
dis_stpt_cz = [];
%9am on the first day
percent = zeros(145,3);
ctr = 1;

for n=1:145
    file = [path, num2str(n)];
    if exist(file, 'file')~=2
        continue
    end
    stp = load(file);
    stp = (stp(1,2)-32)*5/9; 
    file_ = [file, '_s'];
    data = load(file_);
    ts = data(:,1);
    tmp = (data(:,2)-32)*5/9;
    t = [];%t reading in wkhr of the current room
    start = 1327654800;
    for i = 1:7
        t = [t tmp(ts>=start & ts<=(start+8*3600))'];
        start = start + 24*3600;
    end
%     fprintf('len of t %d\n', length(t));   
    
    p = 1-sum(abs(t-stp)<=3)/length(t);
    percent(ctr,1) = n;
    percent(ctr,2) = p;
    p = sum(t<19.4|t>24.4)/length(t);
    percent(ctr,3) = p;

    dis_t_stpt = [dis_t_stpt t-stp];
    dis_t_cz = [dis_t_cz zeros(1,sum(t>=19.4&t<=24.4))];
    dis_t_cz = [dis_t_cz (t(t<19.4)-19.4)];
    dis_t_cz = [dis_t_cz (t(t>24.4)-24.4)];
    if stp>=19.4 && stp<=24.4
        dis_stpt_cz = [dis_stpt_cz 0];
    elseif stp<19.4
        dis_stpt_cz = [dis_stpt_cz stp-19.4];
    else
        dis_stpt_cz = [dis_stpt_cz stp-24.4];
    end
    ctr = ctr+1;
end
figure
hold on
grid on
% [p, v] = ecdf(dis_t_stpt);
% plot(v,p,'g', 'LineWidth', 2)
% [p, v] = ecdf(dis_t_cz);
% plot(v,p,'y', 'LineWidth', 2)
[p, v] = ecdf(dis_stpt_cz);
plot(v,p,'k', 'LineWidth', 2)

%% uncomfort distribution SODA (summer)
clear
clc
path1 = 'C:\Users\dh5gm\Desktop\data\Soda\art\';
path1_ = 'C:\Users\dh5gm\Desktop\data\Soda\stpt\';
fd1 = dir(path1);
num = size(fd1,1);
dis_t_stpt = [];
dis_t_cz = [];
dis_stpt_cz = [];
percent = [];
p2 = [];
room = [];

for n=3:num
    path2 = [path1, fd1(n).name];
    file = dir(strcat(path2, '\', '*200906M.dat'));%path 
    if size(file,1) == 0
        continue
    end
    fd2 = strrep(fd1(n).name,'T','S');
    file_ = dir(strcat(path1_, fd2, '\', '*200906M.dat'));%path 
    if size(file_,1) == 0
        continue
    end
    filename = strcat(path2, '\', file.name);
    filename_ = strcat(path1_, fd2, '\', file_.name);
    data = load(filename);
    ts = data(:,1);
    tmp = data(:,end);
    stp = load(filename_);
    stp = stp(:,end);
    stp = stp(~isnan(stp));
    stp = (stp(1)-32)*5/9;
    t = [];
    start = 1243846800;
    for i = 1:29
        var = tmp(ts>=start & ts<=(start+8*3600))';
        var = var(~isnan(var));
        var = (var-32)*5/9;
        t = [t var];
        start = start + 24*3600;
    end
    
    if min(t'-stp)<-20 || max(t'-stp)>20
        continue
    end
    
%     fprintf('len of t %d\n', length(t));
    p = 1-sum(abs(t-stp)<=3)/length(t);
    percent = [percent; p];
    p = sum(t<22.8|t>27.2)/length(t);
    p2 = [p2; p];
    room = [room; fd1(n).name];
    dis_t_stpt = [dis_t_stpt t-stp];
    dis_t_cz = [dis_t_cz zeros(1,sum(t>=22.8&t<=27.2))];
    dis_t_cz = [dis_t_cz (t(t<22.8)-22.8)];
    dis_t_cz = [dis_t_cz (t(t>27.2)-27.2)];
    if stp>=22.8 && stp<=27.2
        dis_stpt_cz = [dis_stpt_cz 0];
    elseif stp<22.8
        dis_stpt_cz = [dis_stpt_cz stp-22.8];
    else
        dis_stpt_cz = [dis_stpt_cz stp-27.2];
    end
end
figure
hold on
grid on
[p, v] = ecdf(dis_t_stpt);
plot(v,p,'r', 'LineWidth', 2)
[p, v] = ecdf(dis_t_cz);
plot(v,p,'r--', 'LineWidth', 2)
% [p, v] = ecdf(dis_stpt_cz);
% plot(v,p,'k', 'LineWidth', 2)

%% occupied uncomfort
clear
clc
path1 = 'C:\Users\dh5gm\Dropbox\SDB\KETI_oneweek\';
fd1 = dir(path1);
num = size(fd1,1);
dis = [];
for n=3:num
    path2 = strcat(path1, fd1(n).name,'\');%path to each room
    files = dir(strcat(path2,'*.csv'));%path 
    k = size(files,1);
    
    i=4;%4-pir, find the occupied periods, coarse
    filename = [path2, files(i).name];
    in = csvread(filename); %only read the second col
    pir_ts = in(:,1);
    pir_data = in(:,2);
    o_time = pir_ts(pir_data>10);
    
    i=5;%5-temp
    filename = [path2, files(i).name];
    in = csvread(filename); %only read the second col
    temp_ts = in(:,1);
    temp_data = in(:,2);
    for j=1:length(o_time)
        start = o_time(j);
        stop = start+5;
        tmp = temp_data(temp_ts>=start&temp_ts<=stop);
        dis = [dis zeros(1,sum(tmp>=22.8&tmp<=27.2))];
        dis = [dis (22.8-tmp(tmp<22.8))'];
        dis = [dis (tmp(tmp>27.2)-27.2)'];
    end
    
end
figure
[p, v] = ecdf(dis);
plot(v , p, 'b', 'LineWidth', 2)

%% bar
figure
% subplot(1,2,1)
% data = [232, 232, 232, 0.1, 0.1; 243, 243, 243, 0.1, 0.1];%soda-50
%data = [162,176,162,0.1,14;277,304,277,0.1,27]; %sdh-50
%data = [196,224,196,0,28;99,174,98,1,76];%soda_grep
% data = [57,18,11;2,1,1;63,9,1;66,6,0;2,0,0]; %# of rogue rooms in soda by zone
data = [10,2,1;21,6,1;21,8,1;21,8,0;21,10,0;20,3,0;16,9,0];
h = bar(data,'grouped');
set(h,'BarWidth',0.6); % The bars will now touch each other
set(gca,'YGrid','on')
%set(gca,'GridLineStyle','+')
%set(gca,'XTicklabel','LW>=sD&&LH>=sD|LW>=sD&&LH<=sD|LW<=sD&&LH>=sD|LW<=sD&&LH<=sD')
set(gca,'XTicklabel','Zone1|Zone2|Zone3|Zone4|Zone5|Zone6|Zone7');
set(get(gca,'YLabel'),'String','# of Rooms');
lh = legend('# of Room','Once >3','80% >3');
set(lh,'Location','NorthEast','Orientation','vertical');
colormap(summer)
% for i=1:length(data)
%     label = text(i,data(i),num2str(data(i))) ;
% end
% set(label,'Horizontalalignment','center',...
% 'verticalalignment','bottom') ;

%% plot example on curves
clear
clc
key = {1.,   2.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,  12., 13.,  17.,  19.,  21.};
value = 1:numel(key);
mapping = containers.Map(key,value); 
precision = [
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 1.0, 1.0, 0.6666666666666666];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.0, 0.8709677419354839, 0.7105263157894737, 0.8709677419354839, 0.5, 0.6428571428571429, 0.9642857142857143, 0.9285714285714286, 0.6410256410256411, 0.8888888888888888, 0.6944444444444444, 0.8709677419354839, 0.7941176470588235, 0.6216216216216216, 0.7352941176470589, 0.875, 0.88, 0.9615384615384616, 0.92, 0.9545454545454546];
[0.0, 0.32217573221757323, 0.30739299610894943, 1.0, 0.7894736842105263, 0.8181818181818182, 0.8333333333333334, 0.8611111111111112, 0.41361256544502617, 0.7238095238095238, 0.5347222222222222, 0.7623762376237624, 0.8571428571428571, 0.8641975308641975, 0.7087378640776699, 0.8131868131868132, 0.881578947368421, 0.7058823529411765, 0.8488372093023255, 0.9102564102564102, 0.9342105263157895, 0.948051948051948, 0.872093023255814, 0.9538461538461539, 0.948051948051948, 0.8860759493670886, 0.9594594594594594, 0.971830985915493, 0.8674698795180723, 0.8860759493670886];
[0.0, 0.0, 0.625, 0.1111111111111111, 0.3023255813953488, 0.35714285714285715, 0.35135135135135137, 0.3333333333333333, 0.2777777777777778, 0.38095238095238093, 0.2923076923076923, 0.38333333333333336, 0.46875, 0.4594594594594595, 0.3333333333333333, 0.4807692307692308, 0.5405405405405406, 0.7272727272727273, 0.7307692307692307, 0.8823529411764706, 0.8823529411764706, 0.92, 0.6666666666666666, 0.8709677419354839, 0.8529411764705882, 0.92, 0.9166666666666666, 1.0, 0.8260869565217391, 1.0];
[0.24390243902439024, 0.6875, 1.0, 1.0, 1.0, 0.9444444444444444, 1.0, 1.0, 0.9090909090909091, 0.9655172413793104, 1.0, 1.0, 0.825, 0.825, 0.8333333333333334, 0.8125, 0.7317073170731707, 0.6590909090909091, 0.7391304347826086, 0.6407766990291263, 0.6704545454545454, 0.7631578947368421, 0.8909090909090909, 0.7532467532467533, 0.775, 0.6407766990291263, 0.6470588235294118, 0.6571428571428571, 0.6288659793814433, 0.6261682242990654];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 0.0, 0.0, 0.0661764705882353, 0.08633093525179857, 0.09285714285714286, 0.1015625, 0.3, 0.14285714285714285, 0.19607843137254902, 0.11428571428571428, 0.14457831325301204, 0.18333333333333332, 0.38461538461538464, 0.3793103448275862, 0.22857142857142856, 0.5384615384615384, 0.32142857142857145, 0.35135135135135137, 0.3, 0.19230769230769232, 0.21818181818181817, 0.23529411764705882, 0.2972972972972973, 0.32432432432432434, 0.4, 0.40625, 0.2727272727272727, 0.3];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2608695652173913, 0.2608695652173913, 0.2857142857142857, 0.2692307692307692, 0.3, 0.2916666666666667, 0.5, 0.3888888888888889, 0.5384615384615384, 1.0, 0.42857142857142855, 0.5, 1.0, 0.46153846153846156, 0.45454545454545453, 0.2777777777777778, 0.42857142857142855];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.38461538461538464, 0.4];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.1, 0.06666666666666667, 0.2, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
];
ex = [  6.,   4.,   5.,   5.,   8.,   4.,   6.,   5.,   4.,   2.,  21.,...
         2.,   6.,   9.,   1.,   5.,   6.,  21.,   8.,   6.,   6.,   5.,...
         6.,   6.,   5.,   6.,  12.,   6.,  10.,   6.];
color = colormap(jet(size(precision,1)));
hold on
for i=1:size(precision,1)
    plot(precision(i,:), 'color', color(i,:), 'LineStyle', '--', 'LineWidth', 2);
end
legend('co2','humidity','rmt','staus','stpt','flow','HW sup','HW ret','CW sup','CW ret','Air sup','Air ret','Air mix',',Occup')

for i=1:numel(ex)
    plot(i, precision(mapping(ex(i)),i), 'color', color(mapping(ex(i)),:), 'Line', 'o', 'LineWidth', 2);
%     text(x, acc(mapping(ex(i)),x), num2str(ex(i+1)),'VerticalAlignment','bottom','HorizontalAlignment','center');
end
% x = 2;
% for i=3:4:numel(ex)-4
%     plot(x, precision(mapping(ex(i)),x), 'color', color(mapping(ex(i+1)),:), 'Line', 'o', 'LineWidth', 2);
% %     text(x, acc(mapping(ex(i)),x), num2str(ex(i+1)),'VerticalAlignment','bottom','HorizontalAlignment','center');
%     x = x+1;
% end

%%
figure
recall = [
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6455696202531646, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.2571428571428571, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.0, 0.0, 1.0, 0.23333333333333334, 0.14285714285714285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.23333333333333334, 0.3, 0.0, 0.0, 0.75, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.0, 1.0, 1.0, 0.43333333333333335, 0.32857142857142857, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.875, 0.9259259259259259, 0.9493670886075949, 0.26666666666666666, 0.24285714285714285, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.5, 0.8148148148148148, 0.9493670886075949, 0.36666666666666664, 0.5142857142857142, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.75, 0.37037037037037035, 0.9240506329113924, 0.5333333333333333, 0.4857142857142857, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.125, 0.7037037037037037, 0.9873417721518988, 0.5, 0.34285714285714286, 0.0, 0.0, 0.75, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.375, 0.9629629629629629, 0.9620253164556962, 0.3333333333333333, 0.6857142857142857, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.75, 0.7777777777777778, 0.9746835443037974, 0.5333333333333333, 0.5714285714285714, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.875, 0.7777777777777778, 0.9746835443037974, 0.3, 0.42857142857142855, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.5, 0.7407407407407407, 0.9746835443037974, 0.5, 0.5857142857142857, 0.0, 0.0, 0.75, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.7777777777777778, 0.9746835443037974, 0.4666666666666667, 0.5714285714285714, 0.0, 0.0, 0.75, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0];
[0.625, 0.7037037037037037, 0.9873417721518988, 0.5333333333333333, 0.5857142857142857, 0.0, 0.0, 0.875, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.875, 0.8148148148148148, 0.9746835443037974, 0.5666666666666667, 0.5571428571428572, 0.0, 0.0, 0.5, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
[0.5, 0.9259259259259259, 0.9746835443037974, 0.4666666666666667, 0.6142857142857143, 0.0, 0.0, 0.875, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.4166666666666667, 0.0, 0.9259259259259259, 0.9746835443037974, 0.8333333333333334, 0.5428571428571428, 0.0, 0.0, 0.75, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.16666666666666666];
[0.75, 0.9629629629629629, 0.9620253164556962, 0.7666666666666667, 0.6142857142857143, 0.0, 0.0, 0.875, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.4166666666666667, 0.5, 0.9629629629629629, 0.9620253164556962, 0.9333333333333333, 0.5714285714285714, 0.0, 0.0, 0.875, 0.6666666666666666, 0.0, 0.5, 0.0, 0.0, 0.0, 0.75];
[0.5, 0.9259259259259259, 0.9746835443037974, 0.9333333333333333, 0.5857142857142857, 0.0, 0.0, 0.75, 0.6666666666666666, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.25, 0.5, 0.8888888888888888, 0.9620253164556962, 0.8333333333333334, 0.7285714285714285, 0.0, 0.0, 0.5, 0.6666666666666666, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.16666666666666666];
[0.5, 0.8518518518518519, 0.9746835443037974, 0.7333333333333333, 0.9142857142857143, 0.0, 0.0, 0.75, 0.6666666666666666, 0.0, 0.5, 0.0, 1.0, 0.0, 0.6666666666666666, 0.5, 1.0, 0.9493670886075949, 0.9, 0.7857142857142857, 0.0, 0.07692307692307693, 0.875, 0.6666666666666666, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.5];
[0.25, 1.0, 0.9746835443037974, 0.7333333333333333, 0.8857142857142857, 0.0, 0.07692307692307693, 0.75, 0.5555555555555556, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.08333333333333333, 0.5, 1.0, 0.8860759493670886, 0.8333333333333334, 0.8571428571428571, 0.0, 0.0, 0.75, 0.6666666666666666, 0.0, 0.4, 0.0, 0.6666666666666666, 0.0, 0.3333333333333333];
[0.5, 1.0, 0.9746835443037974, 0.8333333333333334, 0.7571428571428571, 0.0, 0.07692307692307693, 0.75, 0.4444444444444444, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.25, 0.25, 0.9629629629629629, 0.9873417721518988, 0.8666666666666667, 0.7, 0.0, 0.0, 0.625, 0.5555555555555556, 0.0, 0.5, 0.0, 0.6666666666666666, 0.0, 0.3333333333333333];
];
color = colormap(lines(size(recall,1)));
hold on
for i=1:size(recall,1)
    plot(recall(i,:), 'color', color(i,:), 'LineStyle', '--', 'LineWidth', 2);
end
legend('co2','humidity','rmt','staus','stpt','flow')
x = 2;
for i=3:4:numel(ex)-4
    plot(x, recall(mapping(ex(i)),x), 'color', color(mapping(ex(i+1)),:), 'Line', 'o', 'LineWidth', 2);
%     text(x, acc(mapping(ex(i)),x), num2str(ex(i+1)),'VerticalAlignment','bottom','HorizontalAlignment','center');
    x = x+1;
end

%%
avSpots = tmp;
meanCycle = 1;
Fs = 1;
Nf = 512;

df = Fs/Nf;
f = 0:df:Fs/2-df;

trSpots = fftshift(fft(avSpots-mean(avSpots),Nf));

dBspots = 20*log10(abs(trSpots(Nf/2+1:Nf)));

yaxis = [20 85];
plot(f,dBspots,1./[meanCycle meanCycle],yaxis)
xlabel('Frequency (year^{-1})')
ylabel('| FFT | (dB)')
axis([0 1/2 yaxis])
text(1/meanCycle + .02,25,['<== 1/' num2str(meanCycle)])

%%
k = length(point);
N = 672;
for i=1:k
        tmp = point{i};
        if length(tmp) > N
            [p,q] = rat(N/length(tmp));
            tmp = resample(tmp, p, q);
            tmp = tmp(1:N);
            point{i} = tmp;
        end
end

%% U-shapelet features
clear
clc
vector = dlmread('resampled-all');
label = vector(:,1);
vector = vector(:,2:end);
load('resampled-all_Many_30.mat');
idx = uShapelets(:,1);
start = uShapelets(:,2)+1;
N = size(vector,1);
k = size(vector,2);
f = zeros(N,length(idx)*4);

for i=1:length(idx)
    shplt = vector(idx(i),start(i)+29);
    dist = zeros(k-29,1);
    for j=1:N
        cur = vector(j,:);
        for k=1:k-29
            tmp = cur(k:k+29);
            dist(k) = norm(shplt-tmp);
        end
        f(j,i*4-3) = min(dist);
        f(j,i*4-2) = max(dist);
        f(j,i*4-1) = median(dist);
        f(j,i*4) = var(dist);
    end
end
dlmwrite('shape_all',f);
dlmwrite('label_all',label);


%%
clear
clc
s = cell(16,1); %set of supervised shapelets
s{1} = [-1.110400 -1.097400 -0.974470 -0.950070 -0.920390 ];
s{2} = [1.544300  1.137400  0.876200  0.667750  0.455820  0.304190  0.510980  0.614810  0.614490  0.607420  0.575680  0.516220  0.479850  0.440850  0.416590  0.320770  0.439820  0.546840  0.582060  0.548030  0.484430  0.453690  0.412900  0.367570  0.287700  0.272380  0.200690  0.125910  0.073337  0.029085  0.020791 -0.008006 -0.004427 -0.052863 -0.016465 -0.057144 -0.015232  0.045603  0.115410 -0.107030 -0.645540 -0.833200 -0.927020 -0.850220 -0.512790 -0.323310 -0.222530 -0.087423 -0.183690 -0.352680 -0.397140 -0.444300 -0.450050 -0.436020 -0.464510 -0.422580 -0.425020 -0.442560 -0.418610 -0.377490 -0.308590 -0.319160 -0.259640 -0.240830 -0.190520 -0.183000 -0.172800 -0.160840 -0.140160 -0.157600 ];
s{3} = [0.000000  0.000000  0.000000  0.000000  0.000000 ];
s{4} = [-0.754330 -0.858930 -0.915070 -0.952220 -0.968200 -1.019500 -1.001700 -1.018000 -0.997120 -1.043900 -1.041700 -1.029400 -1.071200 -1.090500 -1.131300 -1.147900 -1.195800 -1.182300 -1.189600 -1.185700 -1.187200 -1.187500 -1.185300 -1.189400 -1.182200 -1.196500 -1.133900 -1.051100 -1.100800 -1.149200 -1.244300 -1.340500 -1.388800 -1.459000 -1.483700 -1.504400 -1.580500 -1.600300 -1.607000 -1.569700 ];
s{5} = [0.392510  0.400640  0.453300  0.558460  0.610820  0.858380  1.121000  1.341400  1.376700  1.386800  1.403000  1.393900  1.401100  1.376800  1.345800  1.342600  1.306800  1.194500  1.147500  1.153200  1.087600  0.992450  0.907660  0.697840  0.718810  0.685750  0.482690  0.345450  0.196890  0.078073 -0.195410 -0.378340 -0.479900 -0.596370 -0.872610 -0.260360  0.146480  0.239550  0.346190  0.388970  0.436700  0.420930  0.444300  0.385190  0.427260 -0.429140 -1.006900 -1.395100 -1.764900 -2.018000 -2.292300 -1.889900 -1.526900 -1.283600 -1.075300 ];
s{6} = [-1.063000 -1.160900 -1.100000 -1.015800 -1.067500 -1.085300 -1.092700 -1.105500 -1.087200 -1.141900 -1.159400 -1.081700 -1.120200 -1.048700 -1.034600 -1.005900 -0.937860 -1.003200 -1.032400 -1.041800 -0.993430 -1.000300 -1.034700 -1.024000 -0.984020 -0.909320 -0.881820 -0.777160 -0.807530 -0.813690 -0.743840 -0.797970 -0.719730 -0.713950 -0.625230 -0.593200 -0.672240 -0.565380 -0.479340 -0.633660 -0.714070 -0.687860 -0.730030 -0.704220 -0.724260 -0.657740 -0.659360 -0.730640 -0.670250 -0.702280 -0.658450 -0.678600 -0.717930 -0.690670 -0.714310 ];
s{7} = [1.882900  1.923300  1.904000  1.883000  1.840000  1.861300  1.728400  1.485700  1.355200  1.280400  1.181800  1.164800  1.115300  0.936170  0.771670  0.640750  0.544390  0.403900  0.311340  0.177120  0.131960 -0.004227 -0.127620 -0.180710 -0.269680 -0.365610 -0.489420 -0.552160 -0.697100 -0.799040 -0.890310 -0.979540 -1.074000 -1.170500 -1.322400 -1.424800 -1.517500 -1.630900 -1.712400 -1.788000 -1.871200 -1.935100 -1.995900 -2.041100 -2.125900 -2.198900 -2.252600 -2.327600 -2.399300 -2.498000 -2.582900 -2.635800 -2.679000 -2.754800 -2.802300 -2.836000 -2.879800 -2.931200 -2.863400 -2.751400 -2.722300 -2.628700 -2.552400 -2.496700 -2.399900 -2.365200 -2.294500 -2.219100 -2.152700 -2.063900 -2.067300 -2.097800 -2.139600 -2.118500 -2.138800 -2.129500 -2.035800 -2.013400 -1.932300 -1.910700 ];
s{8} = [-0.039810  0.118610  0.919930  1.717700 -0.565630 -0.294410  0.014248  1.080500  1.652800 -0.979040 -0.644100  1.062500  1.707400  1.684000 -0.958220 -1.258300  0.538720  1.549200  1.648300 -1.284000 -0.864160  0.948510  1.649800  1.219300 -1.523900 -0.911640  0.143550  0.522090  0.382040 -1.921500 -1.068900 -0.122550  0.063100 -0.261210 -1.997500 -1.037500 -0.666000 -0.251450 -0.822340 -2.431200 -1.133000 -1.046000 -0.311940 -1.138800 -2.504700 -1.754900 -1.008400 -0.662240 -1.740400 -2.566100 -1.315900 -1.007200 -0.857670 -2.146700 -2.191700 -1.367400 -1.220900 -0.823230 -2.379300 -2.233600 -1.402600 -0.970600 -0.742060 -2.465800 -1.982300 -1.365200 -0.914270 -0.605510 -2.379800 -1.362500 ];
s{9} = [1.183300  1.198200  1.219600  1.211100  1.202500  1.232300  1.215700  1.231300  1.249600  1.254700  1.277800  1.294700  1.312800  1.228300  1.152700  1.152400  1.152500  1.078400  1.003600  0.895020  0.808950  0.729120  0.677390  0.596930  0.516690  0.424900  0.345010  0.252570  0.174710  0.206220  0.232140  0.215740  0.226590  0.171690  0.079237 -0.042860 -0.167450 -0.283550 -0.433180 -0.563410 -0.689920 -0.815780 -0.962610 -1.179400 -1.391700 -1.524600 -1.618900 -1.693900 -1.765100 -1.799400 -1.864000 -1.907000 -1.948100 -1.994400 -2.046600 -2.117600 -2.154700 -2.180300 -2.263600 -2.276200 -2.321000 -2.359100 -2.406200 -2.446600 -2.535500 ];
s{10} = [-0.711680 -0.716480 -0.713940 -0.715170 -0.724590 -0.719900 -0.715910 -0.718400 -0.712310 -0.720640 -0.717160 -0.719040 -0.671290 -0.659130 -0.655230 -0.656990 -0.664430 -0.608990 -0.645420 -0.683630 -0.686220 -0.638890 -0.624030 -0.658170 -0.689910 -0.692370 -0.692950 -0.717740 -0.726280 -0.722090 -0.728960 -0.723630 -0.724730 -0.753590 -0.758580 -0.783960 -0.794230 -0.785060 -0.808100 -0.816350 -0.843110 -0.848260 -0.846500 -0.847190 -0.834980 -0.837330 -0.834270 -0.833730 -0.841260 -0.826750 ];
s{11} = [-1.441800 -1.440000 -1.377100 -1.381400 -1.438100 -1.402300 -1.351600 -1.445800 -1.469100 -1.227500 -1.213700 -1.146000 -1.339500 -1.363000 -1.421200 -1.227400 -1.103000 -1.017400 -1.350400 -1.168500 -1.162600 -1.348000 -1.232500 -1.130700 -1.255000 -1.168300 -1.038700 -1.070700 -1.067400 -0.947980 -0.686590 -0.941040 -0.834840 -0.888220 -0.797440 -1.068400 -0.947040 -0.858980 -0.976530 -0.794060 -0.828790 -0.834160 -0.700440 -1.045100 -0.874140 -1.057100 -0.913990 -0.846690 -0.962850 -1.024900 -0.989780 -0.759710 -0.761330 -0.813340 -0.927360 -0.988720 -0.907220 -0.875180 -0.699690 -0.941800 -0.961670 -1.010900 -1.069800 -1.036400 -1.144000 -0.880090 -0.873290 -0.989280 -1.146800 -1.168800 -1.219000 -1.052400 -0.998300 -1.223700 -1.164900 -1.326000 -1.172500 -1.148200 -1.319000 -1.194000 -1.180500 -1.259300 -1.264200 -1.292000 -1.203200 -0.879910 -1.081800 -1.321400 -1.144400 -1.138200 ];
s{12} = [-0.214850 -0.175090 -0.118340 -0.150820 -0.097759 -0.070602 -0.064485  0.015643 -0.003138 -0.039398 -0.004753  0.069590  0.000877  0.061330  0.032019  0.226620 -0.323940  0.289670  0.890190  0.797980  0.984430  1.026600  1.116900  1.143200  1.206000  1.243600  1.246100  1.209100  1.276100  0.993860 ];
s{13} = [-0.137020 -0.156560 -0.153770 -0.153060 -0.129190 -0.147110 -0.152910 -0.134750 -0.139220 -0.155340 -0.197320 -0.231590 -0.206050 -0.216880 -0.258010 -0.258230 -0.297450 -0.287230 -0.285770 -0.291160 -0.334780 -0.348130 -0.352160 -0.396130 -0.415400 -0.463830 -0.457390 -0.501310 -0.432020 -0.564080 -0.206260 -0.045618 -0.502950 -0.454690 -0.505470 -0.446760 -0.488320 -0.475740 -0.452470 -0.385130 -0.305200 -0.241030 -0.199840 -0.129280 -0.063838 ];
s{14} = [-0.302700 -0.270320 -0.150570  0.069613 -0.055524  0.208840 -0.155820  0.043855  0.340160  0.004372 -0.038700  0.069849  0.229020 -0.052429  0.117010  0.196670  0.365340  0.635830  0.314480  0.362750  0.459890  0.417820  0.218850  0.439500  0.257120  0.414730  0.260190  0.110940  0.380230  0.201200  0.148900  0.132130  0.249690  0.291430 -0.039917 -0.046893 -0.152600 -0.426650  0.041280 -0.305900 -0.300290 -0.332910 -0.203530 -0.461990 -0.582790 -0.356210 -0.390920 -0.219370 -0.298520 -0.464150 -0.359220 -0.392690 -0.655560 -0.427500 -0.567610 ];
s{15} = [2.315800  2.229200  2.158100  2.107800  1.993200  1.754700  1.564200  1.480600  1.378200  1.184700  1.091600  1.074000  0.910600  0.987030  0.846480  0.720030  0.687010  0.583970  0.522990  0.423790  0.286800  0.370050  0.498850  0.422640  0.381480  0.307680  0.301280  0.263070  0.169890  0.162590  0.126120  0.145100  0.051232  0.036983  0.051922  0.008900  0.005928 -0.078171 -0.047808 -0.089263 -0.035129 -0.011481 -0.047963 -0.002632  0.310250  0.444670  0.770540  1.076100  1.015300  1.162600  0.985580  0.645220  0.608020  0.510390  0.393480  0.346360  0.258550  0.233790  0.168600  0.119440  0.081345  0.059727  0.044918  0.000336  0.079277  0.261220  0.389000  0.479150  0.656180  0.925950  1.118900  1.350900  1.554700  1.879000  2.518500  2.920200  3.365600  3.042700  2.665200  2.574900  2.489600  2.666200  2.568700  2.336600  2.043700  1.819200  1.558500  1.427400  1.569000  1.672600  1.646400  1.643900  1.612600  1.574600  1.557900  1.599200  1.624100  1.703100  1.761500  1.805700 ];
s{16} = [0.022776  0.334020  0.597060  0.800590  1.016800  1.107700  1.219000  1.326800  1.420500  1.496800  1.526000  1.559100  1.587200  1.589200  1.618800 ];

vector = dlmread('resampled-all');
label = vector(:,1);
vector = vector(:,2:end);
N = size(vector,1);
M = size(vector,2);
f = zeros(N,length(s)*4);

for i=1:length(s)
    shplt = s{i};
    offset = length(shplt)-1;
    dist = zeros(M-offset,1);
    for j=1:N
        cur = vector(j,:);
        for k=1:M-offset
            tmp = cur(k:k+offset);
            dist(k) = norm(shplt-tmp);
        end
        f(j,i*4-3) = min(dist);
        f(j,i*4-2) = max(dist);
        f(j,i*4-1) = median(dist);
        f(j,i*4) = var(dist);
    end
end
dlmwrite('shape_all',f);
dlmwrite('label_all',label);