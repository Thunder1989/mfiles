% corrcoef and output one file for all rooms over diff time intervals
clear
clc
% timescale = [1,3,5,7,14,21,28];
timescale = 30;
% path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI/';%path to the data folers
path = '/Users/hdz_1989/Downloads/SDB/KETI/';%path to the data folers
folder = dir(path);
FN = 'interior_mat.txt';
FID = fopen(FN,'w');

for k = 1:length(timescale)
    T = timescale(k);
    th1 = 48*T/6;%6h
    th2 = 48*T/(1/3);%0.5h
    th3 = 48*T/24;%6h
    th4 = 48*T/(24*6);%6h

%     data_length = 9999999;
    data_length = 43200;
    %search for the shortest data length
%     for fn=4:size(folder,1)
%         path1 = strcat(path, folder(fn).name,'/');%path to the csv files
%         file = dir(strcat(path1,'*.csv'));
%         filename = [path1, file(1).name];
%         rd = csvread(filename);
%         ts = rd(:,1);
%         data = rd(:,2);
%         data = data(ts<(ts(1)+T*24*60*60));
%         if length(data)<=data_length
%             data_length = length(data);
%         end
%     end
    imf_aggr = zeros(4*(size(folder,1)-3), data_length);
    imf_aggr1 = zeros(4*(size(folder,1)-3), data_length);
    imf_aggr2 = zeros(4*(size(folder,1)-3), data_length);

    %load data => apply emd => aggregate imfs
    for fn=4:size(folder,1)
        path1 = strcat(path, folder(fn).name,'/');%path to the csv files
        file = dir(strcat(path1,'*.csv'));
        num = size(file,1);
%         figure
        for n = 1:num
            filename = [path1, file(n).name];
            data = csvread(filename);
%             ts = rd(:,1);
%             data = rd(:,2);
%             data = data(ts<(ts(1)+T*24*60*60));
            len = length(data);
            if len > data_length
                [p,q] = rat(data_length/len);
                data = resample(data, p, q);
                data = data(1:data_length);
            end
%             outlier removal
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
%             window = ceil(data_length/(T*24*60)*2);%2min moving window, might be removed later
            window = 10;
            data = EWMA(data', window)';%too much noise in the data, so do a exp. moving average
            fprintf('===========computing the imfs of %s/%s===========\n', folder(fn).name, file(n).name)
            imfs = emd(data);%each row is an IMF
            freq = ZCR(imfs');
%             imf_aggr(((fn-4)*3+n),:) = sum(imfs(freq>th1 & freq<th2,:),1);
%             imf_aggr(((fn-4)*3+n),:) = sum(imfs(freq>th2,:),1);
            imf_aggr(((fn-4)*4+n),:) = sum(imfs(freq>th1 & freq<th2,:),1);
            imf_aggr1(((fn-4)*4+n),:) = sum(imfs(freq>th3 & freq<th2,:),1);
            imf_aggr2(((fn-4)*4+n),:) = sum(imfs(freq>th4 & freq<th2,:),1);

%             subplot(4,1,n)
%             hold on
%             plot(data,'r--')
%             plot(imf_aggr(n,:))
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
    HeatMap(xcor_mat)
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>
    xcor_mat = eye(num);
    for i=1:num
      for j=1:(i-1)
        coef = corrcoef(imf_aggr1(i,:),imf_aggr1(j,:));
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
    HeatMap(xcor_mat)
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>    
    xcor_mat = eye(num);
    for i=1:num
      for j=1:(i-1)
        coef = corrcoef(imf_aggr2(i,:),imf_aggr2(j,:));
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
    HeatMap(xcor_mat)
    
end
fclose(FID);