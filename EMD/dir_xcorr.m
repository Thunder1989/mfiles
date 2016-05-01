clear
clc
% timescale = [1,3,5,7,14,21,28];
path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI/';%path to the data folers
folder = dir(path);
FN = 'coef_raw.txt';
FID = fopen(FN,'w');

% for k = 1:length(timescale)
%     T = timescale(k);
%     th1 = 8*T;%6h
%     th2 = 96*T;%0.5h
    data_length = 275770;
    %search for the shortest data length
%     for fn=5:size(folder,1)-1
%         path1 = strcat(path, folder(fn).name,'/');%path to the csv files
%         file = dir(strcat(path1,'*.csv'));
%         filename = [path1, file(1).name];
% %         rd = csvread(filename);
% %         ts = rd(:,1);
% %         data = rd(:,2);
% %         data = data(ts<(ts(1)+T*24*60*60));
%         data = csvread(filename);
%         if length(data)<=data_length
%             data_length = length(data);
%         end
%     end
    data_aggr = zeros(4*(size(folder,1)-5), data_length);

    %load data => direct xocrr as a baseline
    for fn=5:size(folder,1)-1
        path1 = strcat(path, folder(fn).name,'/');%path to the csv files
        file = dir(strcat(path1,'*.csv'));
        num = size(file,1);
        for n = 1:num
            filename = [path1, file(n).name];
%             rd = csvread(filename);
%             ts = rd(:,1);
%             data = rd(:,2);
%             data = data(ts<(ts(1)+T*24*60*60));
            data = csvread(filename);
            len = length(data);
            if len > data_length
                [p,q] = rat(data_length/len);
                data = resample(data, p, q);
                data = data(1:data_length);
            end
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
            data_aggr(((fn-5)*4+n),:) = data;
        end
    end
    
    %compute xcorr
    num = size(data_aggr,1);
    xcor_mat = eye(num);
    for i=1:num
      for j=1:(i-1)
        coef = corrcoef(data_aggr(i,:),data_aggr(j,:));
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
% end
fclose(FID);