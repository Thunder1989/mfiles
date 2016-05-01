%corrcoef and output a file for each room 
clear
clc
path = '/Users/hdz_1989/Downloads/SDB/SDH/';%path to the data folers
folder = dir(path);
for fn=3:size(folder,1)
    path1 = strcat(path, folder(fn).name,'/');%path to the csv files
    file = dir(strcat(path1,'*.csv'));
    num = size(file,1);
    imf_aggr = [];
    tmp = strsplit(path1,'/');
    FN = sprintf('cormat_%s.txt',tmp{end-1});
    FID = fopen(FN,'w');
    timescale = [1,3,5,7,14,21,28];

    for k = 1:length(timescale)
        T = timescale(k);
        th1 = 8*T;%6h
        th2 = 96*T;%0.5h
        for n = 1:num
            filename = [path1, file(n).name];
            rd = csvread(filename);
            if length(rd)>500000
                rd = downsample(rd, ceil(length(rd)/500000));
            end
        %     data = downsample(data, floor(length(data)/4949));
        %     data = data(1:4949);
        %     data = downsample(data, 6);
            ts = rd(:,1);
            data = rd(:,2);
            data = data(ts<(ts(1)+T*24*60*60));
            window = ceil(length(data)/(T*24*60)*2);%2min moving window
            data = EWMA(data', window)';%too much noise in the data, so do a exp. moving average
            fprintf('=============computing the emd of %s/%s=============\n', folder(fn).name, file(n).name)
            imfs = emd(data);%each row is an IMF
            freq = ZCR(imfs');
            imf_aggr(:,n) = sum(imfs([freq>th1 & freq<th2],:),1);
        end

        % figure
        % for i=1:n
        %     subplot(n,1,i)
        %     plot(imf_aggr(:,i))
        % end

        cor_mat = eye(num);
        for i=1:num
          for j=1:(i-1)
            len1 = length(imf_aggr(:,i));
            len2 = length(imf_aggr(:,j));
            if len1 == len2 && len1 > 2
                coef = corrcoef(imf_aggr(:,i),imf_aggr(:,j));
                cor_mat(i, j) = coef(1,2);
            else
                cor_mat(i, j) = 0;
            end
            cor_mat(j, i) = cor_mat(i, j);
          end
        end

        fprintf(FID, '==========cormat of obtained with %d days data==========\n', T);
        for i=1:num
            for j=1:num
                fprintf(FID, '%.2f\t', cor_mat(i, j));
            end
            fprintf(FID,'\n');
        end

        % HeatMap(cor_mat)
        imf_aggr = [];
    end
    fclose(FID);
end