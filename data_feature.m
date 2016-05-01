%% SDH BACnet data
clear
clc
path1 = 'C:\Users\dh5gm\Desktop\data\New\SDH\';
fd1 = dir(path1);
num = size(fd1,1);
% vector = zeros(1474,4);
vector = [];
point = cell(1800,1);
ctr = 1;
for n=3:num
    path2 = strcat(path1, fd1(n).name, '\');%path to each room
    fd2 = dir(path2);
    
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
end
