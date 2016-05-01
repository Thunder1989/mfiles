%% tweak the lumin results
clear
clc
data = importdata('/Users/hdz_1989/Documents/Dropbox/SDB/result/KETI/keti_xcormat_nooutlierremovalonlight_MA.txt');
data = data(3:4:end,3:4:end);
data = abs(data);
% label the orientation of room
% 1-N
% 2-E
% 3-S
% 4-W
% 5-interior
% ori = zeros(size(data,1),1);
% ori([1:5,7,18:20,26:30,40:42,44])=1;
% ori([31:33,35,36,48,50])=3;
% ori([17,47])=4;
% ori(ori==0)=5;
N = [1:5,7,18:20,26:30,40:42,44];
S = [31:33,35,36,48,50];
W = [17,47];
seq = [N S W]';
index = 1:length(data);
index = index';%table recording ID of each row/col

%reorder the corr matrix
for i=1:length(seq)
    if index(i)~=seq(i)
        j = find(index==seq(i));
        data([i,j],:) = data([j,i],:);
        data(:,[i,j]) = data(:,[j,i]);
        index([i,j]) = index([j,i]);
    end
end
HeatMap(data);

%% compute the X-X distribution
figure
%N+W, merge the two because the W case is similar to N
len = numel(seq);
tmp = data(:,1:len);
tmp(19:25,:)=[];
tmp(:,19:25)=[];
len = length(N)+length(W);
intra = tmp(1:len,:);
intra(intra==1)=[];
inter = tmp(len+1:end,:);
inter = reshape(inter,1,numel(inter));
[p1,v1] = ecdf(intra);
[p2,v2] = ecdf(inter);
subplot(3,1,1)
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'N-intra','N-inter','Location','SouthEast')

%S
len1 = length(N);
len2 = length(S);
len = len1+len2;
tmp = data(:,len1+1:len);
intra = tmp(len1+1:len,:);
intra(intra==1)=[];
tmp(len1+1:len,:)=[];
inter = tmp;
inter = reshape(inter,1,numel(inter));
[p1,v1] = ecdf(intra);
[p2,v2] = ecdf(inter);
subplot(3,1,2)
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'S-intra','S-inter','Location','SouthEast')

%Interior
len = numel(seq);
tmp = data(:,len+1:end);
intra = tmp(len+1:end,:);
intra(intra==1)=[];
inter = tmp(1:len,:);
inter = reshape(inter,1,numel(inter));
[p1,v1] = ecdf(intra);
[p2,v2] = ecdf(inter);
subplot(3,1,3)
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'I-intra','I-inter','Location','SouthEast')


%% backup
% figure
% cor(cor==1) = [];%remove X-X pair in room on the diagonal
% [p,v] = ecdf(cor);
% subplot(4,1,i)
% plot(v,p,'k--','LineWidth',1.5)
% grid on
