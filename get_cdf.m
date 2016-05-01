clear
clc
% data = importdata('/Users/hdz_1989/Documents/MATLAB/CEEMDAN/dir_xcormat.txt');
data = importdata('/Users/hdz_1989/Documents/Dropbox/SDB/result/KETI/keti_xcormat_nooutlierremovalonlight_MA.txt');
figure
% room_num = 9;
% % list=[1,3,5,7,14,21,28]
% num=size(data,1)/size(data,2);
% for n=1:room_num%num of rooms
%     intra = [];
%     inter = [];
%     ox=3*n-3;
%     for i=3:num%num of different length days
% %     for i=7:7%num of different length days
%         oy=27*i-27;
%         d=data((oy+1):(oy+27),(ox+1):(ox+3));
%         
%         tmp = d((ox+1):(ox+3),:);
%         tmp(tmp==1) = [];
%         intra = [intra tmp];
%         
%         d((ox+1):(ox+3),:) = [];
%         inter = [inter reshape(d,1,numel(d))];
%     end
%     intra = abs(intra);
%     inter = abs(inter);
%     [p1,v1] = ecdf(intra);
%     [p2,v2] = ecdf(inter);
%     subplot(room_num,1,n)
%     plot(v1,p1,'r',v2,p2,'b*')
% %     plot(v2,p2,'b*')
% end

% %for KETI corrcoef
room = size(data,1)/4;
intra = [];
inter = [];
for n=1:room%num of rooms
    ox=4*n-3;
    d=data(:,ox:ox+3);

    tmp = d(ox:(ox+3),:);
    tmp(tmp==1) = [];
    intra = [intra tmp];

    d(ox:(ox+3),:) = [];
    inter = [inter reshape(d,1,numel(d))];
%     intra = abs(intra);
%     inter = abs(inter);
%     [p1,v1] = ecdf(intra);
%     [p2,v2] = ecdf(inter);
%     subplot(room,1,n)
%     plot(v1,p1,'r',v2,p2,'b*')
%     plot(v2,p2,'b*')
end
intra = abs(intra);
inter = abs(inter);
[p1,v1] = ecdf(intra);
[p2,v2] = ecdf(inter);
plot(v1,p1,'r',v2,p2,'b*')

%for KETI, to see the type pattern across rooms
% for n=1:4
%     intra = data(n:4:end, n:4:end);
%     intra(intra==1) = [];
%     intra = abs(intra);
%     [p1,v1] = ecdf(intra);
%     subplot(4,1,n)
%     plot(v1,p1)
% end