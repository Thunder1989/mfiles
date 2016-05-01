clear
clc
% data = importdata('/Users/hdz_1989/Documents/Dropbox/sdb/result/KETI/keti_xcormat_nooutlierremovalonlight_MA.txt');
% data = importdata('/Users/hdz_1989/Documents/Dropbox/sdb/result/KETI/keti_xcormat_diff_midlen.txt');
% data = data(1:204,:);
% room = length(data)/4;
% tic;

data = importdata('/Users/hdz_1989/Documents/Dropbox/SDB/result/buildsys_corrmat.txt');
%MDS + kmeans
big = 999;
% small = 1;
% th = 0.06;
% data = data(end-14:end,:);
% data = data - th;
data = abs(data);
data = 1./data;
data(logical(eye(size(data))))=0;
data(data==Inf)=big;
% for i=1:numel(data)
%     if data(i)~=0 && data(i)~=big
%         if data(i)>(1/0.06)
%             data(i)=big;
%         else
%             data(i)=small;
%         end
%     end
% end

Y = mdscale(data,3);
% Y(4:9,:)=[];
% Y(7:9,:)=[];
% Y(10:12,:)=[];
[a, C] = kmeans(Y,5)
b = num2str(a);
c = cellstr(b);

figure
hold on
x = Y(:,1)';
y = Y(:,2)';
z = Y(:,3)';
color = 'rgbmk';
for i=1:5
    ox = i*3-3+1;
    plot3(x(ox:ox+2),y(ox:ox+2),z(ox:ox+2), 'color', color(i),'Line','o');
    text(x(ox:ox+2),y(ox:ox+2),z(ox:ox+2), c(i*3-3+1:i*3-3+1+2),'VerticalAlignment','bottom',...
                                                                   'HorizontalAlignment','right');
end

%thresholding-based algo
% data = data(9:end-4, 9:end-4);
% room = length(data);
% NN = 3;
% clx_res = zeros(room,NN);
% count = 0;
% %brute force, assign a point to the cluster it has MAX corrcoef with
% for i = 1:room
%     cor = data(:,i);
%     cor(cor==1) = 0;
%     cor = abs(cor);
%     for k =1:NN
%         [v, clx_res(i,k)] = max(cor);
%         cor(clx_res(i,k)) = 0;
%         clx_res(i,k) = ceil(clx_res(i,k)/4);
%     end
% %     if clx_res(i,1) == ceil(i/4)
% %         count = count+1;
% %     end
%     if sum(clx_res(i,:)==ceil(i/4))~=0
%         count = count+1;
%     end
% end
% accuracy1 = count / room

% %reconstruction acc (upper bound)
% count = 0;
% for i=1:room/4
%     tmp = data(4*i-3:4*i,4*i-3:4*i);
%     tmp(tmp==1) = 0;
%     for j=1:4
%         if max(tmp(:,j))>0.07
%             count = count+1;
%         end
%     end
% end
% recon_acc = count/room

%kNN
% for i = 1:room
%     cor = data(:,i);
%     cor(cor==1) = 0;
% 
%     [v, clx_res(i)] = max(cor);
%     clx_res(i) = floor(clx_res(i)/4);
%     if clx_res(i) == floor(i/4)
%         count = count+1;
%     end
% end

%clear type pattern across rooms and then 
% turn cor_mat into graph connection mat
% data = abs(data);
% for n=1:4
%     data(n:4:end, n:4:end) = 0.01;
% end
% for n=1:length(data)
%     data(n,n) =0;
% end
% con = data;
% con(con>0.07) = 1;
% con(con<=0.07) = 0;
% 
% %ones counting
% indice = nchoosek(1:length(con),4);
% % res = [];
% FN = 'matrix_analysis.txt';
% FID = fopen(FN,'w');
% counter1 = 0;
% counter2 = 0;
% for i=1:length(index)
%     seq = indice(i,:);
%     tmp = con(seq,seq);
%     degree = sum(sum(tmp))/2;
%     if degree>=3
% %         res = [res; seq degree];
%         counter2 = counter2+1;
%         fprintf(FID,'%d %d %d %d %d\n',seq,degree);
%     end
%     counter1 = counter1+1;
%     if mod(counter1,10000)==0
%         fprintf('>>>>>finished %d comb, got %d statisfied<<<<<\n', counter1, counter2)
%     end
% end
% toc