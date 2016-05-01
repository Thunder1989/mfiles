% The script computes the EMDs of a signal and return a matrix,
% each IMF is a colomn vector

clear
clc
data = csvread('/Users/hdz_1989/Downloads/SDB/d/temp_d.csv');
% data = csvread('/Users/hdz_1989/Downloads/SDB/temp.csv');
% data = wavread('/Users/hdz_1989/Documents/Dropbox/sounds/STG_LungS_Asthma_Wz.wav');
% data(data==0) = [];
% figure
% plot(data(:,1), data(:,2))
Nstd = 0.2;
NR = 500;
MaxIter = 5000;
% result = ceemdan(data(:,2), Nstd, NR, MaxIter);
% data = data(:,2);
% data_length = 120000;
% len = length(data);
% [p,q] = rat(data_length/len);
% data = resample(data, p, q);
% data = data(1:data_length);

% data = EWMA(data',10)';
result = emd(data); %emd(x)
% result = emd(data(:,2), 'STOP', [0.05,0.5,0.05]);
%it will take a while to run the script...
figure
len = size(result,1);
row = size(result,2);
for i = 1:len
    subplot(len+1,1,i)
    plot(result(i,:))
end
subplot(len+1,1,len+1)
plot(data(:,1))
%fprintf('got %d IMFs\n', len)
% temp = csvread('/Users/hdz_1989/Downloads/SDB/temp.csv');
% figure
% plot(temp(:,2))