% The script computes the EMDs of a signal and return a matrix,
% each IMF is a colomn vector


clear
clc
data = csvread('/Users/hdz_1989/Downloads/SDB/energy.csv');
%data = csvread('/Users/hdz_1989/Downloads/SDB/temp.csv');
% data = wavread('/Users/hdz_1989/Documents/Dropbox/sounds/STG_LungS_Norm_Tracheal.wav');
% data(data==0) = [];
% figure
% plot(data(:,1), data(:,2))
result = emd(data(:,2), 40, 40, 1);
%it will take a while to run the script...
figure
len = size(result,2);
row = size(result,1);
for i = 1:len
    subplot(len+1,1,i)
    plot(result(:,i))
end
subplot(len+1,1,len+1)
plot(data(:,2))
%fprintf('got %d IMFs\n', len)
temp = csvread('/Users/hdz_1989/Downloads/SDB/temp.csv');
figure
plot(temp(:,2))