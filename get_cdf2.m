clear
clc
timescale = [1,3,5,7,14,21,28];
path = '/Users/hdz_1989/Downloads/SDB/result/intra/';%path to the data folers
file = dir(strcat(path,'*.txt'));
num = size(file,1);
figure

for n=1:size(file,1)
    filename = [path, file(n).name];
    rd = importdata(filename);
    starter = 1;
    d = rd(starter*3-3+1:end,:);
    d(d==1) = [];
    d = abs(d);
    [p,v] = ecdf(d);
    subplot(num,1,n)
    plot(v,p,'r')
end