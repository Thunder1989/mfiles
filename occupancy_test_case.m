%test the occupancy periods
clear
clc
path = '/Users/hdz_1989/Downloads/SDB/KETI/test_unoccupied/';
file = dir(strcat(path,'*.csv'));
num = size(file,1);

figure
data = cell(4,1);
for n=1:num
    filename = [path, file(n).name];
%     input = csvread(filename);
%     input = input(:,2);
    input = importdata(filename);
    input = input.data;
%     outlier removal
    if n~=3
        tmp=input(1);
        for k=2:length(input)
            if abs(abs(input(k))-abs(tmp))>3*min(abs(tmp),abs(input(k))) && tmp~=0
                input(k) = tmp;
                continue;
            end
            tmp = input(k);
        end
        input = EWMA(input', 50);%too much noise in the data, so do a exp. moving average
    end
    data{n} = input;
    subplot(num,1,n)
    grid on
    plot(input,'k')
end