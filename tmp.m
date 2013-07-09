%to compare the baseline from IMF
a = 7;
b = 8;
%data = csvread('/Users/hdz_1989/Downloads/energy.csv');
%result = csvread('/Users/hdz_1989/Downloads/imfs.csv');
figure
hold on
plot(data(:,2),'r--')
res = zeros(size(result,1),1);
for k = a:b
    res = res + result(:,k);
end
plot(res)
%plot((data(:,2)-res), 'b--')