%zone comparsion
clear
clc
data = importdata('/Users/hdz_1989/Documents/Dropbox/SDB/result/KETI/keti_xcormat_nooutlierremovalonlight_MA.txt');

%zone division
zone1 = 1:8;%4F N
zone2 = 9:16;%4F S
zone3 = 17:20;%5F N
zone4 = 21:24;%5F S
zone5 = 25:30;%6F N
zone6 = 31:39;%6F S
zone7 = 40:46;%7F N
zone8 = 47:51;%7F S


