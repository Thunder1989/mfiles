%% plot HeatMap of each type
%and the distribution btw 
% clear
% clc
% data = importdata('/Users/hdz_1989/Documents/Dropbox/SDB/result/KETI/keti_xcormat_nooutlierremovalonlight_MA.txt');
figure
for i=1:4
    cor = data(i:4:end,i:4:end);
    cor = abs(cor);
%     HeatMap(cor);
    cor(cor==1) = [];%remove X-X pair in room on the diagonal
    [p,v] = ecdf(cor);
    subplot(4,1,i)
    plot(v,p,'k--','LineWidth',1.5)
    grid on
end

%% cdf for CO2-Hum pair
CH_cor = data(1:4:end,2:4:end);
CH_intra = diag(CH_cor);%get the pair in the same room
CH_intra = abs(CH_intra);
for i=1:length(CH_cor)
    CH_cor(i,i) = 1;
end
CH_cor(CH_cor==1) = [];
CH_inter = abs(CH_cor);
[p1,v1] = ecdf(CH_intra);
[p2,v2] = ecdf(CH_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'CH-intra','CH-inter','Location','SouthEast')

%% cdf for CO2-Temp pair
CT_cor = data(1:4:end,4:4:end);
CT_intra = diag(CT_cor);%get the pair in the same room
CT_intra = abs(CT_intra);
for i=1:length(CT_cor)
    CT_cor(i,i) = 1;
end
CT_cor(CT_cor==1) = [];
CT_inter = abs(CT_cor);
[p1,v1] = ecdf(CT_intra);
[p2,v2] = ecdf(CT_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'CT-intra','CT-inter','Location','SouthEast')

%% cdf for CO2-Lumin pair
CL_cor = data(1:4:end,3:4:end);
CL_intra = diag(CL_cor);%get the pair in the same room
CL_intra = abs(CL_intra);
for i=1:length(CL_cor)
    CL_cor(i,i) = 1;
end
CL_cor(CL_cor==1) = [];
CL_inter = abs(CL_cor);
[p1,v1] = ecdf(CL_intra);
[p2,v2] = ecdf(CL_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'CL-intra','CL-inter','Location','SouthEast')

%% cdf for Hum-Temp pair
HT_cor = data(2:4:end,4:4:end);
HT_intra = diag(HT_cor);%get the pair in the same room
HT_intra = abs(HT_intra);
for i=1:length(HT_cor)
    HT_cor(i,i) = 1;
end
HT_cor(HT_cor==1) = [];
HT_inter = abs(HT_cor);
[p1,v1] = ecdf(HT_intra);
[p2,v2] = ecdf(HT_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'HT-intra','HT-inter','Location','SouthEast')

%% cdf for Hum-Lumin pair
HL_cor = data(2:4:end,3:4:end);
HL_intra = diag(HL_cor);%get the pair in the same room
HL_intra = abs(HL_intra);
for i=1:length(HL_cor)
    HL_cor(i,i) = 1;
end
HL_cor(HL_cor==1) = [];
HL_inter = abs(HL_cor);
[p1,v1] = ecdf(HL_intra);
[p2,v2] = ecdf(HL_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'HL-intra','HL-inter','Location','SouthEast')

%% cdf for Temp-Lumin pair
TL_cor = data(4:4:end,3:4:end);
TL_intra = diag(TL_cor);%get the pair in the same room
TL_intra = abs(TL_intra);
for i=1:length(TL_cor)
    TL_cor(i,i) = 1;
end
TL_cor(TL_cor==1) = [];
TL_inter = abs(TL_cor);
[p1,v1] = ecdf(TL_intra);
[p2,v2] = ecdf(TL_inter);
figure
hold on
f1 = plot(v1,p1,'k:','LineWidth',1.5);
f2 = plot(v2,p2,'r--','LineWidth',1.5);
grid on
legend([f1 f2],'TL-intra','TL-inter','Location','SouthEast')


%% check error upper-bound for C-L pair
CL_cor = data(1:4:end,3:4:end);
count = 0;
for i=1:length(CL_cor)
  if CL_cor(i,i)==max(CL_cor(i,:))
      count = count + 1;
  end
end
count/length(CL_cor)
