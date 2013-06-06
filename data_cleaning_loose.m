clear
clc
[data, txt] = xlsread('/Users/hdz_1989/Documents/MATLAB/GPS_Data/User_08','Home');
%[data, txt] = xlsread(str,'Home');
length = size(txt,1);

tmp = txt(:,3);
tmp = char(tmp);
t_lh = zeros(length,2);
t_lh(:,1) = str2num(tmp(:,5:8));
t_lh(:,2) = str2num(tmp(:,9:12))/100;
tmp = txt(:,2);
tmp = char(tmp);
t_ah = zeros(length,2);
t_ah(:,1) = str2num(tmp(:,5:8));
t_ah(:,2) = str2num(tmp(:,9:12))/100;

%find all the dates
date = t_ah(1,1);
count = 1;
for i = 1:length
    pos = find(date == t_ah(i,1));
    if size(pos,2)==0
        count = count + 1;
        date(count) = t_ah(i,1);
    end
    pos = find(date == t_lh(i,1));
    if size(pos,2)==0
        count = count + 1;
        date(count) = t_lh(i,1);
    end
%     if t_ah(i) ~= date(count)
%         count = count+1;
%         date(count) = t_ah(i);
%     end
end
date = date';

T_LH = zeros(size(date,1),1);
T_AH = zeros(size(date,1),1);
%T_LH(1) = t_lh(1,2);
%T_AH(1) = t_ah(1,2);

%get T_AH
for i = 1:length
    count = find(date == t_ah(i,1));
    if size(count,1)~=0
        if T_AH(count)==0 && t_ah(i,2)>=11.00 %this value might be changed
            T_AH(count) = t_ah(i,2);
        elseif T_AH(count)~=0 && t_ah(i,2)>T_AH(count) && t_ah(i,2)<=22.00
            T_AH(count) = t_ah(i,2);
        elseif t_ah(i,2)<=3.00
            T_AH(count-1) = t_ah(i,2)+24;
        end
    end
%     if t_ah(i,1)==date(count) && t_ah(i,2)>=11.00
%         T_AH(count) = t_ah(i,2);
%     elseif t_ah(i,1)~=date(count)
%         count = count+1;
%         T_AH(count) = t_ah(i,2);
%     end
%    if T_AH(count)<=8 && T_AH(count)~=0
%         T_AH(count) = T_AH(count) + 24;
%    elseif T_AH>8 && T_AH<12
%       T_AH(count) = 0;
%    end
end
%get T_LH
for i = 1:length
    count = find(date == t_lh(i,1));
    if size(count,1)~=0
        if T_LH(count)==0 && t_lh(i,2)>=5.00
            T_LH(count) = t_lh(i,2);
        elseif t_lh(i,2)<=18.00 && t_lh(i,2)>=6.00
            T_LH(count) = t_lh(i,2);
        end
    end
    
    if T_LH(count)> T_AH(count) && T_AH(count)~=0 %if LH> AH
        T_LH(count) = 0;
    end
end

%get time for Work
[data, txt] = xlsread('/Users/hdz_1989/Documents/MATLAB/GPS_Data/User_08','Work');
%[data, txt] = xlsread(str,'Work');
length = size(txt,1);
tmp = txt(:,2);
tmp = char(tmp);
t_aw = zeros(length,2);
t_aw(:,1) = str2num(tmp(:,5:8));
t_aw(:,2) = str2num(tmp(:,9:12))/100;
tmp = txt(:,3);
tmp = char(tmp);
t_lw = zeros(length,2);
t_lw(:,1) = str2num(tmp(:,5:8));
t_lw(:,2) = str2num(tmp(:,9:12))/100;

T_AW = zeros(size(date,1),1);
T_LW = zeros(size(date,1),1);
%get T_AW
for i = 1:length
    count = find(date == t_aw(i,1));
    if size(count,1)~=0
        if T_AW(count) == 0
            T_AW(count) = t_aw(i,2);
        end
    end
end
%get T_LW
for i = 1:length
    count = find(date == t_lw(i,1));
    if size(count,1)~=0
        T_LW(count) = t_lw(i,2);
    end
end

%get D_H & D_W
length = size(date,1);
D_H = zeros(length,1);
D_W = zeros(length,1);
%D_H, the duration staying home last night
for i = 1:length
    if i == 1
        D_H(i) = toMinute(T_LH(i));
    else
        if T_AH(i-1) >= 24.00
            D_H(i) = toMinute(T_LH(i))-toMinute(T_AH(i-1)-24);
        else
            D_H(i) = toMinute(T_LH(i)+24)-toMinute(T_AH(i-1));
        end
    end
end
%D_W, currently assume no T_LW after 12am, needs check later
for i = 1:length
    D_W(i) = toMinute(T_LW(i))-toMinute(T_AW(i));
end

%PREHEAT prediction algorithm implementation
length = size(T_AH,1);
interval = 5;
t = max(T_AH);
seq_num = (t-mod(t,1))*(60/interval)+floor(mod(t,1)*100/interval)+1;%extend a day to include case when T_AH on next day
day_seq = zeros(length, seq_num);
%generate day occupancy vectors, currently not extend to the next day
for i = 1:length
    t = T_LH(i);
    slot1 = (t-mod(t,1))*(60/interval)+floor(mod(t,1)*100/interval)+1;
    t = T_AH(i);
    slot2 = (t-mod(t,1))*(60/interval)+floor(mod(t,1)*100/interval)+1;
    
    if slot1~=0 && slot2~=0 %neither T_LH nor T_AH is zero
        for j = 1:seq_num
            if j<slot1
                day_seq(i,j) = 1;
            elseif j>slot2
                day_seq(i,j) = 1;
            else
                day_seq(i,j) = 0;
            end
        end
    end
end
%prediction
count = 0;%counter of correct predictions
base = 0;%counter of attempted predictions
t_ahead = 60;
slot_num = t_ahead/interval;
dist_mat = zeros(length,length);
N = 5;%# of kNN
pdt_err = 0;
ctr = 1;
for i = 1:length
    flag = true;
    t = T_LW(i);
    slot1 = (t-mod(t,1))*(60/interval)+floor(mod(t,1)*100/interval)+1;
    t = T_AH(i);
    slot2 = (t-mod(t,1))*(60/interval)+floor(mod(t,1)*100/interval)+1;

    if slot2~=0
        for j = slot1:seq_num-slot_num
            %predict after the person leaves home
            for k = 1:length
                if k~=i
                    %compute M-Distance btw current and historical days
                    dist_mat(i,k) = mDist(day_seq(i,1:j), day_seq(k,1:j));
                end
            end
            dist_mat(i,i) = seq_num;

            index = kNN(dist_mat(i,:),N);%get the index of kNN
            pdt_mat = zeros(1,slot_num);
            %sum up and average
            for k = 1:N
                pdt_mat = pdt_mat + day_seq(index(k),j+1:j+slot_num);
            end
            pdt_mat = round(pdt_mat/N);

            for k = 1:slot_num
                if pdt_mat(k) == day_seq(i,j+k)
                    count = count + 1;
                end
                
                if pdt_mat(k)==1 && flag
                    pdt_err(ctr) = j+k-slot2;
                    flag = false;
                    ctr = ctr + 1;
                end
            end
            base = base + slot_num;
        end
    end
end
accuracy1 = count/base

%LEAVE-ONE-OUT CROSS-VALIDATION
count = 0;
length = size(T_AH,1);
%arriving time 11:00am-2:15am
interval = 5;
t = max(T_AH);
tmin = min(T_AH);
tmin = (tmin-mod(tmin,1))*(60/interval)+ceil((mod(tmin,1)*100-1)/interval)+1;
len = (t-mod(t,1))*(60/interval)+ceil((mod(t,1)*100-1)/interval)+1-tmin+1;%len is the # of arriving time slots
step = 6;
freq_count = zeros(step+1,len);
distribution = zeros(step+1,len*2-1);

for k = 1:length
    t = T_AH(k);
    if t~=0
        tmp = (t-mod(t,1))*(60/interval)+ceil((mod(t,1)*100-1)/interval)+1-tmin+1;
        freq_count(step+1, tmp) = freq_count(step+1, tmp) + 1;
    end
end

for i = 1:length
    sampleT_LH = T_LH(i);
    sampleT_LW = T_LW(i);
    sampleD_W = D_W(i);
    sampleD_H = D_H(i);
    sampleT_AH = T_AH(i);
    index = 0;
    ctr = 1;
    %AH error distribution
    for j = 1:step
        %get count in each slot
        for k = 1:length
           if abs(toMinute(sampleT_LW)-toMinute(T_LW(k)))<=5*j && abs(toMinute(sampleT_LH)-toMinute(T_LH(k)))<=5*j %&& abs(sampleD_H-D_H(k))<=10*j
               %count0 = count0 + 1;
               t = T_AH(k);
               if t~=0 && k~=i
                   tmp = (t-mod(t,1))*(60/interval)+ceil((mod(t,1)*100-1)/interval)+1-tmin+1;
                   freq_count(j, tmp) = freq_count(j, tmp)+1;
               end
               
%                index(ctr) = k;
%                ctr = ctr+1;
           end
        end
        
        t = sampleT_AH;
        if t~=0
            tmp = (t-mod(t,1))*(60/interval)+ceil((mod(t,1)*100-1)/interval)+1-tmin+1;
            time_diff = ((1:len)' - tmp);
            for k = 1:len
                distribution(j, time_diff(k)+len) = distribution(j, time_diff(k)+len) + freq_count(j, k);
            end
        end
        freq_count(j,:) = 0;
    end
    
    %simple overall AH distribution, baseline
    t = sampleT_AH;
    if t~=0
        tmp = (t-mod(t,1))*(60/interval)+ceil((mod(t,1)*100-1)/interval)+1-tmin+1;
        time_diff = ((1:len)' - tmp);
        for k = 1:len
            distribution(step+1, time_diff(k)+len) = distribution(step+1, time_diff(k)+len) + freq_count(step+1, k);
        end
    end
end

%Turn the errors in ETA into energy savings/penalties.
%negative erros bring maintain cost, which is 1 kWh/h
%positive erros bring low-efficiency cost(need to heat quickly),
%which is 1.5 kWh/h
energy_penalty = zeros(step+2,2);
for i = 1:step+1
    distribution(i,:) = distribution(i,:)/max(cumsum(distribution(i,:)));
end

react_length = 1;
maintain_penalty = 1.5;
react_penalty = 4;
for j = 1:step+1
    for i = 1:size(distribution(j,:),2)
        if i<=len
            energy_penalty(j,1) = energy_penalty(j,1) + maintain_penalty*((len-i)*interval/60)*distribution(j,i);
        else
            energy_penalty(j,2) = energy_penalty(j,2) + react_penalty*distribution(j,i);
            %energy_penalty(j,2) = energy_penalty(j,2) + (i-len)*interval*distribution(j,i)/error_max(j);
        end
    end
end
[pr, x] = ecdf(pdt_err);
for i = 1:size(x,1)
    if i==1
        energy_penalty(step+2,1) = energy_penalty(step+2,1) + maintain_penalty*abs(x(i))*interval/60*pr(i);
    else
        if x(i)<0
            energy_penalty(step+2,1) = energy_penalty(step+2,1) + maintain_penalty*abs(x(i))*interval/60*(pr(i)-pr(i-1));
        else
            energy_penalty(step+2,2) = energy_penalty(step+2,2) + react_penalty*(pr(i)-pr(i-1));
        end
    end
end

% s = energy_penalty(:,1) + energy_penalty(:,2);
% mean(energy_penalty(1:6,2))
% std(energy_penalty(1:6,2))
% energy_penalty(7,2)

%time error cdf
% xx = (-(len-1):(len-1))*interval;
% figure
% hold on
% set(gca,'YGrid','on')
% plot(xx, cumsum(distribution(1,:))/max(cumsum(distribution(1,:))),':','Color',[210 105 30]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(2,:))/max(cumsum(distribution(2,:))),'--','Color',[210 105 30]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(3,:))/max(cumsum(distribution(3,:))),'Color',[210 105 30]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(4,:))/max(cumsum(distribution(4,:))),':','Color',[0 100 0]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(5,:))/max(cumsum(distribution(5,:))),'--','Color',[0 100 0]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(6,:))/max(cumsum(distribution(6,:))),'Color',[0 100 0]/255, 'LineWidth', 1.4)
% plot(xx, cumsum(distribution(7,:))/max(cumsum(distribution(7,:))),'Color',[178 34 34]/255, 'LineWidth', 1.5)
% lh = legend('5-min','10-min','15-min','20-min','25-min','30-min','baseline');
% set(lh,'Location','SouthEast','Orientation','vertical')
% set(get(gca,'XLabel'),'String','Prediction Error (Minute)')
% set(get(gca,'YLabel'),'String','CDF of prediction error (%)')
% axis([-600 600 0 1])

%time error pdf
% figure
% hold on
% plot(distribution(1,:)/error_max(1),'r')
% plot(distribution(3,:)/error_max(2),'g')
% plot(distribution(5,:)/error_max(2),'b')
% plot(distribution(13,:)/error_max(13),'k:')
% plot([len,len],[0,max(distribution(13,:))/error_max(13)])
% lh = legend('case1','case2','case3','baseline');
% set(lh,'Location','NorthEast','Orientation','vertical')

%time distribution bar of LH, AH, LW
% figure
% colormap(gray);
% subplot(3,1,1)
% bar(T_LH)
% title('Leave Home')
% subplot(3,1,2)
% bar(T_AH)
% title('Arrive Home')
% % subplot(3,1,3)
% % bar(data(:,4))
% % title('Arrive Work')
% subplot(3,1,3)
% bar(T_LW)
% % datetick('x','d')
% title('Leave Work')

%time distribution with ksdensity
% figure
% set(gca,'YGrid','on')
% hold on
% T_LH(T_LH == 0) = [];%remove 0s
% T_LW(T_LW == 0) = [];
% [f, xi] = ksdensity(T_LH, 'width', 0.75);
% plot(xi, f, ':','Color',[140, 30, 29]/255, 'LineWidth',3);
% [f, xi] = ksdensity(T_LW, 'width', 0.75);
% plot(xi, f, '--','Color',[47, 107, 81]/255,'LineWidth',3);
% lh = legend('Time of Leaving Home','Time of Leaving Work');
% set(lh,'Location','NorthWest','Orientation','vertical')
% axis([0 30 0 max(f)+0.1])

%energy bar
figure
% subplot(1,2,1)
h = bar(energy_penalty(1:step+2,:),'stack');
set(h,'BarWidth',0.6); % The bars will now touch each other
set(gca,'YGrid','on')
%set(gca,'GridLineStyle','+')
%set(gca,'XTicklabel','LW>=sD&&LH>=sD|LW>=sD&&LH<=sD|LW<=sD&&LH>=sD|LW<=sD&&LH<=sD')
set(gca,'XTicklabel','5|10|15|20|25|30|baseline|preheat')
set(get(gca,'YLabel'),'String','Energy Consumption/kWh')
lh = legend('Maintain','React');
set(lh,'Location','SouthEast','Orientation','vertical')
colormap(summer)
% subplot(1,2,2)
% boxplot([T_LH,T_LW])
