close all
clear
clc

T = 7*4; % # of days
D = 7*4; % # of days to plot
Niter = 20;
Nburn = 0;
Nplot = 5;

path_ahu = './data_ahu/';
path_vav = './data_vav/';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

%%
%load ahu data
close all
num = length(ahus);
ahu_dist = cell(num,1);
ahu_list = zeros(num,1);
event_ahu = cell(num,1);
N0_ahu = cell(num,1);
NE_ahu = cell(num,1);
% figure
for n = 2:2
    fn = [path_ahu, ahus(n).name];
    cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
    cur_ahu = cur_ahu(1:4*24*T,:);
    cur_ahu = cur_ahu(:,end);
    cur_ahu = reshape(cur_ahu, 96, []);
    ahu_list(n) = str2double(ahus(n).name(5));

    % AHU col 3-SupplyAirPress,  5-SupplyFanSpeedOutput (the # of cols varies..
%         cur_corr = corrcoef(cur_ahu(:,5), cur_vav(:,1));
%         cur_corr = cur_corr(1,2);
%         vav_corr(n) = abs(cur_corr);
    event_times = zeros(size(cur_ahu));
    res = mmpp(cur_ahu, [Niter,Nburn,Nplot], event_times, [3,3]);
 
    data = res.P_Z(:,:,end);
    event_ahu{n} = reshape(data(2,:),[],1)';
    data = res.X_B(:,:,end);
    N0_ahu{n} = reshape(data,[],1)';
    NE_ahu{n} = reshape(cur_ahu-data,[],1)';

    figure
    hold on; 
%     plot(reshape(cur_ahu(1:4*24*D),1,[]),'k')
%     plot(reshape(N0_ahu{n}(1:4*24*D),1,[]),'b')
%     plot(reshape(NE_ahu{n}(1:4*24*D),1,[]),'r')
%     plot(reshape(event_ahu{n}(1:4*24*D)*20,1,[]),'g')
    plot(res.Mu0,'k')
    plot(res.Sigma0,'r')
%     day_bd = zeros(1,4*24*D);
%     day_bd(1:4*24:end) = 1*max(cur_ahu(:));
%     days = {'Fri','Sat','Sun','Mon','Tue','Wed','Thu'};
%     ctr = 1;
%     for t=1:numel(day_bd)
%         if day_bd(t)==0
%             continue
%         end
%         plot([t t], [-10 day_bd(t)], 'k--', 'LineWidth', 0.5);
%         text(t+30,-10,days{mod(ctr-1,7)+1})
%         ctr = ctr + 1;
%     end
%     fn = sprintf('ahu_%d.png',n);
%     saveas(gcf,fn);    

end


%%
%load vav data
% figure

num = length(vavs);
num = 9;
prob_vav = cell(num,1);
N0_vav = cell(num,1);
NE_vav = cell(num,1);

for m = 1:num
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    cur_vav = csvread(fn,2);
    cur_vav = cur_vav(1:4*24*T,:);
    % VAV col 1 - AirFlowNormalized, 2 - AirValvePosition
    cur_vav = cur_vav(:,1);
    cur_vav = reshape(cur_vav, 96, []);

    event_times = zeros(size(cur_vav));
    res = mmpp(cur_vav, priors, [Niter,Nburn,Nplot], event_times, [3,3]);
    data = res.Z(:,:,end);
    prob_vav{m} = reshape(data,[],1)';
    data = res.N0(:,:,end);
    N0_vav{m} = reshape(data,[],1)';
    data = res.NE(:,:,end);
    NE_vav{m} = reshape(data,[],1)';

%     figure
%     hold on
%     plot(reshape(cur_vav, 1, []),'k--')
%     plot(N0_vav{m},'r--')
%     plot(NE_vav{m},'g--')
    
%     res.D(:,:,end)
%     vav_corr = zeros(length(ahus),1);
end

%% plot
figure
num = 8;
for m=1:num
    subplot(num,2,2*m-1)
    plot(N0_ahu{m},'k')
    subplot(num,2,2*m)
    plot(NE_ahu{m},'r')
end

%% correlating on N0
ahu_list = [1 2 4 5 6 7 8 9];
vav_corr = zeros(length(vavs), length(ahus));
ctr = 0;
for m=1:length(vavs)
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    cur_vav = NE_vav{m};

    for n=1:length(base_ahu)
        fn = [path_ahu, ahus(n).name];
        cur_ahuid = str2double(ahus(n).name(5));
    
        cur_ahu = NE_ahu{n};
        cur_corr = corrcoef(cur_vav, cur_ahu);
        cur_corr = cur_corr(1,2);
        vav_corr(m,n) = abs(cur_corr);
    end
    cur_corr = vav_corr(m,:);
    if ahu_list(cur_corr==max(cur_corr)) == ahuid
        ctr = ctr + 1;
    end
end

% acc = sum( vav_corr(:,1)==max(vav_corr')' )/size(vav_corr,1)
acc = ctr / size(vav_corr,1)

%% gibbs HMM
% close all
clc
num = length(vavs);
ahu_dist = cell(num,1);
ahu_list = zeros(num,1);
event_ahu = cell(num,1);
N0_ahu = cell(num,1);
NE_ahu = cell(num,1);
% figure
for n = 3:3
    fn = [path_vav, vavs(n).name];
    cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
    cur_ahu = cur_ahu(1:4*24*T,:);
    cur_ahu = cur_ahu(:,1);
    res = gibbs_hmm(cur_ahu,0);
    
    figure
    [ax, h1, h2] = plotyy(1:length(res), cur_ahu, 1:length(res), res(:,2));
    hold(ax(1), 'on')
    hold(ax(2), 'on')
    h3 = plot(1:length(res), res(:,1), 'r--', 'Parent', ax(1));

%     plot(cur_ahu,'k','LineWidth',2)
%     plot(res,'b--','LineWidth',2)

end

%% tmp script

vav_tmp = cell2mat(prob_vav);
vav_e = zeros(size(vav_tmp));
% figure
for m = 1:9
    cur = vav_tmp(m,:);
    th = mean(cur) + 1*std(cur);
%     subplot(9,1,m)
%     hold on
%     grid on
%     stem(cur, 'k', 'MarkerSize', 0);
%     plot(cur>=th, 'ro');
    vav_e(m,:) = cur>=th;
end

ahu_tmp = cell2mat(event_ahu);
ahu_e = zeros(size(ahu_tmp));
% figure
for m = 1:8
    cur = ahu_tmp(m,:);
    th = mean(cur) + 1*std(cur);
    ahu_e(m,:) = cur>=th;
%     subplot(8,1,m)
%     hold on
%     grid on
%     stem(cur, 'b', 'MarkerSize', 0);
%     plot(cur>=th, 'ro');
end

% vav_norm = sqrt( sum(vav_e.^2, 2) );
% ahu_norm = sqrt( sum(ahu_e.^2, 2) );
% coef = vav_norm * ahu_norm';
% sim = vav_e * ahu_e' ./ coef;

vav_corr = zeros(size(vav_e,1), size(ahu_e,1));
for i=1:size(vav_e,1)
    for j=1:size(ahu_e,1)
        corr = corrcoef(vav_e(i,:), ahu_e(j,:));
        vav_corr(i,j) = corr(1,2);
    end
end

acc = sum( vav_corr(:,1)==max(vav_corr')' )/size(vav_corr,1)

%%
vav_ep = cell(size(prob_vav));
for m=1:9
    cur = prob_vav{m};
    cur = reshape(cur,96,7,[]);
    cur = sum(cur,3)/size(cur,3);
    cur = permute( reshape(cur',7,[],24), [2 1 3] ); %fucking tricky
    cur = reshape( sum(cur)/size(cur,1),7,[] )';    %even trickier
    vav_ep{m} = cur; %prob per hour for each day
end

%% Movie Test.
 
% Set up some function. 
% Sine between -2*pi and 2*pi.
x = (10*-pi:0.1:10*pi)'; % Note the transpose.
y = sin(x);
fid = figure;
hold on
% The final plot.
plot(x,y, '*');
 
% Set up the movie.
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 30; % How many frames per second.
open(writerObj); 
 
for i=1:size(y)      
    % We just use pause but pretend you have some really complicated thing here...
    pause(0.1);
    figure(fid); % Makes sure you use your desired frame.
    plot(x(i),y(i),'or');
 
    if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    end
 
end
hold off
close(writerObj); % Saves the movie.
