close all
clear
clc

T = 28; % # of days
Niter = 50;
Nburn = 10;
Nplot = 20;

path_ahu = './data_ahu/';
path_vav = './data_vav/';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

%%
%load ahu data
num = length(ahus);
ahu_dist = cell(num,1);
ahu_list = zeros(num,1);
% figure
prob_ahu = cell(num,1);
base_ahu = cell(num,1);
N0_ahu = cell(num,1);
NE_ahu = cell(num,1);
for n = 4:4
    fn = [path_ahu, ahus(n).name];
    cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
    cur_ahu = cur_ahu(1:4*24*T,:);
    cur_ahu = cur_ahu(:,1);
    cur_ahu = reshape(cur_ahu, 96, []);
    ahu_list(n) = str2double(ahus(n).name(5));

    % AHU col 3 - SupplyAirPress,  5 - SupplyFanSpeedOutput (the # of cols varies..
%         cur_corr = corrcoef(cur_ahu(:,5), cur_vav(:,1));
%         cur_corr = cur_corr(1,2);
%         vav_corr(n) = abs(cur_corr);
    event_times = zeros(size(cur_ahu));
    res = mmpp(cur_ahu, priors, [Niter,Nburn,Nplot], event_times, [3,3]);
 
    data = res.Z(:,:,end);
    prob_ahu{n} = reshape(data,[],1)';
    data = res.L(:,:,end);
    base_ahu{n} = reshape(data,[],1)';
    data = res.N0(:,:,end);
    N0_ahu{n} = reshape(data,[],1)';
    data = res.NE(:,:,end);
    NE_ahu{n} = reshape(data,[],1)';
%     subplot(length(ahus),1,n)
%     plotyy(1:numel(data), reshape(cur_ahu,[],1), 1:numel(data), reshape(data,[],1), @plot, @stem);
%     stem( reshape(mean(res.Z, 3),[],1), 'k', 'MarkerSize', 1);
end

%%
%load vav data
% figure

num = length(vavs);
prob_vav = cell(num,1);
base_vav = cell(num,1);
N0_vav = cell(num,1);
NE_vav = cell(num,1);

for m = 1:9
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
    data = res.L(:,:,end);
    base_vav{m} = reshape(data,[],1)';
    data = res.N0(:,:,end);
    N0_vav{m} = reshape(data,[],1)';
    data = res.NE(:,:,end);
    NE_vav{m} = reshape(data,[],1)';

    figure
    hold on
    plot(reshape(cur_vav, 1, []),'k--')
    plot(N0_vav{m},'r--')
    plot(NE_vav{m},'g--')
    
%     res.D(:,:,end)
%     vav_corr = zeros(length(ahus),1);
end

%% plot
figure
num = 8;
for m=1:num
    subplot(num,2,2*m-1)
    plot(N0_ahu{m},'b')
    subplot(num,2,2*m)
    plot(NE_ahu{m},'k')
end

%% correlating on N0
ahu_list = [1 2 4 5 6 7 8 9];
vav_corr = zeros(length(vavs), length(ahus));
ctr = 0;
for m=1:length(vavs)
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    cur_vav = N0_vav{m};

    for n=1:length(base_ahu)
        fn = [path_ahu, ahus(n).name];
        cur_ahuid = str2double(ahus(n).name(5));
    
        cur_ahu = N0_ahu{n};
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

ahu_tmp = cell2mat(prob_ahu);
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

%% MCMC test
trial = 1000000;
ctr = 0;
for i=1:trial
    if rand(1)^2+rand(1)^2<=1
        ctr = ctr+1;
    end
end
est = ctr/trial*4
