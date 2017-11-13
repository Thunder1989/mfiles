close all
clear
clc

T = 7*2; % # of days
D = 7*2; % # of days to plot
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

%% gibbs sampling based
% close all
clc
num = length(ahus);
ahu_dist = cell(num,1);
ahu_list = zeros(num,1);
event_ahu = cell(num,1);
N0_ahu = cell(num,1);
NE_ahu = cell(num,1);
% figure
for n = 2:2
    fn = [path_vav, vavs(n).name];
    cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
    cur_ahu = cur_ahu(1:4*24*T,:);
    cur_ahu = cur_ahu(:,1); %ahu last col, vav 1st col
    tic
    [res, Z, M] = gibbs_sgf_K(cur_ahu, 2, 1, 0);
%     [res, Z, M] = gibbs_sgf_K_para(cur_ahu, 2, 1, 0);
    toc

%     figure
%     hold on
% %     yyaxis left
%     %events in shade
%     stem(1:length(Z), Z*max(cur_ahu), 'Marker','None', 'LineWidth',4, 'Color',[.8 .8 .8]);
%     %original data
%     plot(1:length(res), cur_ahu,'k-')
%     %filtered data
%     plot(1:length(res), res(:,1),'g-')
% %     yyaxis right
%     %velocity
%     plot(1:length(res), res(:,2),'r')
% %     area(1:length(Z), Z*max(cur_ahu), 'EdgeColor', 'none', 'FaceColor', [.8 .8 .8]);
%     legend('event','original','filtered','vel')

end

Z_ = Z;

%% self-defined colors
Z = mean(Z_(:,11:2:end),2);
figure
hold on
colors = containers.Map(1:4,{[54/255,160/255,204/255],[211/255,142/255,194/255],[80/255,180/255,110/255],[.8 .8 .455]});
y_lim = max(cur_ahu);
for i=1:length(Z)
    stem(i, y_lim, 'Marker','None', 'LineWidth', 4, 'Color', colors(ceil(Z(i))) );
end
plot(cur_ahu,'k')
ylim([min(cur_ahu)-5 max(cur_ahu)+5])

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

%%
% Z = mode(Z_(:,end-10:end),2);
Z = mode(Z_(:,11:3:end),2);
figure
hold on
%     yyaxis left
%events in shade
stem(1:length(Z), Z*max(cur_ahu), 'Marker','None', 'LineWidth', 4, 'Color', [.8 .8 .8]);
% stem(1:length(Z), (1-Z)*max(cur_ahu), 'Marker','None', 'LineWidth',4, 'Color',[.8 .8 .8]);
% stem(-10*event_times(:), 'LineWidth',1, 'Color',[.3 .3 .3]);
%original data
plot(1:length(res), cur_ahu,'k-')
%filtered data
plot(1:length(res), res(:,1),'g-')
%     yyaxis right
%velocity
plot(1:length(res), res(:,2),'r')
%     area(1:length(Z), Z*max(cur_ahu), 'EdgeColor', 'none', 'FaceColor', [.8 .8 .8]);
legend('event','original','filtered','vel')

%% clustering
k = 6;
data = cellfun(@transpose, vav, 'UniformOutput', false);
% [idx,C] = kmeans(cell2mat(data),8);
options = [];
options.NeighborMode = 'KNN';
options.k = 10;
options.WeightMode = 'Cosine';
% options.t = 1;
W = constructW(cell2mat(data),options);
W_ = max(W, W'); %symmetrize w_, N by N
imagesc(W_);
pause
D = diag(sum(W_,2));
L = D - W_; %unormalized Laplacian
[evc, evl] = eigs(L); %each column of evc is an eigenvector
idx = find(diag(evl)>=0);
input = evc(:,idx(1:k));
%                 input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
c_idx = kmeans(input,k);

%% acc
bid = 320;
week = 4;
path_ahu = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\ahu_common\');
path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav_common\');
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

N = 5*2; %Kalman Filter lookback window size
M = 8; %EWMA window size
T = 7*week; % # of days

ahu_event = cell(length(ahus),1);
ahu_kf_res = cell(length(ahus),1);
ahu_list = zeros(length(ahus),1);
num = length(ahus);
for n = 1:num
%     fprintf('processing %s\n',ahus(n).name)
    fn = [path_ahu, ahus(n).name];
    fn = regexprep(fn,'_cut','_all');
    str = regexp(ahus(n).name,'[0-9]+','match');
    cur_ahuid = str2double(str(1));
    ahu_list(n) = cur_ahuid;
end

num = length(vavs);
vav_edge = cell(num,1);
vav_event = cell(num,1);
vav_kf_res = cell(num,1);
correct = [];
wrong = [];
wrong_test = [];

ctr = 0;
for m = 1:num
%     fprintf('processed %s\n',vavs(m).name);

    fn = [path_vav, vavs(m).name];
    fn = regexprep(fn,'_cut','_all');
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahuid = str2double(str(1));
    e_vav = vav{m};
    
    vav_sim = zeros(length(ahus),1);
    for n = 1:length(ahus)
        fn = [path_ahu, ahus(n).name];
        fn = regexprep(fn,'_cut','_all');
        e_ahu = ahu{n};

        if sum(e_ahu)==0 || sum(e_vav)==0
            vav_sim(n) = 0;
        else
            vav_sim(n) = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav));
        end
    end
    
    true = find(ahu_list==ahuid);
    if ismember( true, find(vav_sim==max(vav_sim)) ) && length( find(vav_sim==max(vav_sim)) )<length(ahus);
        ctr = ctr + 1;
    end
end

fprintf('acc on GloME event corr is %.4f\n', ctr/num);