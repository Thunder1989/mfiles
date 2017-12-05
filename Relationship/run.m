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

%% ne_n0 model
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

% self-defined colors
Z = mode(Z_(:,11:2:end),2);
figure
hold on
colors = containers.Map(1:4,{[54/255,160/255,204/255],[211/255,142/255,194/255],[80/255,180/255,110/255],[.8 .8 .455]});
y_lim = max(cur_ahu);
for i=1:length(Z)
    stem(i, y_lim, 'Marker','None', 'LineWidth', 4, 'Color', colors(ceil(Z(i))) );
end
plot(cur_ahu,'k')
ylim([min(cur_ahu)-5 max(cur_ahu)+5])

%% grouped sgf
% close all
clc
load('320_output.mat');
%process covar
Q = cellfun(@cell2mat, Q_vav, 'UniformOutput', false);
Q = cellfun(@(x) reshape(x,1,[]), Q, 'UniformOutput', false);
Q = cellfun(@sort, Q, 'UniformOutput', false);
Q = cell2mat(Q);

K = 10;
c_id = kmeans(Q, K);
% map = [c_id'; 1:length(Q_vav)]';
% map = sortrows(map,1)';
colors = containers.Map(1:4,{[54/255,160/255,204/255],[211/255,142/255,194/255],[80/255,180/255,110/255],[.8 .8 .455]});
%%
for k = 2:K
    cur_data = [];
    ctr  = 0;
    for n = 1:length(Q_vav)
        if c_id(n) ~= k
            continue
        end
        fn = [path_vav, vavs(n).name];
        cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
        cur_ahu = cur_ahu(1:4*24*T,:);
        cur_ahu = cur_ahu(:,1); %ahu last col, vav 1st col
        cur_data = [cur_data; cur_ahu];
        ctr = ctr + 1;
    end
    fprintf('finish concatenating %d for cluster %d...\n', ctr, k);
    [res, Z, M] = gibbs_sgf_K(cur_data, 2, 1, 0);

    Z_ = Z;
    Z = mode(Z_(:,11:2:end),2);
    figure
    hold on
    y_lim = max(cur_data);
    for i=1:length(Z)
        stem(i, y_lim, 'Marker','None', 'LineWidth', 4, 'Color', colors(ceil(Z(i))) );
    end
    plot(cur_data,'k')
    ylim([min(cur_data)-5 max(cur_data)+5]) 
end

%%
vav_ep = cell(size(prob_vav));
for m=1:9
    vav_cur = prob_vav{m};
    vav_cur = reshape(vav_cur,96,7,[]);
    vav_cur = sum(vav_cur,3)/size(vav_cur,3);
    vav_cur = permute( reshape(vav_cur',7,[],24), [2 1 3] ); %fucking tricky
    vav_cur = reshape( sum(vav_cur)/size(vav_cur,1),7,[] )';    %even trickier
    vav_ep{m} = vav_cur; %prob per hour for each day
end

%% Movie Test
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

%% tfidf alone acc
load('320_events.mat');
ahu_ = cellfun(@transpose,ahu,'UniformOutput',false);
vav_ = cellfun(@transpose,vav,'UniformOutput',false);

% align vav with corresponding ahu
num = size(vav_,1);
for m = 1:num
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahu_id = str2double(str(1));
    idx = find(ahu_list == ahu_id);
    f1 = vav_{m};

    f2 = ahu_{idx};
    f2 = f2 | [false f2(1:end-1)];
    vav_{m} = double(f1 & f2);
end

ahu_res_ = ceil(2*cell2mat(ahu_));
vav_res_ = ceil(2*cell2mat(vav_));

fea_ahu = tfidf(ahu_res_);
fea_vav = tfidf(vav_res_);

% fea_ahu = cell2mat(ahu_);
% fea_vav = cell2mat(vav_);

% merge ahu and vav together for tf-idf
% fea = tfidf([cell2mat(ahu_kf_res); cell2mat(vav_kf_res)]);
% fea_ahu = fea(1:size(ahu_kf_res,1), :);
% fea_vav = fea(size(ahu_kf_res,1)+1:end, :);

ctr = 0;
num = size(fea_vav,1);
for m = 1:num
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahu_id = str2double(str(1));
    f1 = fea_vav(m,:);

    vav_sim = zeros(length(ahus),1);
    for n = 1:length(ahus)
        f2 = fea_ahu(n,:);
        cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
        vav_sim(n) = abs(cur_sim);
    end
    
    if ismember( ahu_id, ahu_list(vav_sim==max(vav_sim)) ) && length( find(vav_sim==max(vav_sim)) ) < length(vav_sim)
        ctr = ctr + 1;
    else
        %TBD: some debugging
    end
        
end

fprintf('acc on tfidf cossim is %.4f\n', ctr/num);

%% global-local check
clc
load('320_events.mat');
ahu_ = cellfun(@transpose,ahu,'UniformOutput',false);
vav_ = cellfun(@transpose,vav,'UniformOutput',false);
ahu_conf = cell2mat(ahu_);
vav_conf = cell2mat(vav_);
ahu_ = round(cell2mat(ahu_));
vav_ = round(cell2mat(vav_));
subset = [1,7];

%take ahu1 and ahu7 data
num = size(ahu_,1);
ahu_tmp = [];
ahu_list = [];
for m = 1:num
    str = regexp(ahus(m).name,'[0-9]+','match');
    ahu_id = str2double(str(1));
    if ismember(ahu_id, subset)
        ahu_tmp = [ahu_tmp; ahu_(m,:)];
        ahu_list = [ahu_list; ahu_id];
    end
end

num = size(vav_,1);
vav_sub = [];
vav_label = [];
vav_conf_sub = [];
for m = 1:num
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahu_id = str2double(str(1));
    if ismember(ahu_id, subset)
        vav_sub = [vav_sub; vav_(m,:)];
        vav_conf_sub = [vav_conf_sub; vav_conf(m,:)];
        vav_label = [vav_label; ahu_id];
    end
end

%tfidf
% fea_ahu = tfidf(ahu_tmp);
% fea_vav = tfidf(vav_tmp);

fea_ahu = ahu_tmp;
fea_vav = vav_sub;

%acc eval
ctr = 0;
num = size(vav_sub,1);
for m = 1:num
    ahu_id = vav_label(m);
    f1 = fea_vav(m,:);

    vav_sim = zeros(length(ahu_tmp),1);
    for n = 1:size(ahu_tmp,1)
        f2 = fea_ahu(n,:);

        cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
        vav_sim(n) = abs(cur_sim);
    end
    
    if ismember( ahu_id, ahu_list(vav_sim==max(vav_sim)) ) && length( find(vav_sim==max(vav_sim)) ) < length(vav_sim)
        ctr = ctr + 1;
    end
end
fprintf('acc before correction is %.4f\n', ctr/num);

%vertical comparison - kmeans
K = 2; %num of topics
N = 2; %num of states for KF output
c_idx = kmeans(vav_sub', K);
assert(length(c_idx) == size(vav_sub,2));

%horizontal comparison - MLE
num = size(vav_sub,1);
TH = 0.6;
vav_sub_updated = vav_sub;
num_updated = zeros(size(vav_sub_updated,1),1);
assert ( isequal(vav_sub_updated, vav_sub) );
for m = 1:num
    ctr = 0;
    vav_cur = vav_sub(m,:);
    conf_cur = vav_conf_sub(m,:);
    for k = 1:K
        z_tmp = vav_cur(c_idx==k);
        p_tmp = zeros(N,1);
        for n = 1:N
            p_tmp(n) = sum(z_tmp==n-1);
        end
        assert( sum(p_tmp) == length(z_tmp) );
        p_tmp = p_tmp/sum(p_tmp); %p(z|topic=k)
        [~, Z] = max(p_tmp);
        
        %updating based on p(z|topic)
        for i = 1:length(vav_sub_updated(m,:))
            conf_tmp = conf_cur(i);
            if conf_tmp <= 0.6 && conf_tmp >= 1-TH && c_idx(i)==k
                if vav_sub_updated(m,i) ~= Z-1
                    vav_sub_updated(m,i) = Z-1;
                    ctr = ctr + 1;
                end
            end
        end
    end
    num_updated(m) = ctr;
end
assert ( isequal(num_updated, sum(vav_sub_updated~=vav_sub, 2) ) )

%acc eval
fea_vav = vav_sub_updated;
ctr = 0;
num = size(vav_sub,1);
for m = 1:num
    ahu_id = vav_label(m);
    f1 = fea_vav(m,:);

    vav_sim = zeros(length(ahu_tmp),1);
    for n = 1:size(ahu_tmp,1)
        f2 = fea_ahu(n,:);

        cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
        vav_sim(n) = abs(cur_sim);
    end
    
    if ismember( ahu_id, ahu_list(vav_sim==max(vav_sim)) ) && length( find(vav_sim==max(vav_sim)) ) < length(vav_sim)
        ctr = ctr + 1;
    end
end
fprintf('acc after correction is %.4f\n', ctr/num);
