clear
clc

T = 7*4; % # of days

path_ahu = 'D:\TraneData\cut\ahu_property_file_10320_cut\ahu\';
path_vav = 'D:\TraneData\cut\ahu_property_file_10320_cut\vav\';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

%% correlation accuracy
close all
clc
wrong = [];
correct = [];
ctr = 0;
ctr1 = 0;
k = 3;
topk = 0;
N = 5*2; %Kalman Filter lookback window size
% num = 21;

ahu_event = cell(length(ahus),1);
ahu_kf_res = cell(length(ahus),1);
ahu_sgf = cell(length(ahus),1);
ahu_sgf_res = cell(length(ahus),1);
ahu_mle_res = cell(length(ahus),1);
ahu_event_ondiff = cell(length(ahus),1);
ahu_list = zeros(length(ahus),1);
num = length(ahus);
for n = 1:num
%     fprintf('processing %s\n',ahus(n).name)
    fn = [path_ahu, ahus(n).name];
    str = regexp(ahus(n).name,'[0-9]+','match');
    cur_ahuid = str2double(str(1));
    data_ahu = csvread(fn,1);
    data_ahu = data_ahu(1:4*24*T,end);

    delta = [0 diff(data_ahu)'];
    e_ahu = edge(repmat(data_ahu',3,1),0.4); %th = 0.4
    e_ahu = e_ahu(1,:);
    e_ahu = e_ahu | [false e_ahu(1:end-1)] | [e_ahu(2:end) false];    
%     e_ahu = delta .* double(e_ahu);
    assert(sum(e_ahu)~=0);

    ahu_event{n} = double(e_ahu);
    ahu_list(n) = cur_ahuid;
    
    kf = StandardKalmanFilter(data_ahu',8,N,'EWMA'); 
%     kf = gibbs_hmm_uni(data_ahu,0);
    diff2 = abs(data_ahu' - kf);
    diff2(isnan(diff2)) = 0;
    ahu_kf_res{n} = diff2(2:end-1); %TBD: make the manual period self-deciding

%     [mle, x] = data_mle(data_ahu',1); 
%     diff2 = abs(data_ahu - mle(:));
%     diff2(isnan(diff2)) = 0;
%     ahu_mle_res{n} = diff2;

%     [y,z,p] = gibbs_sgf(data_ahu,0);
%     ahu_sgf{n} = z;
%     diff2 = abs(data_ahu - y(:,1));
%     diff2(isnan(diff2)) = 0;
%     ahu_sgf_res{n} = diff2;

%     ahu_event_ondiff{n} = get_z_hmm(data_ahu);

end

num = length(vavs);
vav_edge = cell(num,1);
vav_event = cell(num,1);
vav_kf_res = cell(num,1);
vav_sgf = cell(num,1);
vav_sgf_res = cell(num,1);
vav_mle_res = cell(num,1);
vav_event_ondiff = cell(length(ahus),1);
res = zeros(num,length(ahus)+1);
wrong_test = [];
score_tmp = [];
w = 0.5;
debug = 0;
for m = 1:num
    fn = [path_vav, vavs(m).name];
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahuid = str2double(str(1));
    data_vav = csvread(fn,1);
    data_vav = data_vav(1:4*24*T,1);
    
    delta = [0 diff(data_vav)'];
    e_vav = edge(repmat(data_vav',3,1),0.6); %th = 1.25 for 320 596, 0.6 for 642
    e_vav = e_vav(1,:);
    vav_edge{m} = double(e_vav);
    e_vav = e_vav | [false e_vav(1:end-1)] | [e_vav(2:end) false];
%     e_vav = delta .* double(e_vav);
    e_vav = double(e_vav);
    assert(sum(e_vav)~=0);
    vav_event{m} = e_vav;

    vav_corr = zeros(length(ahus),1);
    vav_sim = zeros(length(ahus),1);
    vav_score = zeros(length(ahus),1);
    vav_score1 = zeros(length(ahus),1);
    
    kf = StandardKalmanFilter(data_vav',8,N,'EWMA'); 
%     kf = gibbs_hmm_uni(data_vav,0);
    diff1 = abs(data_vav' - kf);
    diff1(isnan(diff1)) = 0;
    diff1 = diff1(2:end-1);
    vav_kf_res{m} = diff1;

%     [y,z,p] = gibbs_sgf(data_vav,0);
% %     vav_sgf{m} = z;
%     diff1 = abs(data_vav - y(:,1));
%     diff1(isnan(diff1)) = 0;
%     vav_sgf_res{n} = diff1;
%     z_vav = z;

%     [mle, x] = data_mle(data_vav',1); 
%     diff1 = abs(data_vav - mle(:));
%     diff1(isnan(diff2)) = 0;
%     vav_mle_res{n} = diff1;

%     z_vav = get_z_hmm(data_vav);
%     vav_event_ondiff{n} = z_vav;

    for n = 1:length(ahus)
        fn = [path_ahu, ahus(n).name];
        cur_ahuid = ahu_list(n);
        data_ahu = csvread(fn,1);
        data_ahu = data_ahu(1:4*24*T,end);
        e_ahu = ahu_event{n};
        diff2 = ahu_kf_res{n};
%         z_ahu = ahu_event_ondiff{n};
%         z_ahu = ahu_sgf{n};
        
        %debugging block
        if debug == 1
            figure
            hold on
            plot(data_vav,'k','LineWidth',1.5)
            plot(data_ahu,'r','LineWidth',1.5)
%             stem(z_vav*max(data_vav),'b','LineWidth',0.5,'Marker','none')
%             stem(z_ahu*max(data_ahu),'g','LineWidth',0.5,'Marker','none')
            plot(diff1,'b','LineWidth',1.5)
            plot(diff2,'g','LineWidth',1.5)
%             plot([zeros(1,15) diff1],'b','LineWidth',1.5)
%             stem(50*([zeros(1,15) diff1]>=mean(diff1)+1*std(diff1)),'g','LineWidth',0.5,'Marker','none')
            title(sprintf('%s vs %s',ahus(n).name(1:end-4),vavs(m).name(1:end-4)));
            legend('vav','ahu','vav\_kf\_res','ahu\_kf\_res','FontSize', 12);
            pause
        end
        
        % AHU col 3 - SupplyAirPress,  5 - SupplyFanSpeedOutput
        % VAV col 1 - AirFlowNormalized, 2 - AirValvePosition
        cur_corr = corrcoef(data_ahu, data_vav);
        cur_corr = cur_corr(1,2);
        vav_corr(n) = abs(cur_corr);

        cur_sim = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav));
        vav_sim(n) = cur_sim;

%         vav_score(n) = matched_power_score(14, vav_edge{m}, data_vav, data_ahu); %k=14
        vav_score(n) = dot(diff1, diff2)/(norm(diff1)*norm(diff2));
%         vav_score1(n) = abs( dot(z_ahu, z_vav) ) / (norm(z_ahu)*norm(z_vav)); 
    end
    
    if ahu_list(vav_corr==max(vav_corr)) == ahuid
        ctr = ctr + 1;
    end
    
%     vav_sim = vav_score;

    %TBD:weighted sum of cossim and mps
%     vav_score = vav_score/max(vav_score);
%     vav_sim = w * vav_sim + (1-w) * vav_score;
    
    if ahu_list(vav_sim==max(vav_sim)) == ahuid
        ctr1 = ctr1 + 1;
        correct = [correct; vav_sim', m, find(ahu_list==ahuid), entropy(vav_sim), find(vav_score==max(vav_score),1)];
    else
        [v,i] = sort(vav_sim,'descend');
        flag = 0;
        if ~isempty( find(ahu_list(i(1:k))==ahuid,1) )
            topk = topk + 1;
            flag = 1;
        end
        true = find(ahu_list==ahuid);
        predicted = find(vav_sim==max(vav_sim));
        wrong = [wrong; vav_sim', m, true, predicted, entropy(vav_sim), flag, vav_sim(predicted)-vav_sim(true)];

        max_in_all = ismember( true, find(vav_score==max(vav_score)) );
%         max_in_all1 = ismember( true, find(vav_score1==max(vav_score1),1) );
        
        topk_score = vav_score(i(1:k));
        max_in_topk = vav_score(true)==max(topk_score);
        
        wrong_test = [wrong_test; vav_score', true, predicted, find(vav_score==max(vav_score),1), max_in_all, max_in_topk];
    end
    
    res(m,1:end-1) = vav_sim;
    res(m,end) = find(ahu_list==ahuid);

    fprintf('processed %s\n',vavs(m).name)

end

fprintf('--------------------------------------\n');
fprintf('acc on simple cc is %.4f\n', ctr/num);
fprintf('acc on canny edge seq cossim is %.4f\n', ctr1/num);
fprintf('top_%d rate for miss is %.4f\n', k, topk/num);

% tmp = sum(wrong_test(:,12)| wrong_test(:,13)) + size(correct,1);
tmp = sum(wrong_test(:,end-1)) + size(correct,1);
fprintf('acc on combined is %.4f\n', tmp/num);
% tmp = sum(wrong_test(:,13)) + size(correct,1);
% fprintf('acc on combined1 is %.4f\n', tmp/num);
% tmp = sum(wrong_test(:,13)) + size(correct,1);
% fprintf('acc__ on combined is %.4f\n', tmp/num);

% for i=1:size(ahu_event,1)
%     tmp = ahu_event{i};
%     tmp = reshape(tmp,4*24*7,[])';
%     tmp = sum(tmp,1);
%     ahu_event{i} = tmp;
% end
% for i=1:size(vav_event,1)
%     tmp = vav_event{i};
%     tmp = reshape(tmp,4*24*7,[])';
%     tmp = sum(tmp,1);
%     vav_event{i} = tmp;
% end

% fea_ahu = tfidf(cell2mat(ahu_event));
% fea_vav = tfidf(cell2mat(vav_event));
fea_ahu = tfidf(cell2mat(ahu_kf_res));
fea_vav = tfidf(cell2mat(vav_kf_res));
% fea_ahu = tfidf(cell2mat(ahu_mle_res));
% fea_vav = tfidf(cell2mat(vav_mle_res));

% events = [cell2mat(ahu_event);cell2mat(vav_event)];
% fea = tfidf(events);
% fea_ahu = fea(1:length(ahu_event),:);
% fea_vav = fea(1:length(vav_event),:);

ctr2 = 0;
topk = 0;
for m = 1:num
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    f1 = fea_vav(m,:);

    vav_sim = zeros(length(ahus),1);
    for n = 1:length(ahus)
        fn = [path_ahu, ahus(n).name];
        f2 = fea_ahu(n,:);

        cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
        vav_sim(n) = abs(cur_sim);
    end
    
    if ahu_list(vav_sim==max(vav_sim)) == ahuid
        ctr2 = ctr2 + 1;
    else
        [v,i] = sort(vav_sim,'descend');
        if ~isempty( find(ahu_list(i(1:k))==ahuid,1) )
            topk = topk + 1;
        end      
    end    
end

fprintf('acc on tfidf cossim is %.4f\n', ctr2/num);
fprintf('top_%d rate for miss is %.4f\n', k, topk/num);

%%
figure
hold on
grid on
[p, v] = ecdf(correct(:,end));
plot(v,p,'r', 'LineWidth', 2)
[p, v] = ecdf(wrong(:,end));
plot(v,p,'k', 'LineWidth', 2)
legend('correct\_sim','wrong\_sim', 'Location','southeast', 'FontSize', 12);
