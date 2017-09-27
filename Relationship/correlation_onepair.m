function acc = correlation_onepair(bid, week)

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
    data_ahu = csvread(fn,1);
    data_ahu = data_ahu(1:4*24*T,end);

%     delta = [0 diff(data_ahu)'];
    e_ahu = edge(repmat(data_ahu',3,1), 0.4); %th = 0.4
    e_ahu = e_ahu(1,:);
    e_ahu = e_ahu | [false e_ahu(1:end-1)] | [e_ahu(2:end) false];    
%     e_ahu = delta .* double(e_ahu);
%     assert(sum(e_ahu)~=0);

    ahu_event{n} = double(e_ahu);
    ahu_list(n) = cur_ahuid;
  
%     res = mmpp(reshape(data_ahu,4*24,[]), []);
%     tmp = res.Z(:,:,end);
%     ahu_kf_res{n} = tmp(:);
    kf = StandardKalmanFilter(data_ahu',M,N,'EWMA'); 
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
correct = [];
wrong = [];
wrong_test = [];

ctr = 0;
ctr1 = 0;
k = 3;
topk = 0;
w = 0.5;
for m = 1:num
%     fprintf('processed %s\n',vavs(m).name);

    fn = [path_vav, vavs(m).name];
    fn = regexprep(fn,'_cut','_all');
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahuid = str2double(str(1));
    data_vav = csvread(fn,1);
    if bid ~= 320
        data_vav = data_vav(1:4*24*T,1); %for longer period, 3rd col is airflow
    else
        data_vav = data_vav(1:4*24*T,3); %for longer period, 3rd col is airflow
    end
%     delta = [0 diff(data_vav)'];
    e_vav = edge(repmat(data_vav',3,1), 1.25); %th = 1.25 for 320 596, 0.6 for 642
    e_vav = e_vav(1,:);
    vav_edge{m} = double(e_vav);
    e_vav = e_vav | [false e_vav(1:end-1)] | [e_vav(2:end) false];
%     e_vav = delta .* double(e_vav);
    e_vav = double(e_vav);
%     assert(sum(e_vav)~=0);
    vav_event{m} = e_vav;

    vav_corr = zeros(length(ahus),1);
    vav_sim = zeros(length(ahus),1);
    vav_score = zeros(length(ahus),1);
    
    kf = StandardKalmanFilter(data_vav',M,N,'EWMA'); 
%     kf = gibbs_hmm_uni(data_vav,0);
    diff1 = abs(data_vav' - kf);
    diff1(isnan(diff1)) = 0;
    diff1 = diff1(2:end-1);
    vav_kf_res{m} = diff1;

%     res = mmpp(reshape(data_vav,4*24,[]), []);
%     tmp = res.Z(:,:,end);
%     vav_kf_res{n} = tmp(:);
%     diff1 = tmp(:);

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
        fn = regexprep(fn,'_cut','_all');
        data_ahu = csvread(fn,1);
        data_ahu = data_ahu(1:4*24*T,end);
        e_ahu = ahu_event{n};
        diff2 = ahu_kf_res{n};
%         z_ahu = ahu_event_ondiff{n};
%         z_ahu = ahu_sgf{n};
              
        % AHU col 3 - SupplyAirPress,  5 - SupplyFanSpeedOutput
        % VAV col 1 - AirFlowNormalized, 2 - AirValvePosition
        cur_corr = corrcoef(data_ahu, data_vav);
        cur_corr = cur_corr(1,2);
        vav_corr(n) = abs(cur_corr);

        if sum(e_ahu)==0 || sum(e_vav)==0
            vav_sim(n) = 0;
        else
            cur_sim = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav));
            vav_sim(n) = cur_sim;
        end

        if sum(diff1)==0 || sum(diff2)==0
            vav_score(n) = 0;
        else
            vav_score(n) = dot(diff1,diff2) / ( norm(diff1) * norm(diff2) );
%         vav_score(n) = matched_power_score(14, vav_edge{m}, data_vav, data_ahu); %k=14
%         vav_score1(n) = abs( dot(z_ahu, z_vav) ) / (norm(z_ahu)*norm(z_vav)); 
        end

    end
    
    if ahu_list(vav_corr==max(vav_corr)) == ahuid
        ctr = ctr + 1;
    end
    
    if ahu_list(vav_sim==max(vav_sim)) == ahuid
        ctr1 = ctr1 + 1;
        correct = [correct; vav_sim', m, find(ahu_list==ahuid), find(vav_score==max(vav_score),1)];
    else
        [v,i] = sort(vav_sim,'descend');
        flag = 0;
        if ~isempty( find(ahu_list(i(1:k))==ahuid,1) )
            topk = topk + 1;
            flag = 1;
        end
        true = find(ahu_list==ahuid);
        predicted = find(vav_sim==max(vav_sim),1);
        wrong = [wrong; vav_sim', m, true, predicted, flag, vav_sim(predicted)-vav_sim(true)];

        max_in_all = ismember( true, find(vav_score==max(vav_score)) ) & length( find(vav_score==max(vav_score)) ) < length(ahus);
%         max_in_all1 = ismember( true, find(vav_score1==max(vav_score1),1) );
        
        topk_score = vav_score(i(1:k));
        max_in_topk = vav_score(true)==max(topk_score);
        
        wrong_test = [wrong_test; vav_score', true, predicted, find(vav_score==max(vav_score),1), max_in_all, max_in_topk];
%         if max_in_all == 1
%             fprintf('paused\n')
%             pause
%         end
%         
    end   
end

fprintf('----output of %d weeks data on 10%d-----------------\n',week, bid);
fprintf('acc on simple cc is %.4f\n', ctr/num);
fprintf('acc on canny edge seq cossim is %.4f\n', ctr1/num);
fprintf('top_%d rate for miss is %.4f\n', k, topk/num);

% tmp = sum(wrong_test(:,12)| wrong_test(:,13)) + size(correct,1);
tmp = sum(wrong_test(:,end-1)) + size(correct,1);
fprintf('acc on combined is %.4f\n', tmp/num);
acc = tmp/num;

fid = fopen('C:\Users\dzhon\Dropbox\acc_vs_time.txt','a');
fprintf(fid, 'acc on %d with %d weeks data - %0.4f\n', bid, week, acc);
fclose(fid);

% fea_ahu = tfidf(cell2mat(ahu_kf_res));
% fea_vav = tfidf(cell2mat(vav_kf_res));
% 
% ctr2 = 0;
% topk = 0;
% for m = 1:num
%     fn = [path_vav, vavs(m).name];
%     ahuid = str2double(vavs(m).name(5));
%     f1 = fea_vav(m,:);
% 
%     vav_sim = zeros(length(ahus),1);
%     for n = 1:length(ahus)
%         fn = [path_ahu, ahus(n).name];
%         f2 = fea_ahu(n,:);
% 
%         cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
%         vav_sim(n) = abs(cur_sim);
%     end
%     
%     if ahu_list(vav_sim==max(vav_sim)) == ahuid
%         ctr2 = ctr2 + 1;
%     else
%         [v,i] = sort(vav_sim,'descend');
%         if ~isempty( find(ahu_list(i(1:k))==ahuid,1) )
%             topk = topk + 1;
%         end      
%     end    
% end
% 
% fprintf('acc on tfidf cossim is %.4f\n', ctr2/num);
% fprintf('top_%d rate for miss is %.4f\n', k, topk/num);

function cdf(correct, wrong)
    figure
    hold on
    grid on
    [p, v] = ecdf(correct(:,end));
    plot(v,p,'r', 'LineWidth', 2)
    [p, v] = ecdf(wrong(:,end));
    plot(v,p,'k', 'LineWidth', 2)
    legend('correct\_sim','wrong\_sim', 'Location','southeast', 'FontSize', 12);
