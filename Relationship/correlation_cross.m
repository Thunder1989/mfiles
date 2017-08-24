clear
clc

T = 7*4; % # of days

path_ahu = 'D:\TraneData\cut\ahu_property_file_10606_cut\ahu_common\';
path_vav = 'D:\TraneData\cut\ahu_property_file_10606_cut\vav_common\';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

fn = [path_ahu, ahus(1).name];
cur_data = csvread(fn,1);
ahu_measure_num = size(cur_data,2);
fn = [path_vav, vavs(1).name];
cur_data = csvread(fn,1);
vav_measure_num = size(cur_data,2);

N = 5*2; %Kalman Filter lookback window size
% num = 21;
num = length(ahus);
ahu_event = cell(num,ahu_measure_num);
ahu_kf_res = cell(num,ahu_measure_num);
ahu_sgf = cell(num,ahu_measure_num);
ahu_sgf_res = cell(num,ahu_measure_num);
ahu_mle_res = cell(num,ahu_measure_num);
ahu_list = zeros(num,1);
for n = 1:num
    fprintf('processing %s\n',ahus(n).name)
    str = regexp(ahus(n).name,'[0-9]+','match');
    cur_ahuid = str2double(str(1));
    ahu_list(n) = cur_ahuid;
    
    fn = [path_ahu, ahus(n).name];
    ahu_data = csvread(fn,1);
    for m = 1:ahu_measure_num
        cur_data = ahu_data(1:4*24*T,m);

        delta = [0 diff(cur_data)'];
        e_ahu = edge( repmat(cur_data',3,1), 0.4 ); %th = 0.4
        e_ahu = e_ahu(1,:);
        e_ahu = e_ahu | [false e_ahu(1:end-1)] | [e_ahu(2:end) false];
    %     e_ahu = delta .* double(e_ahu);
%         assert(sum(e_ahu)~=0);
        ahu_event{n,m} = double(e_ahu);

        kf = StandardKalmanFilter(cur_data',8,N,'EWMA'); 
    %     kf = gibbs_hmm_uni(data_ahu,0);
        diff2 = abs(cur_data' - kf);
        diff2(isnan(diff2)) = 0;
        ahu_kf_res{n,m} = diff2(2:end-1); %TBD: make the manual period self-deciding

    %     [mle, x] = data_mle(data_ahu',1); 
    %     diff2 = abs(data_ahu - mle(:));
    %     diff2(isnan(diff2)) = 0;
    %     ahu_mle_res{n} = diff2;

    %     [y,z,p] = gibbs_sgf(data_ahu,0);
    %     ahu_sgf{n} = z;
    %     diff2 = abs(data_ahu - y(:,1));
    %     diff2(isnan(diff2)) = 0;
    %     ahu_sgf_res{n} = diff2;

    end
end

num = length(vavs);
vav_edge = cell(num,vav_measure_num);
vav_event = cell(num,vav_measure_num);
vav_kf_res = cell(num,vav_measure_num);
vav_sgf = cell(num,vav_measure_num);
vav_sgf_res = cell(num,vav_measure_num);
vav_mle_res = cell(num,vav_measure_num);
wrong_test = [];
score_tmp = [];
w = 0.5;
debug = 0;
vav_to_ahu_mapping = zeros(num,1);
for n = 1:num
    fprintf('processing %s\n',vavs(n).name)
    str = regexp(vavs(n).name,'[0-9]+','match');
    vav_to_ahu_mapping(n) = str2double(str(1));
    
    fn = [path_vav, vavs(n).name];
    vav_data = csvread(fn,1);
    for m = 1:vav_measure_num
        cur_data = vav_data(1:4*24*T,m);

        delta = [0 diff(cur_data)'];
        e_vav = edge( repmat(cur_data',3,1), 1.25 ); %th = 1.25 for 320 596, 0.6 for 642
        e_vav = e_vav(1,:);
        vav_edge{n,m} = double(e_vav);
        e_vav = e_vav | [false e_vav(1:end-1)] | [e_vav(2:end) false];
    %     e_vav = delta .* double(e_vav);
        e_vav = double(e_vav);
%         assert(sum(e_vav)~=0);
        vav_event{n,m} = e_vav;

        kf = StandardKalmanFilter(cur_data',8,N,'EWMA'); 
    %     kf = gibbs_hmm_uni(data_vav,0);
        diff1 = abs(cur_data' - kf);
        diff1(isnan(diff1)) = 0;
        diff1 = diff1(2:end-1);
        vav_kf_res{n,m} = diff1;

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
    
    end
end

k = 3;
num = length(vavs);
acc = zeros(vav_measure_num, ahu_measure_num);
for vav_measure_id = 1:vav_measure_num

    for ahu_measure_id = 1:ahu_measure_num
%     for ahu_measure_id = 1:1

        correct = [];
        wrong = [];
        wrong_test = [];
        ctr = 0;
        topk = 0;
        for vav_id = 1:length(vavs)
        
            ahuid = vav_to_ahu_mapping(vav_id);
            e_vav = vav_event{vav_id, vav_measure_id};
            diff1 = vav_kf_res{vav_id, vav_measure_id};

            vav_sim = zeros(length(ahus),1);
            vav_score = zeros(length(ahus),1);
            for ahu_id = 1:length(ahus)
                e_ahu = ahu_event{ahu_id, ahu_measure_id};
                diff2 = ahu_kf_res{ahu_id, ahu_measure_id};

                if sum(e_ahu)==0 || sum(e_vav)==0
                    vav_sim(ahu_id) = 0;
                else
                    cur_sim = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav));
                    vav_sim(ahu_id) = cur_sim;
                end

%                 vav_score(n) = matched_power_score(14, vav_edge{m}, data_vav, data_ahu); %k=14
                if sum(diff1)==0 || sum(diff2)==0
                    vav_score(ahu_id) = 0;
                else
                    vav_score(ahu_id) = dot(diff1,diff2) / ( norm(diff1) * norm(diff2) );
%                 vav_score1(n) = abs( dot(z_ahu, z_vav) ) / (norm(z_ahu)*norm(z_vav)); 
                end
            end            
            
            if ismember(ahuid, ahu_list(vav_sim==max(vav_sim))) && max(vav_sim)~=0
                ctr = ctr + 1;
                correct = [correct; vav_sim', vav_id, find(ahu_list==ahuid), find(vav_score==max(vav_score),1)];
            else
                [v,i] = sort(vav_sim,'descend');
                flag = 0;
                if ~isempty( find(ahu_list(i(1:k))==ahuid,1) )
                    topk = topk + 1;
                    flag = 1;
                end
                true = find(ahu_list==ahuid);
                predicted = find(vav_sim==max(vav_sim),1);
                wrong = [wrong; vav_sim', vav_id, true, predicted, flag, vav_sim(predicted)-vav_sim(true)];

                max_in_all = ismember( true, find(vav_score==max(vav_score)) ) & length( find(vav_score==max(vav_score)) ) < length(ahus);
        %         max_in_all1 = ismember( true, find(vav_score1==max(vav_score1),1) );

                topk_score = vav_score(i(1:k));
                max_in_topk = vav_score(true)==max(topk_score);

                wrong_test = [wrong_test; vav_score', true, predicted, find(vav_score==max(vav_score),1), max_in_all, max_in_topk];
                
            end
        end
        
%             fprintf('--------------------------------------\n');
%             fprintf('acc on simple cc is %.4f\n', ctr/num);
%             fprintf('acc on canny edge seq cossim is %.4f\n', ctr/num);
%             fprintf('top_%d rate for miss is %.4f\n', k, topk/num);

%             tmp = sum(wrong_test(:,12)| wrong_test(:,13)) + size(correct,1);
        if isempty(wrong_test)
            acc(vav_measure_id, ahu_measure_id) = size(correct,1)/num;
        else
            tmp = sum(wrong_test(:,end-1)) + size(correct,1);
            acc(vav_measure_id, ahu_measure_id) = tmp/num;
        end
%             fprintf('acc on combined is %.4f\n', tmp/num);
%             tmp = sum(wrong_test(:,13)) + size(correct,1);
%             fprintf('acc on combined1 is %.4f\n', tmp/num);

    end
end
