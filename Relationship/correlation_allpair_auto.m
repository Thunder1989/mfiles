clear
clc

T = 7*4; % # of days
bid = 320;
path_ahu = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\ahu_common\');
path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav_common\');
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

fn = [path_ahu, ahus(1).name];
cur_data = csvread(fn,1);
ahu_measure_num = size(cur_data,2);

N = 5*2; %Kalman Filter lookback window size
num = length(ahus);
ahu_event = cell(num,ahu_measure_num);
ahu_kf_res = cell(num,ahu_measure_num);
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
%         assert(sum(e_ahu)~=0);
        ahu_event{n,m} = double(e_ahu);

        kf = StandardKalmanFilter(cur_data',8,N,'EWMA'); 
    %     kf = gibbs_hmm_uni(data_ahu,0);
        diff2 = abs(cur_data' - kf);
        diff2(isnan(diff2)) = 0;
        ahu_kf_res{n,m} = diff2(2:end-1); %TBD: make the manual period self-deciding

    end
end

num = length(vavs);
vav_event = cell(num,1);
vav_kf_res = cell(num,1);
vav_assignment = zeros(num,1);
for n = 1:num

    fprintf('processing %s\n',vavs(n).name)
    str = regexp(vavs(n).name,'[0-9]+','match');
    vav_assignment(n) = str2double(str(1));
    
    fn = [path_vav, vavs(n).name];
    vav_data = csvread(fn,1);
    cur_event = zeros(4*24*T, size(vav_data,2));
    cur_kf_res = zeros(4*24*T-2, size(vav_data,2));
    for m = 1:size(vav_data,2)

        cur_data = vav_data(1:4*24*T, m);
        delta = [0 diff(cur_data)'];
        e_vav = edge( repmat(cur_data',3,1), 1.25 ); %th = 1.25 for 320 596, 0.6 for 642
        e_vav = e_vav(1,:);
        e_vav = e_vav | [false e_vav(1:end-1)] | [e_vav(2:end) false];
        e_vav = double(e_vav);
    %         assert(sum(e_vav)~=0);
        cur_event(:,m) = e_vav;

        kf = StandardKalmanFilter(cur_data',8,N,'EWMA'); 
    %     kf = gibbs_hmm_uni(data_vav,0);
        diff1 = abs(cur_data' - kf);
        diff1(isnan(diff1)) = 0;        
        cur_kf_res(:,m) = diff1(2:end-1);
    end
    vav_event{n} = cur_event;
    vav_kf_res{n} = cur_kf_res;

end

fn = strcat('feature_10', num2str(bid));
save(fn, 'num', 'vav_assignment', 'ahu_list', 'ahu_event', 'ahu_kf_res', 'vav_event', 'vav_kf_res');

%% auto pick a pair for each vav and compute accuracy
clc
bid = 320;
path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav_common\');
vavs = dir(strcat(path_vav, '*.csv'));
fn = strcat('feature_10', num2str(bid));
load(fn);
if bid == 596 %take out column for enthalpy 
    ahu_event(:,2) = []; 
    ahu_kf_res(:,2) = []; 
else
    ahu_event(:,1) = []; 
    ahu_kf_res(:,1) = []; 
end
ahu_num = size(ahu_kf_res,1);
ahu_measure_num = size(ahu_kf_res,2);
ctr = 0;
for vav_id = 1:num
    ahuid = vav_assignment(vav_id);

    cur_vav_kf_res = vav_kf_res{vav_id};
    cur_vav_event = vav_event{vav_id};
    raw1 = cell(size(cur_vav_kf_res,2), ahu_measure_num);
    raw2 = cell(size(cur_vav_kf_res,2), ahu_measure_num);
    for vav_measure_id = 1:size(cur_vav_kf_res,2)

        e_vav = cur_vav_event(:, vav_measure_id);
        diff1 = cur_vav_kf_res(:, vav_measure_id);
        for ahu_measure_id = 1:ahu_measure_num
        
            vav_sim = zeros(ahu_num,1);
            vav_score = zeros(ahu_num,1);
            for ahu_id = 1:ahu_num
            
                e_ahu = ahu_event{ahu_id, ahu_measure_id};
                diff2 = ahu_kf_res{ahu_id, ahu_measure_id};

                if sum(e_ahu)==0 || sum(e_vav)==0
                    vav_sim(ahu_id) = 0;
                else
                    cur_sim = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav));
                    vav_sim(ahu_id) = cur_sim;
                end

                if sum(diff1)==0 || sum(diff2)==0
                    vav_score(ahu_id) = 0;
                else
                    vav_score(ahu_id) = dot(diff1,diff2) / ( norm(diff1) * norm(diff2) );
                end
            end
            raw1{vav_measure_id, ahu_measure_id} = vav_sim;
            raw2{vav_measure_id, ahu_measure_id} = vav_score;
        end
    end
    
    [i,j] = find_cell2(raw1);
    vav_sim = raw1{i,j};
    if ismember(ahuid, ahu_list(vav_sim==max(vav_sim))) && max(vav_sim)~=0 && length( find(vav_sim==max(vav_sim)) ) < ahu_num
%         fprintf('bingo1! - %d\n', length(find(vav_sim==max(vav_sim))) );
        ctr = ctr + 1;
    else
        [i,j] = find_cell2(raw2);
        vav_score = raw2{i,j};
        true = find(ahu_list==ahuid);
        max_in_all = ismember( true, find(vav_score==max(vav_score)) ) & length( find(vav_score==max(vav_score)) ) < ahu_num;
        if max_in_all == 1
%             fprintf('bingo2! - %d\n', length(find(vav_score==max(vav_score))) );
            ctr = ctr + 1;
        else
            fprintf('failed on %s (%d) with pair (%d,%d)\n', vavs(vav_id).name, true, i, j);
            disp(vav_score')
        end
    end
end
fprintf('acc is %.4f\n', ctr/num);
