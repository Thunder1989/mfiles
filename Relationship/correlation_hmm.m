function acc = correlation_hmm(bid, week)

path_ahu = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\ahu\');
path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav\');
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

N = 5*2; %Kalman Filter lookback window size
T = 7*week; % # of days
Q = 2;

ahu_event = cell(length(ahus),1);
ahu_kf_res = cell(length(ahus),1);
ahu_list = zeros(length(ahus),1);
num = length(ahus);
for n = 1:num
%     fprintf('processing %s\n',ahus(n).name)
    fn = [path_ahu, ahus(n).name];
    str = regexp(ahus(n).name,'[0-9]+','match');
    cur_ahuid = str2double(str(1));
    ahu_list(n) = cur_ahuid;
    data_ahu = csvread(fn,1);
    data_ahu = data_ahu(1:4*24*T,end);
    
    data = {data_ahu};
    try
        [p_start, A, phi, loglik] = ChmmGauss(data, Q);
        logp_xn_given_zn = Gauss_logp_xn_given_zn(data{1}, phi);
    %     [~,~, loglik] = LogForwardBackward(logp_xn_given_zn, p_start, A);
        path = LogViterbiDecode(logp_xn_given_zn, p_start, A);
        ahu_event{n} = path-1;
    catch
        fprintf('exception on %s\n',ahus(n).name);
        continue
    end
end

num = length(vavs);
% vav_event = cell(num,1);
ctr = 0;
for m = 1:num
    fprintf('processed %s\n',vavs(m).name);

    fn = [path_vav, vavs(m).name];
    str = regexp(vavs(m).name,'[0-9]+','match');
    ahuid = str2double(str(1));
    data_vav = csvread(fn,1);
    data_vav = data_vav(1:4*24*T,3); %for longer period, 3rd col is airflow

    data = {data_vav};
    try
        [p_start, A, phi, loglik] = ChmmGauss(data, Q);
        logp_xn_given_zn = Gauss_logp_xn_given_zn(data{1}, phi);
    %     [~,~, loglik] = LogForwardBackward(logp_xn_given_zn, p_start, A);
        path = LogViterbiDecode(logp_xn_given_zn, p_start, A);
        e_vav = path-1;
%         figure
%         plot(data_vav,'b')
%         hold on
%         plot(e_vav*max(data_vav) ,'r')       
    catch
        fprintf('exception on %s\n',vavs(m).name);
        continue
    end
    
    vav_sim = zeros(length(ahus),1);
    for n = 1:length(ahus)
        e_ahu = ahu_event{n};
        
        if sum(e_ahu)==0 || sum(e_vav)==0
            vav_sim(n) = 0;
        else
            cur_sim = corrcoef(e_ahu,e_vav);
            vav_sim(n) = abs(cur_sim(1,2));
        end

    end
    
    if ismember( ahuid, ahu_list( vav_sim==max(vav_sim) ) ) && length( find(vav_sim==max(vav_sim)) ) < length(ahus);
        ctr = ctr + 1;
    end

end

acc = ctr/num;
fprintf('----output of %d weeks data-----------------\n',week);
fprintf('acc is %.4f\n', acc);


