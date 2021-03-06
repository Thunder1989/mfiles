function [event_ahu, event_vav, Q_ahu, R_ahu, Q_vav, R_vav] = get_events_all(bid)
    close all
    clc
   
    T = 7*4; % # of days
    
    path_ahu = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\ahu_common\');
    path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav_common\');
    ahus = dir(strcat(path_ahu, '*.csv'));
    vavs = dir(strcat(path_vav, '*.csv'));

    num = length(ahus);
    event_ahu = cell(num,1);
    conf_ahu = cell(num,1);
    Q_ahu = cell(num,1);
    R_ahu = cell(num,1);
    num = length(vavs);
    event_vav = cell(num,1);
    conf_vav = cell(num,1);
    Q_vav = cell(num,1);
    R_vav = cell(num,1);
    for n = 1:length(ahus);
        fprintf('processing %s ...\n',ahus(n).name)
        fn = [path_ahu, ahus(n).name];
        cur_ahu = csvread(fn,1); %skip the 1st row, which is the headers
        cur_ahu = cur_ahu(1:4*24*T,:);
        cur_ahu = cur_ahu(:,end); %ahu last col, vav 1st col
        [~,Z,~,Q,R] = gibbs_sgf(cur_ahu,1,0);
        
        event_ahu{n} = mean(Z(:,11:3:end),2);
        Q_ahu{n} = Q;
        R_ahu{n} = R;
    end
    
    for n = 1:length(vavs);
        fprintf('processed %s ...\n',vavs(n).name);
        fn = [path_vav, vavs(n).name];
        cur_vav = csvread(fn,1); %skip the 1st row, which is the headers
        cur_vav = cur_vav(1:4*24*T,:);
        cur_vav = cur_vav(:,1); %ahu last col, vav 1st col
        [~,Z,~,Q,R] = gibbs_sgf(cur_vav,1,0);
        
        event_vav{n} = mean(Z(:,11:3:end),2);
        Q_vav{n} = Q;
        R_vav{n} = R;
    end