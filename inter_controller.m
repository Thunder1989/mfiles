clear
clc
path_ahu = '/Users/hdz1989/Downloads/code/inference/data_ahu/';
path_vav = '/Users/hdz1989/Downloads/code/inference/data_vav/';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

inter_corr = [];
intra_corr = [];
ctr = 0;
for m = 1:9
% for m = 1:length(vavs)
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    data_vav = csvread(fn,2);
    
    vav_corr = zeros(length(ahus),1);
    ahu_list = zeros(length(ahus),1);
    for n = 1:length(ahus)
        fn = [path_ahu, ahus(n).name];
        cur_ahuid = str2double(ahus(n).name(5));
        data_ahu = csvread(fn,2);
%         plot(data(:,3))
%         plot(data(:,5))
%         plot(data_vav(:,1))
%         plot(data_vav(:,2))
        % AHU col 3 - SupplyAirPress,  5 - SupplyFanSpeedOutput
        % VAV col 1 - AirFlowNormalized, 2 - AirValvePosition
        cur_corr = corrcoef(data_ahu(:,5), data_vav(:,1));
        cur_corr = cur_corr(1,2);
        vav_corr(n) = abs(cur_corr);
        ahu_list(n) = cur_ahuid;
        if cur_ahuid == ahuid
            intra_corr = [intra_corr; cur_corr];
        else
            inter_corr = [inter_corr; cur_corr];
        end
    end
    if ahu_list(vav_corr==max(vav_corr)) == ahuid
        ctr = ctr + 1;
    else
        vavs(m).name
%         figure
%         hold on
%         grid on
%         yyaxis left
%         plot(data_vav(:,1), 'LineWidth',2)
% 
%         yyaxis right
%         [id, n] = max(vav_corr);
%         fn = [path_ahu, ahus(n).name];
%         data_ahu = csvread(fn,2);
%         fprintf('wrong AHU: %s\n', ahus(n).name);
%         plot(data_ahu(:,5), 'r', 'LineWidth',2)
%         
%         n = ahuid;
%         fn = [path_ahu, ahus(n).name];
%         data_ahu = csvread(fn,2);
%         fprintf('correct AHU: %s\n', ahus(n).name);
%         plot(data_ahu(:,5), 'k', 'LineWidth',2)
%         legend({'vav\_AirFlow', 'maxAHU\_supplyFanSpeed', 'trueAHU\_supplyFanSpeed'}, 'Location','southeast', 'FontSize',12)
%         pause
    end
end

fprintf('acc is %.4f\n', ctr/length(vavs));
%%
figure
hold on
grid on
[p, v] = ecdf(intra_corr);
plot(v,p,'r', 'LineWidth', 2)
[p, v] = ecdf(inter_corr);
plot(v,p,'k', 'LineWidth', 2)
legend('intra\_ahu\_corr','inter\_ahu\_corr', 'Location','southeast', 'FontSize',12);

%%
figure
yyaxis left
plot(data(:,5))
yyaxis right
plot(data_vav(:,1))