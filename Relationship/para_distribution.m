function para_distribution(bid, Q_vav)
    % get Q R distribution
    path_ahu = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\ahu_common\');
    path_vav = strcat('D:\TraneData\cut\ahu_property_file_10', num2str(bid) , '_cut\vav_common\');
    ahus = dir(strcat(path_ahu, '*.csv'));
    vavs = dir(strcat(path_vav, '*.csv'));
    ahu_id = zeros(length(vavs),1);
    for n = 1:length(vavs);
        str = regexp(vavs(n).name,'[0-9]+','match');
        ahu_id(n) = str2double(str(1));
    end
    unique(ahu_id)

    cur = ahu_id(1);
    Q_0 = cell(4,1);
    Q_1 = cell(4,1);
    R_0 = cell(4,1);
    R_1 = cell(4,1);
    for n = 2:length(vavs)
        if ahu_id(n) ~= cur
            cur = ahu_id(n);
            plot_all_dim(Q_0);
            plot_all_dim(Q_1);
            pause
            
            Q_0 = cell(4,1);
            Q_1 = cell(4,1);
            continue;
        else
            cur_M = Q_vav{n};
            idx = get_max(cur_M);
            
            %Q for event
            v = cur_M{idx};
            for dim = 1:4
                Q_1{dim} = [Q_1{dim} v(dim)];
            end
            %Q for non-event
            v = cur_M{3-idx};
            for dim = 1:4
                Q_0{dim} = [Q_0{dim} v(dim)];
            end
        end
    end

function idx = get_max(M)
    M_max = max(M{:});
    M_max = M_max(end);
    
    if M{1}(end) == M_max %larger r for event
        idx = 1;
    else
        idx = 2;
    end
    
function plot_all_dim(M)
    figure
    hold on
    for i=1:4
        subplot(2,2,i)
        hist(M{i},20)
    end
