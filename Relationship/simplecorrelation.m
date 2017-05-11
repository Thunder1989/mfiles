clear
clc

T = 7*4; % # of days

path_ahu = './data_ahu/';
path_vav = './data_vav/';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

%%

inter_corr = [];
intra_corr = [];
ctr = 0;
ctr1 = 0;
k = 3;
topk = 0;
num = length(vavs);
% num = 21;
ahu_event = cell(length(ahus),1);
vav_event = cell(num,1);
ahu_list = zeros(length(ahus),1);

for n = 1:length(ahus)
    fn = [path_ahu, ahus(n).name];
    cur_ahuid = str2double(ahus(n).name(5));
    data_ahu = csvread(fn,1);
    data_ahu = data_ahu(1:4*24*T,end);
    e_ahu = edge(repmat(data_ahu',3,1),0.4);
    e_ahu = e_ahu(1,:);
    e_ahu = double(e_ahu | [false e_ahu(1:end-1)] | [e_ahu(2:end) false]);    
    ahu_event{n} = e_ahu;
    ahu_list(n) = cur_ahuid;
end

for m = 1:num
    fn = [path_vav, vavs(m).name];
    ahuid = str2double(vavs(m).name(5));
    data_vav = csvread(fn,1);
    data_vav = data_vav(1:4*24*T,1);
    e_vav = edge(repmat(data_vav',3,1),1);
    e_vav = e_vav(1,:);
    e_vav = double(e_vav | [false e_vav(1:end-1)] | [e_vav(2:end) false]);
    vav_event{m} = e_vav;

    vav_corr = zeros(length(ahus),1);
    vav_sim = zeros(length(ahus),1);
    for n = 1:length(ahus)
        fn = [path_ahu, ahus(n).name];
        cur_ahuid = str2double(ahus(n).name(5));
        data_ahu = csvread(fn,1);
        data_ahu = data_ahu(1:4*24*T,end);
        e_ahu = ahu_event{n};

%         figure
%         hold on
%         %         plot(data(:,3))
%         plot(data_ahu,'r','LineWidth',1.5)
%         plot(data_vav,'k','LineWidth',1.5)
%         stem(e_vav*max(data_vav),'b','LineWidth',0.5,'Marker','none')
%         stem(e_ahu*max(data_ahu),'g','LineWidth',0.5,'Marker','none')
% %         plot(data_vav(:,2))
%         title(sprintf('%s vs %s',ahus(n).name(1:end-4),vavs(m).name(1:end-4)));
%         pause
        
        % AHU col 3 - SupplyAirPress,  5 - SupplyFanSpeedOutput
        % VAV col 1 - AirFlowNormalized, 2 - AirValvePosition
        cur_corr = corrcoef(data_ahu, data_vav);
        cur_corr = cur_corr(1,2);
        vav_corr(n) = abs(cur_corr);

        cur_sim = dot(e_ahu, e_vav)/(norm(e_ahu)*norm(e_vav)); 
        vav_sim(n) = abs(cur_sim);

        if cur_ahuid == ahuid
            intra_corr = [intra_corr; cur_corr];
        else
            inter_corr = [inter_corr; cur_corr];
        end
    end
    if ahu_list(vav_corr==max(vav_corr)) == ahuid
        ctr = ctr + 1;
    end
    
    if ahu_list(vav_sim==max(vav_sim)) == ahuid
        ctr1 = ctr1 + 1;
    else
        [b,i] = sort(vav_sim,'descend');
        if ~isempty( find(ahu_list(i(1:k))==ahuid) )
            topk = topk + 1;
        end      
%         m
%             vav_sim
%         vavs(m).name
%         ahu_list(vav_sim==max(vav_sim))
    end
    
end

fprintf('acc on simple cc is %.4f\n', ctr/num);
fprintf('acc on canny edge seq cossim is %.4f\n', ctr1/num);
fprintf('top_%d rate for miss is %.4f\n', k, topk/num);


for i=1:size(ahu_event,1)
    tmp = ahu_event{i};
    tmp = reshape(tmp,4*24*7,[])';
    tmp = sum(tmp,1);
    ahu_event{i} = tmp;
end
for i=1:size(vav_event,1)
    tmp = vav_event{i};
    tmp = reshape(tmp,4*24*7,[])';
    tmp = sum(tmp,1);
    vav_event{i} = tmp;
end
fea_ahu = tfidf(cell2mat(ahu_event));
fea_vav = tfidf(cell2mat(vav_event));
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
        [b,i] = sort(vav_sim,'descend');
        if ~isempty( find(ahu_list(i(1:k))==ahuid) )
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
