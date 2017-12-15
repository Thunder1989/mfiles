% global-local check
close all
clc
load('320_events.mat');
ahu_ = cellfun(@transpose,ahu,'UniformOutput',false);
vav_ = cellfun(@transpose,vav,'UniformOutput',false);
ahu_conf = cell2mat(ahu_);
vav_conf = cell2mat(vav_);
ahu_ = cellfun(@round,ahu_,'UniformOutput',false);
vav_ = cellfun(@round,vav_,'UniformOutput',false);
ahu_ = cellfun(@remap_event,ahu_,'UniformOutput',false);
vav_ = cellfun(@remap_event,vav_,'UniformOutput',false);
ahu_ = cell2mat(ahu_);
vav_ = cell2mat(vav_);

%take ahu1 and ahu7 data
subset = [7,9];
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

    vav_sim = zeros(size(ahu_tmp,1),1);
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
K = 3; %num of topics
N = 2; %num of states for KF output
[c_idx,~,~,D] = kmeans(vav_sub', K);
assert(length(c_idx) == size(vav_sub,2));
%-visualize-
[c,idx] = sort(c_idx);
figure
imagesc(vav_sub(:,idx))
hold on
ts = [0; diff(c)];
x = find(ts~=0);
% stem(x-0.5, ones(size(x))*num+0.5, 'r','Marker','None','LineWidth',4)

D = min(D,2);
tmp = [c_idx, D];
[tmp,idx_] = sortrows(tmp);
figure
imagesc(vav_sub(:,idx_))
hold on
stem(x-0.5, ones(size(x))*num+0.5, 'r','Marker','None','LineWidth',4)

%horizontal comparison - MLE
TH = 0.7;
vav_sub_updated = vav_sub;
num_updated = zeros(size(vav_sub_updated,1),1);
assert( isequal(vav_sub_updated, vav_sub) );
distribution = zeros(num, K*N);
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
        distribution(m,(k-1)*N+1:k*N) = p_tmp;
        
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
assert( isequal(num_updated, sum(vav_sub_updated~=vav_sub, 2) ) );
figure
imagesc(distribution)

%acc eval
fea_vav = vav_sub_updated;
ctr = 0;
for m = 1:num
    ahu_id = vav_label(m);
    f1 = fea_vav(m,:);

    vav_sim = zeros(size(ahu_tmp,1),1);
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

%acc eval
% fea_vav = fea_vav(:,idx(c<K)); %reordered by cluster id
% fea_ahu = fea_ahu(:,idx(c<K));
fea_vav = fea_vav(:,c~=K); %keep the original order
fea_ahu = fea_ahu(:,c~=K);
ctr = 0;
for m = 1:num
    ahu_id = vav_label(m);
    f1 = fea_vav(m,:);

    vav_sim = zeros(size(ahu_tmp,1),1);
    for n = 1:size(ahu_tmp,1)
        f2 = fea_ahu(n,:);

        cur_sim = dot(f1, f2)/(norm(f1)*norm(f2)); 
        vav_sim(n) = abs(cur_sim);
    end
    
    if ismember( ahu_id, ahu_list(vav_sim==max(vav_sim)) ) && length( find(vav_sim==max(vav_sim)) ) < length(vav_sim)
        ctr = ctr + 1;
    end
end
fprintf('acc after taking out random period is %.4f\n', ctr/num);