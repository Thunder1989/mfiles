%% multi-masking using each ahu for each vav
load('320_events.mat');

path_ahu = './data_ahu/';
path_vav = './data_vav/';
ahus = dir(strcat(path_ahu, '*.csv'));
vavs = dir(strcat(path_vav, '*.csv'));

num = length(ahus);
ahu_list = zeros(num,1);
for n = 1:num
    str = regexp(ahus(n).name,'[0-9]+','match');
    ahu_list(n) = str2double(str(1));
end
num = size(vavs,1);
vav_list = zeros(num,1);
for m = 1:num
    str = regexp(vavs(m).name,'[0-9]+','match');
    vav_list(m) = str2double(str(1));
end

ahu_ = cellfun(@transpose,ahu,'UniformOutput',false);
vav_ = cellfun(@transpose,vav,'UniformOutput',false);
ahu_ = cellfun(@round,ahu_,'UniformOutput',false);
vav_ = cellfun(@round,vav_,'UniformOutput',false);
ahu_ = cellfun(@remap_event,ahu_,'UniformOutput',false);
vav_ = cellfun(@remap_event,vav_,'UniformOutput',false);

fea_ahu = cell2mat( ahu_ );
fea_vav = cell2mat( vav_ );

%plot ahu events
% figure
% num = size(ahu_,1);
% for i =1:num
%     subplot(num,1,i)
%     plot(ahu_{i})
% end

%masking vav with each ahu
num = size(fea_vav,1);
ctr = 0;
ahu_list_copy = repmat(ahu_list, 1, length(ahu_list));
assignment = zeros(size(fea_vav,1),1);
for m = 1:num
    ahu_id = vav_list(m);
    idx = find(ahu_list == ahu_id);
    f1 = fea_vav(m, :);

    vav_sim = zeros(length(ahu_list), length(ahu_list));
    for n = 1:length(ahu_list)
%         if n ~= idx
%             continue
%         end
        
        mask = fea_ahu(n,:);
        mask = mask | [false mask(1:end-1)];
        vav_tmp = double(f1 & mask);
        for k = 1:length(ahus)
            f2 = fea_ahu(k,:);
            cur_sim = dot(vav_tmp, f2)/(norm(vav_tmp)*norm(f2)); 
            vav_sim(n, k) = abs(cur_sim);       
        end
    end
    
    assignment(m) = ahu_list_copy( vav_sim==max(max(vav_sim)) );
    pred = ahu_list_copy( vav_sim==max(max(vav_sim)) );
    pred = ahu_list_copy(vav_sim == max(diag(vav_sim)));
    assert ( max( vav_sim(mod( find(vav_sim == max(diag(vav_sim)))-1, length(ahu_list) )+1, :) ) == max(diag(vav_sim)) ); %complicated indexing, lol
    if ismember(ahu_id, pred) && length( find(vav_sim==max(max(vav_sim))) )==1
        ctr = ctr + 1;
    end
    
end

fprintf('acc on tfidf cossim is %.4f\n', ctr/num);

%re-masking using the assigned ahu
% for m = 1:num
%    f1 = fea_vav(m, :);        
%    mask = fea_ahu(ahu_list==assignment(m),:);
%    mask = mask | [false mask(1:end-1)];
%    fea_vav(m, :) = double(f1 & mask);
% end

%% test distribution within each groups from multi-masking
assign_map = cell(max(assignment),1);
for i = 1:length(assignment)
    assign_map{assignment(i)} = [assign_map{assignment(i)} i];
end

corr_map = cell(max(assignment),1);
for i = 1:length(assign_map)
    cur = assign_map{i};
    num = length(cur);
    tmp = zeros(num, num);
    
    for j = 1:num
        
        vav_tmp = fea_vav(cur(j),:);
        for k = 1:j-1
            cor = corrcoef(vav_tmp, fea_vav(cur(k),:));
            tmp(j,k) = cor(1,2);
        end
    
    end
    corr_map{i} = max(tmp, tmp');
end

corr_sum = cellfun(@sum, corr_map, 'UniformOutput', false);
corr_sum = cellfun(@transpose, corr_sum, 'UniformOutput', false);
for i = 1:length(corr_sum)
    corr_sum{i} = sortrows([ corr_sum{i}, vav_list(assign_map{i}) ]);
end

%cal intersection rate for top-k vavs in each ahu group
k = 4;
rate = zeros(length(assign_map),1);
for i = 1:length(assign_map)
    cur = assign_map{i};
    vav_tmp = fea_vav(cur,:);
    ctr = 0;
    for j = 1:size(vav_tmp,2)
        if length( unique(vav_tmp(:,j)) )==1
            ctr = ctr + 1;
        end
    end
    rate(i) = ctr / size(vav_tmp,2);
end
rate

%% re-mask vav with the ahu assignment from multimasking
for m = 1:num
    f1 = fea_vav(m, :);        
    mask = fea_ahu(ahu_list==assignment(m),:);
    mask = mask | [false mask(1:end-1)];
    fea_vav(m, :) = double(f1 & mask);
end

%vertical comparison - kmeans
K=5;
N=2;
vav_sub = fea_vav;
[c_idx,~,~,D] = kmeans(vav_sub', K, 'EmptyAction','drop');
assert(length(c_idx) == size(vav_sub,2));

% horizontal comparison and udpating z - MLE
TH = 0.6;
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
            if conf_tmp <= TH && conf_tmp >= 1-TH && c_idx(i)==k
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
