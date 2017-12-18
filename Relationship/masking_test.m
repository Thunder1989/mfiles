%% multi-masking using each ahu for each vav
load('320_events.mat');

num = length(ahus);
ahu_list = zeros(num,1);
for n = 1:num
    str = regexp(ahus(n).name,'[0-9]+','match');
    ahu_list(n) = str2double(str(1));
end
num = size(vav_,1);
vav_list = zeros(num,1);
for m = 1:num
    str = regexp(vavs(m).name,'[0-9]+','match');
    vav_list(m) = str2double(str(1));
end

%ahu masks
% figure
% num = size(ahu_,1);
% for i =1:num
%     subplot(num,1,i)
%     plot(ahu_{i})
% end

ahu_ = cellfun(@transpose,ahu,'UniformOutput',false);
vav_ = cellfun(@transpose,vav,'UniformOutput',false);
ahu_ = cellfun(@round,ahu_,'UniformOutput',false);
vav_ = cellfun(@round,vav_,'UniformOutput',false);
ahu_ = cellfun(@remap_event,ahu_,'UniformOutput',false);
vav_ = cellfun(@remap_event,vav_,'UniformOutput',false);

fea_ahu = cell2mat( ahu_ );
fea_vav = cell2mat( vav_ );

%masking vav with each ahu
num = size(fea_vav,1);
ctr = 0;
ahu_list_copy = repmat(ahu_list, 1, length(ahu_list));
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
    
    if ismember( ahu_id, ahu_list_copy( vav_sim==max(max(vav_sim)) ) ) && length( find(vav_sim==max(max(vav_sim))) )==1
        ctr = ctr + 1;
    end
        
end

fprintf('acc on tfidf cossim is %.4f\n', ctr/num);