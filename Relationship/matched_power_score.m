function score = matched_power_score(K, edge_vav, data_vav, data_ahu)

    K = min(K, sum(edge_vav));
    edge_idx = find(edge_vav==1);
    %handle edges at the beginning and last
    while edge_idx(end)+1>length(edge_vav)
        edge_idx = edge_idx(1:end-1);
    end
    
    if edge_idx(1)==1
        edge_idx = edge_idx(2:end);
    end

    change = [0 diff(data_vav)'];
    edge_mag = abs(change(edge_idx+1));
    [v, i] = sort(edge_mag,'descend');
    edge_mag_vav = change(edge_idx(i(1:K))+1); %indexing from change to keep the sign of change
    
    edge_mag_ahu = zeros(1,K);
    for k=1:K
        pos = i(k);
        pos = edge_idx(pos);
        edge_mag_ahu(k) = data_ahu(pos+1) - data_ahu(pos-1); %seems pos+1 also works, meaning no time shift allowed
    end
    
    score = sum(edge_mag_vav .* edge_mag_ahu);
    score = score/K;
    if score == 0
        score = -9999;
    end
    
