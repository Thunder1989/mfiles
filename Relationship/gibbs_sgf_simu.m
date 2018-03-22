function [Y, Z, M, Q, R] = gibbs_sgf_simu(Q_ahu, R_ahu, F_switch)
    
    clc
    
    %----generate data----
    rep = 20;
    len = repmat(10,1,rep);
    vel = repmat([0 30 -30 0],1,rep/4);
    state = repmat([1 2 2 1],1,rep/4);

    F = [1 1; 0 1];
    H = [1 0; 0 1];
    
    Q_star = Q_ahu{2};
    R_star = R_ahu{2};
    
    Z_star = zeros(sum(len),1);
    Z_star(11:10:end) = 1;
    Z_star(1:40:end) = 0; %checked

    y_prev = rand*100; %initial Y_0
    
    Y = [];
    X = [];

    %eq: Y = N(FY,Q), X = N(HX,R)
    for i = 1:length(len)
        count = len(i);
        for t = 1:count
            tmp_y = mvnrnd( F*[y_prev vel(i)]',  Q_star{state(i)} );
            Y = [Y; tmp_y];
            y_cur = tmp_y(1);
            vy_cur = tmp_y(2);

            tmp_x = mvnrnd( H*[y_cur vy_cur]', R_star{state(i)} );
            X = [X; tmp_x];
            y_prev = y_cur;
        end

    end
    
    figure
    plot([Y Z_star*max(Y(:))])
    

    %------initialization------
    rnd = rand(size(X,1),1);
    Z = zeros(size(X,1),1);
    Z(rnd>0.5) = 1;
	
    M = [.9, .5; .1, .5]; %[p00, p10; p01, p11]
    F_ = {[1 1; 0 1], [1 0; 0 1]}; %1-steady, 2-transition
    F = F_{1};
    H = [1 0; 0 1];
    
    Q = {rand(2), rand(2)};
 	R = {rand(2), rand(2)};
%     fprintf('initial Q and R:\n');
    class_mapped = map_i(Q);
    for i = 0:1
        if F_switch
            F = F_{class_mapped(i)};
        end
        [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z,F,H);
    end
   
    Z_err = [];
    Q_err = [];
    R_err = [];

    for iter = 1:6
        fprintf('----iter #%d----\n',iter);
        
        %------EM------
        K = 10; % iter# for EM
        N = 200; % samples per iter
        p_data = zeros(K,1);
        for k = 1:K
    %         fprintf('--------------EM iter# %d--------------\n',k);

            %---E step---
            class_mapped = map_i(Q);

            %sample Z
            Z_sample = repmat(Z,1,N);
            for t = 2:length(Z)-1 %todo: shuffle the order
                p_tmp = zeros(2,1);
                for i = 1:length(p_tmp)
                    p_tmp(i) = M( i, Z(t-1)+1 ) * M( Z(t+1)+1, i ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i});% * mvnpdf(Y(t,:)', F*Y(t-1,:)', Q{i});
                end
                p_tmp = p_tmp/sum(p_tmp);

                for n = 1:N
                    if p_tmp(1) > rand(1) %todo: change to find min on cumsum
                        Z_sample(t,n) = 0;
                    else
                        Z_sample(t,n) = 1;
                    end
                end
            end

            %sample Y
            Y_sample = repmat(Y,1,1,N);
            p_tmp = zeros(N,1);
            for t = 2:size(Y,1)-1

                Q_prev = Q{Z(t-1)+1};
                Q_next = Q{Z(t+1)+1};
                Q_cur = Q{Z(t)+1};
                R_cur = R{Z(t)+1};
    %             celldisp(R)
    %             Z(t)
                if F_switch
                    F = F_{class_mapped(Z(t))};
                end

                if class_mapped(Z(t-1))==1 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==1 ...
                    || class_mapped(Z(t-1))==2 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==2 
                    %010 or 101: p(x_t|y_t)
                    Sigma_e = H^-1*R_cur*(H^-1)';
                    Mu_e = H^-1*X(t,:)';
                elseif class_mapped(Z(t-1))==2 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==1 ...
                    || class_mapped(Z(t-1))==1 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==2
                    %110 or 001: p(x_t|y_t) p(y_t|y_t-1)
                    Sigma_e = ( Q_cur^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * Q_cur^-1 * (F*Y(t-1,:)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
                elseif class_mapped(Z(t-1))==2 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==1 ...
                    || class_mapped(Z(t-1))==1 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==2
                    %100 or 011: p(x_t|y_t) p(y_t|y_t+1)
                    Sigma_e = ( (F^-1*Q_next*(F^-1)')^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y(t+1,:)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');    
                end

                %normal case - product of 3 gaussians, always used for sampling position 
                Sigma_ = ( Q_cur^-1 + (F^-1*Q_next*(F^-1)')^-1 )^-1;
                Mu_ = Sigma_* Q_cur^-1 * (F*Y(t-1,:)') + Sigma_ * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y(t+1,:)');
                Sigma = ( Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');

                if Z(t)~=Z(t-1) || Z(t)~=Z(t+1)
                    for n = 1:N
                        Y_sample(t,1,n) = normrnd(Mu(1), Sigma(1,1));
                        Y_sample(t,2,n) = normrnd(Mu_e(2), Sigma_e(2,2));
                    end
                else
                    for n = 1:N
                        Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
                    end
                end

                %data likelihood
                p = 0;
                % for t=2:size(Y,1)-1 
                %     p = p + log( mvnpdf(Y_sample(t,:,n), Y_sample(t-1,:,n)*F', Q) )...
                %         + log( mvnpdf(X(t,:), Y_sample(t,:,n)*H', R) );
                % end
                p_tmp(n) = p;
            end

            p_data(k) = mean(p_tmp);

            %---M step---
            Y = mean(Y_sample,3);
            Z = mode(Z_sample(:,end-10:end),2);
            M = get_M(Z_sample);

            for i = 0:1
                if F_switch
                    F = F_{class_mapped(i)};
                end
                [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z_sample,F,H); %check: use Y_sample?
            end

        end
    
    %     Z = mean(Z_sample(:,end-20:end),2);
        Z = Z_sample;
        
        %calculate errors
        Z = mode(Z(:,11:2:end),2);
        
        assert( length(Z) == length(Z_star) );
        Z_err = [Z_err sum(abs(Z - Z_star))/length(Z)];
        
        tmp = cellfun(@minus, Q, Q_star, 'UniformOutput', false);
        tmp = cellfun(@abs, tmp, 'UniformOutput', false);
        tmp = cell2mat( cellfun(@(x) sum(x(:)), tmp, 'UniformOutput', false) );
        Q_err = [Q_err; tmp(:)'];
        
        tmp = cellfun(@minus, R, R_star, 'UniformOutput', false);
        tmp = cellfun(@abs, tmp, 'UniformOutput', false);
        tmp = cell2mat( cellfun(@(x) sum(x(:)), tmp, 'UniformOutput', false) );
        R_err = [R_err; tmp(:)'];
        
        Z_ = remap_event(Z);
        flip = ~isequal(Z, Z_);
        for i = 1:length(Z)
            if rand>=0.2
                Z(i) = abs(1*flip - Z_star(i));
            end
        end
        
        M = get_M(Z(:));
        for i = 0:1
            if F_switch
                F = F_{class_mapped(i)};
            end
            [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z(:),F,H); %check: use Y_sample?
        end
        
    end
    
    data = X(:,1);
    res = Y;
%     Z = mode(Z(:,11:2:end),2);
    Z = remap_event(Z);

    figure
    hold on
    yyaxis left
    %events in shade
    stem(1:length(Z), Z*max(data), 'Marker','None', 'LineWidth', 5, 'Color', [.8 .8 .8]);
    %original data
    plot(1:length(res), data, 'k--', 'LineWidth', 1.5)
    %filtered data
    %     plot(1:length(res), res(:,1),'g-')
    yyaxis right
    %velocity
    plot(1:length(res), res(:,2), 'r-','LineWidth', 2)

    figure
    subplot(3,1,1)
    plot(Z_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)
    subplot(3,1,2)
    plot(Q_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)
    subplot(3,1,3)
    plot(R_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)

%     M
%     celldisp(Q)
%     celldisp(R)
    

function [Q, R] = get_Q_R(i,X,Y,Z_sample,F,H)
    Q = zeros(2);
    R = zeros(2);
    
    ctr = 0;
    for n = 1:size(Z_sample,2)
        Z_cur = Z_sample(:,n);
        idx = find(Z_cur==i);
        for k = 1:length(idx)
            t = idx(k);
            if t==1
                continue;
            end
            tmp = Y(t,:)' - F*Y(t-1,:)';
            Q = Q + tmp * tmp';

            tmp = X(t,:)' - H*Y(t,:)';
            R = R + tmp * tmp';

            ctr = ctr + 1;
        end
    end
    
    Q = Q/ctr;
    R = R/ctr;


function M = get_M(Z_sample)
    M = zeros(2);
    for n = 1:size(Z_sample,2)
        Z = Z_sample(:,n);
        n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0));
        n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
        % z0 = betarnd(n01+prior.z01, n0-n01+prior.z00);
        % z1 = betarnd(n10+prior.z10, n1-n10+prior.z11);
        z0 = n01/n0;
        z1 = n10/n1;
    	M = M + [1-z0, z1; z0, 1-z1];
    end
    M = M / size(Z_sample,2);
    
    
function class_mapped = map_i(R)

%     fprintf('re-mapping z-0/1 to 1-steady/2-transition...\n')
%     celldisp(R)

    k = [0 1];
    R_max = max(R{:});
    R_max = R_max(1);
    
    if R{1}(1) == R_max %larger r for event
        v = [2 1];
    else
        v = [1 2];
    end

    class_mapped = containers.Map(k,v);
%     class_mapped(0)
%     class_mapped(1)
%     pause

    
function idx = get_event_i(R)

    R_max = max(R{:});
    R_max = R_max(end);
    
    if R{1}(end) == R_max %larger r for event
        idx = 0;
    else
        idx = 1;
    end

    