function [Y,Z,M,Q,R] = gibbs_sgf_K(data, Ks, F_switch, debug)

    %------initialization------
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,5);
    Y2 = EWMA(X2,4);
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];
    Z = randi(Ks, size(X,1), 1); %initialize z seq
	
    M = rand(Ks);
    M = bsxfun(@rdivide, M, sum(M)); %transition matrix
    F_ = { [1 1; 0 1], [1 0; 0 1] }; %1-steady, 2-event
    F = F_{1}; %F will be estimated
    H = [1 0; 0 1];
    
    Q = mat2cell(rand(2*Ks,2), 2*ones(1,Ks), 2);
 	R = mat2cell(rand(2*Ks,2), 2*ones(1,Ks), 2);
    non_event_idx = get_non_event_i(R);
    for i = 1:Ks
        if F_switch
            if i==non_event_idx
                F = F_{1};
            else
                F = F_{2};
            end
        end
        [Q{i}, R{i}] = get_Q_R(i,X,Y,Z,F,H);
    end
    
    %------EM------
    K = 10; % iters for EM
    N = 100; % samples per iter
    p_data = zeros(K,1);
    for k = 1:K
        %fprintf('--------------iter# %d--------------\n',k);
        
        %---E step---
        non_event_idx = get_non_event_i(R);
        
        %sample Z
        Z_sample = repmat(Z,1,N);
        p_tmp = zeros(Ks,1);
        for n = 2:N
            for t = 2:length(Z)-1
                if F_switch
                    if Z_sample(t,n-1)==non_event_idx
                        F = F_{1};
                    else
                        F = F_{2};
                    end
                end

                for i = 1:length(p_tmp)
                    p_tmp(i) = M( i, Z_sample(t-1,n) ) * M( Z_sample(t+1,n-1), i ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i}); %* mvnpdf(Y(t,:)', F*Y(t-1,:)', Q{i});
                end
                p_tmp = p_tmp/sum(p_tmp);
                Z_sample(t,n) = find(cumsum(p_tmp) >= rand(1), 1);
            end
        end
        
        %sample Y
        Y_sample = repmat(Y,1,1,N);
%         p_tmp = zeros(N,1); %for data likelihood
        for n = 2:N
            for t = 2:size(Y,1)-1 %todo: shuffle the order

%                 Q_prev = Q{Z(t-1)+1};
                Q_next = Q{Z(t+1)}; %Question: use latest Z?
                Q_cur = Q{Z(t)};
                R_cur = R{Z(t)};

                if F_switch
                    if Z(t)==non_event_idx
                        F = F_{1};
                    else
                        F = F_{2};
                    end
                end

                if Z(t-1)==non_event_idx && Z(t)~=non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)~=non_event_idx && Z(t)==non_event_idx && Z(t+1)~=non_event_idx
                    %010 or 101: p(x_t|y_t)
                    Sigma_e = H^-1*R_cur*(H^-1)';
                    Mu_e = H^-1*X(t,:)';
                elseif Z(t-1)~=non_event_idx && Z(t)~=non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)==non_event_idx && Z(t)==non_event_idx && Z(t+1)~=non_event_idx
                    %110 or 001: p(x_t|y_t) p(y_t|y_t-1)
                    Sigma_e = ( Q_cur^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * Q_cur^-1 * (F*Y_sample(t-1,:,n)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
                elseif Z(t-1)~=non_event_idx && Z(t)==non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)==non_event_idx && Z(t)~=non_event_idx && Z(t+1)~=non_event_idx
                    %100 or 011: p(x_t|y_t) p(y_t|y_t+1)
                    Sigma_e = ( (F^-1*Q_next*(F^-1)')^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y_sample(t+1,:,n-1)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');    
                end

                %normal case - product of 3 gaussians, always used for sampling position 
                Sigma_ = ( Q_cur^-1 + (F^-1*Q_next*(F^-1)')^-1 )^-1;
                Mu_ = Sigma_* Q_cur^-1 * (F*Y_sample(t-1,:,n)') + Sigma_ * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y_sample(t+1,:,n-1)');
                Sigma = ( Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');

                if ( Z(t)~=Z(t-1) || Z(t)~=Z(t+1) ) && ...
                ( Z(t)==non_event_idx || Z(t-1)==non_event_idx || Z(t+1)==non_event_idx )
                    Y_sample(t,1,n) = normrnd(Mu(1), Sigma(1,1));
                    Y_sample(t,2,n) = normrnd(Mu_e(2), Sigma_e(2,2));
                else
                    Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
                end

                %data likelihood
%                 p = 0;
%                 for t=2:size(Y,1)-1 
%                     p = p + log( mvnpdf(Y_sample(t,:,n), Y_sample(t-1,:,n)*F', Q) )...
%                         + log( mvnpdf(X(t,:), Y_sample(t,:,n)*H', R) );
%                 end
%                 p_tmp(n) = p;
            end
        end
        p_data(k) = mean(p_tmp);
        
        %---M step---
        
        Y = mean(Y_sample, 3);
        Z = mode(Z_sample(:,end-10:end), 2);
        M = get_M(Z_sample, Ks);
	
        for i = 1:Ks
            if F_switch
                if i==non_event_idx
                    F = F_{1};
                else
                    F = F_{2};
                end
            end
            [Q{i}, R{i}] = get_Q_R(i,X,Y,Z_sample,F,H); %check: use Y_sample?
        end
        
        if debug==1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            stem(Z*50,'r--','LineWidth',1,'Marker','None')
            pause(0.1)
        end
        
    end
    
%     Z = mean(Z_sample(:,end-20:end),2);
    Z = Z_sample;


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


function M = get_M(Z_sample,K) %K-class transition matrix, each col sums to 1
    M = zeros(K);
    N = zeros(K,1);
    for n = 1:size(Z_sample,2)
        Z = Z_sample(:,n);
        for i = 1:K
            N(i) = N(i) + length( find(Z(1:end-1)==i) );
            for j = 1:K
                N_ij = length(find(Z(1:end-1)==i & Z(2:end)==j)); 
                M(i,j) = M(i,j) + N_ij;
            end
        end 
    end
    assert (sum(N) == numel(Z_sample) - size(Z_sample,2));
    M = bsxfun(@rdivide,M,N)';
    
    
function class_mapped = map_i(R)

    k = [0 1];
    R_max = min(R{:});
    R_max = R_max(end);
    
    if R{1}(end) == R_max %larger r for event
        v = [2 1];
    else
        v = [1 2];
    end

    class_mapped = containers.Map(k,v);


function idx = get_non_event_i(R)

    R_min = min(cell2mat(cellfun(@(x) x(end), R, 'UniformOutput', false))); 
   
    for i=1:length(R)
        if R{i}(end) == R_min 
            idx = i;
        end
    end    
