function [Y,Z,M,Q,R] = gibbs_sgf_K_para(data, Ks, F_switch, debug)

    %------initialization------
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,5);
    Y2 = EWMA(X2,4);
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];
    Z = randi(Ks, size(X,1), 1); %initialize z seq
	
    M = rand(3);
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
    K = 20; % iters for EM
    N = 200; % samples per iter
    p_data = zeros(K,1);
    for k = 1:K
%         fprintf('--------------iter# %d--------------\n',k);

        %---E step---
        non_event_idx = get_non_event_i(R);
        
        %sample Z
        Z_sample = repmat(Z,1,N);
        for t = 2:length(Z)-1 %todo: shuffle the order
            p_tmp = zeros(Ks,1);
            for i = 1:length(p_tmp)
                p_tmp(i) = M( i, Z(t-1) ) * M( Z(t+1), i ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i});% * mvnpdf(Y(t,:)', F*Y(t-1,:)', Q{i});
            end
            p_tmp = p_tmp/sum(p_tmp);

            for n = 1:N
                Z_sample(t,n) = find( cumsum(p_tmp)>=rand(1), 1 );
            end
        end
        
        %sample Y
%         p_tmp = zeros(N,1); %data likelihood
        F = mat2cell( repmat([1 1; 0 1],size(Y,1),1), 2*ones(1,size(Y,1)), 2 );
        if F_switch
            for t=1:size(Y,1)
                if Z(t)==non_event_idx
                    F{t} = F_{1};
                else
                    F{t} = F_{2};
                end
            end
        end
        
        F_even = F(2:2:end);
        X_even = X(2:2:end,:);
        Y_even = Y(2:2:end,:);
        Z_even = Z(2:2:end);

        F_odd  = F(1:2:end);
        X_odd  = X(1:2:end,:);
        Y_odd  = Y(1:2:end,:);
        Z_odd  = Z(1:2:end);
        
        %1 3 5 7
        % 2 4 6
        num_even = size(Y_even,1);
        Y_prev = Y_odd(1:num_even,:);
        Y_next = zeros(num_even,2);
        Y_next(1:size(Y_odd)-1,:) = Y_odd(2:end,:);
        Z_prev = Z_odd(1:num_even);
        Z_next = randi(Ks,num_even,1);
        Z_next(1:size(Z_odd)-1) = Z_odd(2:end);
        Y_even_sample = sample_Y(non_event_idx,N,F_even,H,Q,R,X_even,Y_prev,Y_next,Y_even,Z_prev,Z_next,Z_even);

        % 2 4 6 
        %1 3 5 7
        num_odd = size(Y_odd,1);
        Y_prev = zeros(num_odd,2);
        Y_prev(2:num_odd,:) = Y_even(1:num_odd-1,:);
        Y_next = zeros(num_odd,2);
        Y_next(1:size(Y_even),:) = Y_even;
        Z_prev = randi(Ks,num_odd,1);
        Z_prev(2:end) = Z_even(1:num_odd-1);
        Z_next = randi(Ks,num_odd,1);
        Z_next(1:size(Z_even)) = Z_even; 
        Y_odd_sample = sample_Y(non_event_idx,N,F_odd,H,Q,R,X_odd,Y_prev,Y_next,Y_odd,Z_prev,Z_next,Z_odd);
        
        Y_sample = zeros([size(Y), N]);
        Y_sample(1:2:end,:,:) = Y_odd_sample;
        Y_sample(2:2:end,:,:) = Y_even_sample;
       
%         p_data(k) = mean(p_tmp);
        
        %---M step---
        Y = mean(Y_sample,3);
        Z = mode(Z_sample(:,end-10:end),2);
        M = get_M(Z_sample, Ks);
	
        for i = 1:Ks %TBD: parallelize
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

function Y_ = sample_Y(non_event_idx,N,F,H,Q,R,X,Y_prev,Y_next,Y,Z_prev,Z_next,Z)
    %FXYZ will be sliced, i.e, either evens or odds
    
    Y_ = zeros([size(Y),N]);

    parfor t = 1:size(Y,1)
        
        Q_next = Q{Z_next(t)};
        Q_cur = Q{Z(t)};
        R_cur = R{Z(t)};

        %transition cases
        if Z_prev(t)==non_event_idx && Z(t)~=non_event_idx && Z_next(t)==non_event_idx ...
            || Z_prev(t)~=non_event_idx && Z(t)==non_event_idx && Z_next(t)~=non_event_idx
            %010 or 101: p(x_t|y_t)
            Sigma_e = H^-1*R_cur*(H^-1)';
            Mu_e = H^-1*X(t,:)';
            
        elseif Z_prev(t)~=non_event_idx && Z(t)~=non_event_idx && Z_next(t)==non_event_idx ...
            || Z_prev(t)==non_event_idx && Z(t)==non_event_idx && Z_next(t)~=non_event_idx
            %110 or 001: p(x_t|y_t) p(y_t|y_t-1)
            Sigma_e = ( Q_cur^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
            Mu_e = Sigma_e * Q_cur^-1 * (F{t}*Y_prev(t,:)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
            
        elseif Z_prev(t)~=non_event_idx && Z(t)==non_event_idx && Z_next(t)==non_event_idx ...
            || Z_prev(t)==non_event_idx && Z(t)~=non_event_idx && Z_next(t)~=non_event_idx
            %100 or 011: p(x_t|y_t) p(y_t|y_t+1)
            Sigma_e = ( (F{t}^-1*Q_next*(F{t}^-1)')^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
            Mu_e = Sigma_e * (F{t}^-1*Q_next*(F{t}^-1)')^-1 * (F{t}^-1*Y_next(t,:)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');    
        end

        %normal case - product of 3 gaussians, always used for sampling position 
        Sigma_ = ( Q_cur^-1 + (F{t}^-1*Q_next*(F{t}^-1)')^-1 )^-1;
        Mu_ = Sigma_* Q_cur^-1 * (F{t}*Y_prev(t,:)') + Sigma_ * (F{t}^-1*Q_next*(F{t}^-1)')^-1 * (F{t}^-1*Y_next(t,:)');
        Sigma = ( Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
        Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');

        if ( Z(t)~=Z_prev(t) || Z(t)~=Z_next(t) ) && ...
            ( Z(t)==non_event_idx || Z_prev(t)==non_event_idx || Z_next(t)==non_event_idx )
            Y_(t,:,:) = [ normrnd(Mu(1), Sigma(1,1), N, 1), normrnd(Mu_e(2), Sigma_e(2,2), N, 1) ]';
        else
            Y_(t,:,:) = mvnrnd(Mu, Sigma, N)';
        end

        %data likelihood
%             p = 0;
%             for t=2:size(Y,1)-1 
%                 p = p + log( mvnpdf(Y_sample(t,:,n), Y_sample(t-1,:,n)*F', Q) )...
%                     + log( mvnpdf(X(t,:), Y_sample(t,:,n)*H', R) );
%             end
%             p_tmp(t) = p;
    end

    
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
            break;
        end
    end
