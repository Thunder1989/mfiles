function [Y,Z,M] = gibbs_sgf(data, F_switch, debug)

    %------initialization------
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,5);
    Y2 = EWMA(X2,3);
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];
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
    class_mapped = map_i(R);
    for i = 0:1
        if F_switch
            F = F_{class_mapped(i)};
        end
        [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z,F,H);
    end

%     [g_mu, g_sigma] = get_global_para(get_event_i(R),Y,Z);
    
	%------EM------
    K = 10; % iters for EM
    N = 100; % samples per iter
    p_data = zeros(K,1);
    for k = 1:K
        fprintf('--------------iter# %d--------------\n',k);

        %E step
        class_mapped = map_i(R);
        
        %sample Z
        Z_sample = repmat(Z,1,N);
        for t = 2:length(Z)-1 %todo: shuffle the order
            p_tmp = zeros(2,1);
            for i = 1:length(p_tmp)
                p_tmp(i) = M( i, Z(t-1)+1 ) * M( Z(t+1)+1, i ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i});
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
            
%             if Z(t)~=Z(t-1) || Z(t)~=Z(t+1)
            if class_mapped(Z(t-1))==1 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==1 ...
                || class_mapped(Z(t-1))==2 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==2
                %010 or 101: p(x_t|y_t)
                Sigma = H^-1*R_cur*(H^-1)';
                Mu = H^-1*X(t,:)';
            elseif class_mapped(Z(t-1))==2 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==1 ...
                || class_mapped(Z(t-1))==1 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==2
                %110 or 001: p(x_t|y_t) p(y_t|y_t-1)
                Sigma = (Q_cur^-1 + (H^-1*R_cur*(H^-1)')^-1)^-1;
                Mu = Sigma * Q_cur^-1 * (F*Y(t-1,:)') + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
            elseif class_mapped(Z(t-1))==2 && class_mapped(Z(t))==1 && class_mapped(Z(t+1))==1 ...
                || class_mapped(Z(t-1))==1 && class_mapped(Z(t))==2 && class_mapped(Z(t+1))==2
                %100 or 011: p(x_t|y_t) p(y_t|y_t+1)
                Sigma = ((F^-1*Q_next*(F^-1)')^-1 + (H^-1*R_cur*(H^-1)')^-1)^-1;
                Mu = Sigma * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y(t+1,:)') + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');    
            else
                %normal case - product of 3 gaussians
                Sigma_ = (Q_cur^-1 + (F^-1*Q_next*(F^-1)')^-1)^-1;
                Mu_ = Sigma_* Q_cur^-1 * (F*Y(t-1,:)') + Sigma_ * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y(t+1,:)');
                Sigma = (Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1)^-1;
                Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
            end
            
            for n = 1:N
                Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
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
        
        %M step
        Y = mean(Y_sample,3);
        Z = mode(Z_sample(:,end-10:end),2);
        M = get_M(Z_sample);
	
        for i = 0:1
            if F_switch
                F = F_{class_mapped(i)};
            end
            [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z_sample,F,H); %check: use Y_sample?
        end
        
%         [g_mu, g_sigma] = get_global_para(get_event_i(R),Y,Z);

        if debug==1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            stem(Z*50,'r--','LineWidth',1,'Marker','None')
            pause(0.1)
        end
        
    end
    
%     figure
%     plot(p_data,'k','LineWidth',2)

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


function [g_mu, g_sigma] = get_global_para(i,Y,Z)
    
%     Z_prev = [Z(1); Z(1:end-1)];
%     Z_next = [Z(2:end); Z(end)];
%     y_temp = Y( Z_prev~=Z | Z_next~=Z, 2);
    y_temp = Y(Z==i,2);
    g_mu = mean(y_temp(:));
    g_sigma = std(y_temp(:));
    
    
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
    R_max = R_max(end);
    
    if R{1}(end) == R_max %larger r for event
        v = [2 1];
    else
        v = [1 2];
    end

    class_mapped = containers.Map(k,v);

    
function idx = get_event_i(R)

    R_max = max(R{:});
    R_max = R_max(end);
    
    if R{1}(end) == R_max %larger r for event
        idx = 0;
    else
        idx = 1;
    end

    