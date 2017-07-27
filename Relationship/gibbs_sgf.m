function [Y,Z] = gibbs_sgf(data, debug)

    %------initialization------
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,8);
    Y2 = EWMA(X2,5);
    F = [1 1; 0 1];
    H = [1 0; 0 1];
    
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];
    rnd = rand(size(X,1),1);
    Z = zeros(size(X,1),1);
    Z(rnd>0.5) = 1;
	M = [.9, .5; .1, .5]; %[p00, p10; p01, p11]
    
    Q = {ones(2) + 10e-10, ones(2) + 10e-10};
 	R = {ones(2) + 10e-10, ones(2) + 10e-10};
    fprintf('initial Q and R:\n')
    for i = 0:1
	    idx = find(Z==i);
        if ~isempty(idx)
            [Q{i+1}, R{i+1}] = get_Q_R(idx,X,Y,F,H);
        end
    end

	%------EM------
    K = 8; % iters for EM
    N = 100; % samples per iter
    p_data = zeros(K,1);
    for k = 1:K
        fprintf('--------------iter #%d--------------\n',k);

        %E step
        %sample Z
        Z_sample = repmat(Z,1,N);
        for t = 2:length(Z)-1 %todo: shuffle the order
        	p_tmp = zeros(2,1);
        	for i = 0:length(p_tmp)-1
        		p_tmp(i+1) = M( i+1, Z(t-1)+1 ) * M( Z(t+1)+1, i+1 ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i+1});
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
        	Q_cur = Q{Z(t)+1};
        	R_cur = R{Z(t)+1};
            Sigma_ = (Q_prev^-1 + (F^-1*Q_cur*(F^-1)')^-1)^-1;
            Mu_ = Sigma_* Q_prev^-1 * (F*Y(t-1,:)') + Sigma_ * (F^-1*Q_cur*(F^-1)')^-1 * (F^-1*Y(t+1,:)');
            Sigma = (Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1)^-1;
            Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
            
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
        Z = mode(Z_sample,2); %todo: better update Z
        M = get_M(Z);
	
        for i = 0:1
            [Q{i+1}, R{i+1}] = get_Q_R(i,X,Y,Z_sample,F,H);
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
    
    figure
    plot(p_data,'k','LineWidth',2)


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
    Q = Q/ctr
    R = R/ctr

    
function M = get_M(Z)
	n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0));
	n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
	% z0 = betarnd(n01+prior.z01, n0-n01+prior.z00);
	% z1 = betarnd(n10+prior.z10, n1-n10+prior.z11);
	z0 = n01/n0;
	z1 = n10/n1;
	M = [1-z0, z1; z0, 1-z1];
    