function output = gibbs_hmm(data, debug)
    %initialization
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,8);
    Y2 = EWMA(X2,5);
    F = [1 1; 0 1];
    H = [1 0; 0 1];
    
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];

    fprintf('initial Q and R:\n')
%     Q = cov(Y(2:end,:) - Y(1:end-1,:)*F')
%     R = cov(X - Y*H')
    Q = zeros(2);
    for t=2:length(X)
        tmp = Y(t,:)' - F*Y(t-1,:)';
        Q = Q + tmp * tmp';
    end
    Q = Q/(length(X)-1)

    R = zeros(2);
    for t=2:length(X)
        tmp = X(t,:)' - H*Y(t,:)';
        R = R + tmp * tmp';
    end
    R = R/(length(X)-1)
    
    K = 30; %# of iters for EM
    N = 200; %samples per iter
    
    p_data = zeros(K,1);
    for k=1:K
        fprintf('--------------iter #%d--------------\n',k);
        %E step
        Y_sample = repmat(Y,1,1,N);
        p_tmp = zeros(N,1);
        for n=1:N
            for t=2:size(Y,1)-1 %TBD: shuttle the order
%                 mu_ = ( Y(t+1) + Y(t-1) ) / 2;
%                 sigma_ = sqrt( Q^2 / 2 );
%                 mu = (X(t)*sigma_^2 + mu_*R^2) / (sigma_^2 + R^2);
%                 sigma = sqrt( Q^2/2 * R^2/(Q^2/2 + R^2) );
%                 Y_sample(n,t) = normrnd(mu, sigma); %sampling y_t
                Sigma_ = (Q^-1 + (F^-1*Q*(F^-1)')^-1)^-1;
                Mu_ = Sigma_* Q^-1 * (F*Y(t-1,:)') + Sigma_ * (F^-1*Q*(F^-1)')^-1 * (F^-1*Y(t+1,:)');
                Sigma = (Sigma_^-1 + (H^-1*R*(H^-1)')^-1)^-1;
                Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R*(H^-1)')^-1 * (H^-1*X(t,:)');
                Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
            end
            
            %data likelihood
            p = 0;
            for t=2:size(Y,1)-1 
                p = p + log( mvnpdf(Y_sample(t,:,n), Y_sample(t-1,:,n)*F', Q) )...
                    + log( mvnpdf(X(t,:), Y_sample(t,:,n)*H', R) );
            end
            p_tmp(n) = p;
        end

        p_data(k) = mean(p_tmp);
        
        %M step
        Y = mean(Y_sample,3);
%         Q = cov(Y(2:end,:) - Y(1:end-1,:)*F')
%         R = cov(X - Y*H')
        Q = zeros(2);
        for t=2:length(X)
            tmp = Y(t,:)' - F*Y(t-1,:)';
            Q = Q + tmp * tmp';
        end
        Q = Q/(length(X)-1)

        R = zeros(2);
        for t=2:length(X)
            tmp = X(t,:)' - H*Y(t,:)';
            R = R + tmp * tmp';
        end
        R = R/(length(X)-1)

        if debug==1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            pause(0.1)
        end
    end
    
    figure
    plot(p_data,'k','LineWidth',2)
    
    output = Y;