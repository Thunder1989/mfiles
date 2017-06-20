function output = gibbs_hmm(data, debug)
    %initialization
    X1 = data(:)';
    X2 = [0 diff(X1)];
    Y1 = EWMA(X1,10);
    Y2 = EWMA(X2,10);
    F = [1 1; 0 1/4];
    H = [1 0];
    
    fprintf('initial Q and R:\n')
    Q = cov(diff(Y1),diff(Y2))
    R = cov(X1 - Y1, X2 - Y2)
    
    X = [X1(:) X2(:)];
    Y = [Y1(:) Y2(:)];
    
    K = 20; %# of iters for EM
    N = 200; %samples per iter
    
    p_data = zeros(K,1);
    for k=1:K
        fprintf('iter #%d...\n',k);
        %E step
        Y_sample = repmat(Y,1,1,N);
        p_tmp = zeros(N,1);
        for n=1:N
            for t=2:size(Y,2)-1 %TBD: shuttle the order
%                 mu_ = ( Y(t+1) + Y(t-1) ) / 2;
%                 sigma_ = sqrt( Q^2 / 2 );
%                 mu = (X(t)*sigma_^2 + mu_*R^2) / (sigma_^2 + R^2);
%                 sigma = sqrt( Q^2/2 * R^2/(Q^2/2 + R^2) );
%                 Y_sample(n,t) = normrnd(mu, sigma); %sampling y_t
                Sigma_ = (Q^-1 + Q^-1)^-1;
                Mu_ = Sigma_*Q^-1*(F*Y(t-1,:)') + Sigma_*Q^-1*(F*Y(t+1,:)');
                Sigma = (Sigma_^-1 + R^-1)^-1;
                Mu = Sigma*Sigma_^-1*Mu_ + Sigma*R^-1*(H*X(t,:)');
                Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
            end
            
%             %data likelihood
%             p = 0;
%             for t=2:size(Y,2)-1
%                 p = p + log( normpdf(Y_sample(n,t), Y_sample(n,t-1), Q) )...
%                     + log( normpdf(X(t), Y_sample(n,t), R) );
%             end
%             p_tmp(n) = p;
        end
        
%         p_data(k) = mean(p_tmp);
        
        %M step
        Y = mean(Y_sample,3);
        Q = cov(diff(Y))
        R = cov(X - Y)
        
        if debug == 1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            pause(0.1)
        end
    end
    
%     figure
%     plot(p_data, 'k','LineWidth',2)
    
    output = Y;