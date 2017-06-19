function output = gibbs_hmm(data, debug)
    %initialization
    X = data(:)';
    X = [0 diff(X)];
    Y = EWMA(X,10);
%     Y = X + randn(size(X));
    fprintf('initial Q and R:\n')
    Q = std(diff(Y));
    R = std(X - Y);
    
    K = 10; %# of iters for EM
    N = 200; %samples per iter
    
    p_data = zeros(K,1);
    for k=1:K
        fprintf('iter #%d...\n',k);
        %E step
        Y_sample = repmat(Y,N,1);
        p_tmp = zeros(N,1);
        for n=1:N
            for t=2:size(Y,2)-1 %TBD: shuttle the order
                mu_ = ( Y(t+1) + Y(t-1) ) / 2;
                sigma_ = sqrt( Q^2 / 2 );
                mu = (X(t)*sigma_^2 + mu_*R^2) / (sigma_^2 + R^2);
                sigma = sqrt( Q^2/2 * R^2/(Q^2/2 + R^2) );
                Y_sample(n,t) = normrnd(mu, sigma); %sampling y_t
            end
            
            %data likelihood
            p = 0;
            for t=2:size(Y,2)-1
                p = p + log( normpdf(Y_sample(n,t), Y_sample(n,t-1), Q) )...
                    + log( normpdf(X(t), Y_sample(n,t), R) );
            end
            p_tmp(n) = p;
        end
        
        p_data(k) = mean(p_tmp);
        
        %M step
        Y = mean(Y_sample);
        Q = std(diff(Y));
        R = std(X - Y);   
        
        if debug == 1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            pause(0.1)
        end
    end
    
    figure
    plot(p_data, 'k','LineWidth',2)
    
    output = Y;