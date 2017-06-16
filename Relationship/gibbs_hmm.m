function output = gibbs_hmm(data, debug)
    %initialization
    X = data(:)';
    Y = EWMA(X,10);
    Q = std(Y);
    R = std(X - Y);
    
    K = 20; %# of iters for EM
    N = 200; %samples per iter
    
    for k=1:K
        fprintf('iter #%d...\n',k);
        %E step
        Y_sample = repmat(Y,N,1);
        for n=1:N
            for t=2:size(Y,2)-1
                mu_ = ( Y(t) + Y(t-1) ) / 2;
                sigma_ = sqrt(Q^2 / 2);
                mu = (X(t)*sigma_^2 + mu_*R^2) / (sigma_^2 + R^2);
                sigma = sqrt( Q^2/2 * R^2/(Q^2/2 + R^2) );
                Y_sample(n,t) = normrnd(mu, sigma);%sampling here
            end
        end
        
        %M step
        Y = mean(Y_sample);
%         Y_ = EWMA(Y,10);
        Q = std(Y);
        R = std(X - Y);
        
        if debug == 1
            figure
            hold on
            plot(X,'k','LineWidth',2)
            plot(Y,'r--','LineWidth',1)
            pause(0.1)
        end
    end
    
    output = Y;