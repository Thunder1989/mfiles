function [Bmu,Bsigma] = data_mle(X, MODE)
 
    X = reshape(X(:), 96, []);
    Nd=7; Nh=size(X,1); 

    %------priors--------
    %X_B prior hyperparameter
    if MODE == 1 %96*7
        priors.mu_h = zeros(Nh,Nd)+50; %>40
        priors.sigma_h = zeros(Nh,Nd)+2; %>10
    else %96+7
        priors.mu_d = zeros(1,Nd)+50; %>40
        priors.sigma_d = zeros(1,Nd)+2; %>10
        priors.mu_h = zeros(Nh,1)+10; %>40
        priors.sigma_h = zeros(Nh,1)+2; %>10
    end
    %---------------------

    if MODE %96*7 parameters
        mu_h = zeros(Nh,Nd);
        sigma_h = zeros(Nh,Nd);
        for d=1:size(mu_h,2) %day
            for h=1:size(mu_h,1)   %interval
                hour_sample = X(h,d:7:end);
                [mu, sigma] = get_post_para(hour_sample, priors.mu_h(h,d), priors.sigma_h(h,d));
                mu_h(h,d) = mu;
                sigma_h(h,d) = sigma;
            end
        end
        Bmu = repmat(mu_h,1,size(X,2)/7);
        Bsigma = repmat(sigma_h,1,size(X,2)/7);

    else
        mu_d = zeros(1,Nd);
        sigma_d = zeros(1,Nd);
        mu_h = zeros(Nh,1);
        sigma_h = zeros(Nh,1);
        for d=1:size(mu_d,2)	%day
            day_sample = mean(X(:,d:7:end));
            [mu, sigma] = get_post_para(day_sample, priors.mu_d(d), priors.sigma_d(d));
            mu_d(d) = mu;
            sigma_d(d) = sigma;
        end
        residue = bsxfun(@minus, X, mean(X));
        for h=1:size(mu_h,1)	%interval
            hour_sample = residue(h,:);
            [mu, sigma] = get_post_para(hour_sample, priors.mu_h(h), priors.sigma_h(h));
            mu_h(h) = mu;
            sigma_h(h) = sigma;
        end
        d1 = repmat(mu_d, Nh, 1);
        d2 = repmat(mu_h, 1, Nd);
        d = d1 + d2;
        Bmu = repmat(d,1,size(X,2)/7);
        h1 = repmat(sigma_d, Nh, 1);
        h2 = repmat(sigma_h, 1, Nd);
        h = h1 + h2;
        Bsigma = repmat(h,1,size(X,2)/7);
        size(Bmu)
        size(Bsigma)

    end 
    
    
function [mu, sigma] = get_post_para(X, mu_0, sigma_0)
    ep = 10e-20;
    var_0 = sigma_0 ^ 2;
    var_n = var(X(:)) + ep;
    N = numel(X);
    var_ = ( 1/var_0 + N/var_n )^-1;
    mu = ( mu_0/var_0 + sum(X(:))/var_n ) * var_;
    sigma = sqrt(var_);
