function samples = sensorMCMC(X,ITERS,events,EQUIV)
% samples = sensorMCMC(Data, priors, [Niter Nburn Nplot], events, EQUIV)
%    Data   : (Ntimes x 7*Nweeks) matrix of count data (assumed starting Sunday)
%    Priors : structure with parameter values of prior distributions
%    [...]  : Iteration controls: total # of iterations, # used for burn-in,
%               and # of iterations between plots
%    Events : binary matrix with known events marked (used for plotting)
%    Sharing: used to force parameter sharing, = [S1,S2]
%              S1 = force sharing of delta (day effect) among days
%              S2 = force sharing of eta (time of day) among *ays
%              Values: 1 (all days are shared),  2 (weekdays/weekends),  3 (none)

if (nargin < 4), events=zeros(size(X)); end;
if (numel(events)==0), events = zeros(size(X)); end;

if (nargin < 5) EQUIV=[3,3]; end;

switch length(ITERS),
    case 0, Nplot=25;       Nburn=10;       Niter = 100;
    case 1, Nplot=25;       Nburn=10;       Niter = ITERS(1);
    case 2, Nplot=25;       Nburn=ITERS(2); Niter = ITERS(1);
    case 3, Nplot=ITERS(3); Nburn=ITERS(2); Niter = ITERS(1);
end;

fprintf('Running the Markov-modulated Gaussian model...\n');
fprintf('Data is %d days long with %d measurements per day;\n',size(X,2),size(X,1));
switch EQUIV(1) % d(t)
  case 1, fprintf('All days share total (per day) rate; ');
  case 2, fprintf('Weekend/weekdays share total (per day) rate; ');
  case 3, fprintf('Total (per day) rate unshared; ');
end
switch EQUIV(1) % h(t)
  case 1, fprintf('All days share time profile.\n');
  case 2, fprintf('Weekend/weekdays share time profile.\n');
  case 3, fprintf('Time profile unshared.\n');
end
fprintf('Running for %d iterations, with %d for burn-in and plotting every %d.\n',Niter,Nburn,Nplot);

Z=zeros(size(X)); X_B=max(X,1); 
M=[.99,.5;.01,.5]; %[p00, p10; p01, p11]
Nd=7; Nh=size(X,1); Nw=size(X,2)/7;
D=7*4; 

%------priors--------
priors.MODE = 0;

priors.sigma_B = 2;	%X_B sigma
priors.sigma_E = 3;	%X_E sigma

%X_B prior hyperparameter
if priors.MODE == 1 %96*7
    priors.mu_h = zeros(Nh,Nd)+50; %>40
    priors.sigma_h = zeros(Nh,Nd)+2; %>10
else %96+7
    priors.mu_d = zeros(1,Nd)+50; %>40
    priors.sigma_d = zeros(1,Nd)+2; %>10
    priors.mu_h = zeros(Nh,1)+10; %>40
    priors.sigma_h = zeros(Nh,1)+2; %>10
end

%X_E prior hyperparameter
priors.mu_0 = 1;
priors.sigma_0 = 2;

%transition prior hyperparameter
priors.z01 = .1*1000; priors.z00 = .9*1000;	% z(t) event process
priors.z10 = .25*1000; priors.z11 = .75*1000;

%---------------------

samples.Z = zeros([size(Z),Niter]);
samples.M  = zeros([size(M),Niter]);
samples.X_B = zeros([size(X_B),Niter]);
samples.P_data = zeros(1,Niter);
samples.Mu0 = zeros(1,Niter);
samples.Sigma0 = zeros(1,Niter);
samples.P_Z = zeros([2,numel(X),Niter]);
samples.P_E = zeros([2,numel(X),Niter]);
samples.prior = priors;

% MAIN LOOP: MCMC FOR INFERENCE
if priors.MODE
    Bmu = repmat(priors.mu_h,1,Nw);
    Bsigma = repmat(priors.sigma_h,1,Nw);
else
    d1 = repmat(priors.mu_d, Nh, 1);
    d2 = repmat(priors.mu_h, 1, Nd);
    d = d1 + d2;
    Bmu = repmat(d,1,size(X_B,2)/7);
    h1 = repmat(priors.sigma_d, Nh, 1);
    h2 = repmat(priors.sigma_h, 1, Nd);
    h = h1 + h2;
    Bsigma = repmat(h,1,size(X_B,2)/7);
end
Mu0 = priors.mu_0;
Sigma0 = priors.sigma_0;
A0 = log( M^100 * [1;0] );
for iter=1:Niter+Nburn
    fprintf('iter #%d\n',iter);
    [Z,X_B,P_data,P_Z,P_E,A0] = draw_Z_Para(X,Bmu,Bsigma,Mu0,Sigma0,M,priors,A0); %E step
    [Bmu,Bsigma,Mu0,Sigma0,priors] = draw_Para_SData(X,X_B,Z,Mu0,Sigma0,priors,EQUIV); %M step
    M = draw_M_Z(Z,priors);
    samples.P_data(iter) = P_data;
    samples.Mu0(iter) = Mu0;
    samples.Sigma0(iter) = Sigma0;
    samples.P_Z(:,:,iter) = P_Z;
    samples.P_E(:,:,iter) = P_E;
    samples.Z(:,:,iter-Nburn) = Z;
    samples.M(:,:,iter-Nburn) = M;
    samples.X_B(:,:,iter-Nburn) = X_B;    

    if (mod(iter,5)==0)
        figure
        D = 7*4; % # of days to plot
        hold on
        plot(X(1:4*24*D),'r','LineWidth',1.5)
        plot(Bmu(1:4*24*D),'k','LineWidth',1.5)
        plot(Bsigma(1:4*24*D),'b','LineWidth',1)
        tmp = mean(samples.X_B(:,:,iter-Nplot+1:iter),3);
        plot( tmp(1:4*24*D),'g')
        tmp = mean(samples.P_Z(:,:,iter-Nplot+1:iter),3);
        plot( tmp(2,1:4*24*D),'m')
        title('average results')
        legend('Original','Bmu','Bsigma','X\_B','P(Z)')
        day_bd = zeros(1, 4*24*D);
        day_bd(1:4*24:end) = 1*max(X(1:4*24*D));
        days = {'Fri','Sat','Sun','Mon','Tue','Wed','Thu'};
        ctr = 1;
        for t=1:numel(day_bd)
            if day_bd(t)==0
                continue
            end
            plot([t t], [-10 day_bd(t)], 'k--', 'LineWidth', 0.5);
            text(t+30,-10,days{mod(ctr-1,7)+1})
            ctr = ctr + 1;
        end
    end
end

samples.prior = priors;

% gen_movie(samples.P_Z,'event');
% gen_movie(samples.P_E,'emission');

% log p_data
figure
plot(samples.P_data,'k','LineWidth',2);

% E step
function [Z,X_B,P_data,P_Z,P_E,A0] = draw_Z_Para(X,Bmu,Bsigma,Mu0,Sigma0,M,prior,A_start)
    X_B=X; Z=0*X; ep=1e-50;

    a0 = A_start; %starting point for alpha
    M = log(M);
    pe=zeros(2,numel(X)); a=zeros(2,numel(X));

    %emission probability
    for t=1:numel(X)
        Bmu_t = Bmu(t);
        Bsigma_t = sqrt( Bsigma(t)^2 + prior.sigma_B^2 );
        pe(1,t) = log( normpdf(X(t), Bmu_t, Bsigma_t) ); %z_t=0
        
        sigma1 = Bsigma_t;
        sigma2 = sqrt( prior.sigma_E^2 + Sigma0^2 );
        pe(2,t) = log( 1/sqrt(2*pi*(sigma1^2+sigma2^2)) ) ...
            + -(X(t)-Bmu_t-Mu0)^2 / ( 2*(sigma1^2+sigma2^2) );
    end
    
    % forward
    a(:,1) = a0 + pe(:,1); 
    % for t=2,...N, we compute \alpha(t,k) = pe(t,k) * \sum_{j} M_jk * \alpha(t-1,j) 
    for t=2:numel(X)
        a(1,t) = logsumexp( M(1,:) + a(:,t-1)' ) + pe(1,t); 
        a(2,t) = logsumexp( M(2,:) + a(:,t-1)' ) + pe(2,t); 
    end
    P_data = logsumexp( a(:,end) );

    % backward
    b=zeros(2,numel(X));
    b(:,end) = [0;0]; 
    % for t=N-1,...1, we compute \beta(t,k) = \sum_{j} M_kj * \beta(t+1,j) * pe(t+1,j)
    for t=numel(X)-1:-1:1
        b(1,t) = logsumexp( M(:,1) + b(:,t+1) + pe(:,t+1) ); 
        b(2,t) = logsumexp( M(:,2) + b(:,t+1) + pe(:,t+1) ); 
    end
    
    % merge and normalize
    p = a + b;
%     norm = logsumexp(p);
%     assert ( isequal(bsxfun(@minus, p, norm), bsxfun(@minus, p, p_data) ));
    p = bsxfun(@minus, p, P_data); %normalization
    A0 = p(:,1);

    % sampling
    for t=numel(X):-1:1
        if ( log(rand(1)) > p(1,t) )	% if event at time t
            Z(t)=1;
            Bmu_t = Bmu(t);
            Bsigma_t = sqrt( Bsigma(t)^2 + prior.sigma_B^2);
            sigma1 = Bsigma_t; sigma2 = sqrt( prior.sigma_E^2 + Sigma0^2);
            mu12 = ( (X(t)-Bmu_t)*sigma2^2 + Mu0*sigma1^2 ) / (sigma1^2 + sigma2^2);
            sigma12 = sqrt( sigma1^2*sigma2^2 / (sigma1^2+sigma2^2) );
            X_B(t) = normrnd(mu12, sigma12); % sampling X_B
        else
            Z(t)=0; X_B(t)=X(t); % no event
        end
    end
    P_Z = p;
    P_E = pe;

% Given Z, Sample M
function [M] = draw_M_Z(Z,prior)
    n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0)); %checked
    n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
    z0 = betarnd(n01+prior.z01, n0-n01+prior.z00);
    z1 = betarnd(n10+prior.z10, n1-n10+prior.z11);
    M = [1-z0, z1; z0, 1-z1];

% M step
function [Bmu,Bsigma,Mu0,Sigma0,prior] = draw_Para_SData(X,X_B,Z,Mu0_,Sigma0_,prior,EQUIV)
    Nd=7; Nh=size(X_B,1); Nw=size(X_B,2)/7;
    X_E = X - X_B;
    
    %compute sigma_E and posterior hyperparameters for mu_E, only if X_E exists
    if ~isempty( find(Z~=0,1) )
        data = X_E(Z==1);
        [mu, sigma] = get_post_para(data, prior.mu_0, prior.sigma_0);
        Mu0 = mu;
        Sigma0 = sigma;
        prior.sigma_E = sqrt(var(data(:)));
    else
        Mu0 = Mu0_;
        Sigma0 = Sigma0_;
    end
    
    %compute posterior hyperparameters for mu_B(t)
    if prior.MODE %96*7 parameters
        mu_h = zeros(Nh,Nd);
        sigma_h = zeros(Nh,Nd);
        for d=1:size(mu_h,2) %day
            for h=1:size(mu_h,1)   %interval
                hour_sample = X_B(h,d:7:end);
                [mu, sigma] = get_post_para(hour_sample, prior.mu_h(h,d), prior.sigma_h(h,d));
                mu_h(h,d) = mu;
                sigma_h(h,d) = sigma;
            end
        end
        Bmu = repmat(mu_h,1,size(X_B,2)/7);
        Bsigma = repmat(sigma_h,1,size(X_B,2)/7);

    else
        mu_d = zeros(1,Nd);
        sigma_d = zeros(1,Nd);
        mu_h = zeros(Nh,1);
        sigma_h = zeros(Nh,1);
        for d=1:size(mu_d,2)	%day
            day_sample = mean(X_B(:,d:7:end));
            [mu, sigma] = get_post_para(day_sample, prior.mu_d(d), prior.sigma_d(d));
            mu_d(d) = mu;
            sigma_d(d) = sigma;
        end
        residue = bsxfun(@minus, X, mean(X));
        for h=1:size(mu_h,1)	%interval
            hour_sample = residue(h,:);
            [mu, sigma] = get_post_para(hour_sample, prior.mu_h(h), prior.sigma_h(h));
            mu_h(h) = mu;
            sigma_h(h) = sigma;
        end
        d1 = repmat(mu_d, Nh, 1);
        d2 = repmat(mu_h, 1, Nd);
        d = d1 + d2;
        Bmu = repmat(d,1,size(X_B,2)/7);
        h1 = repmat(sigma_d, Nh, 1);
        h2 = repmat(sigma_h, 1, Nd);
        h = h1 + h2;
        Bsigma = repmat(h,1,size(X_B,2)/7);
        size(Bmu)
        size(Bsigma)
    end

    prior.sigma_B = sqrt(var(X_B(:)));
    prior
    
    %TBD: enforce paramter sharing between days
    switch EQUIV(1)
        case 1,
        case 2,
        case 3,
    end
    
    switch EQUIV(2)
        case 1,
        case 2,
        case 3,
    end    
    
    %debugging block
%     figure
%     D = 7*4; % # of days to plot
%     hold on
%     plot(X(1:4*24*D),'r','LineWidth',1.5)
%     plot(Bmu(1:4*24*D),'k','LineWidth',1.5)
%     plot(X(1:4*24*D) - Bmu(1:4*24*D),'b','LineWidth',1.5)
%     plot(Bsigma(1:4*24*D),'g','LineWidth',1)
%     legend('Original','Bmu','Event','Bsigma')
%     pause
%     day_bd = zeros(1, numel(N));
%     day_bd(1:4*24:end) = 1*max(N(:));
%     for t=1:numel(day_bd)
%         if day_bd(t)==0
%             continue
%         end
%         plot([t t], [-10 day_bd(t)], 'k--', 'LineWidth', 0.5);
%     end
%     Bmu_test = repmat(Dmu,Nh,1) + Hmu;
%     Bmu_test = repmat(Bmu_test,1,4);
%     Bsigma_test = sqrt( repmat(Dsigma,Nh,1).^2 + Hsigma.^2 );
%     Bsigma_test = repmat(Bsigma_test,1,4);
%     assert ( isequal(Bmu,Bmu_test) );
%     assert ( isequal(Bsigma,Bsigma_test) );

%compute posterior hyperparameters
function [mu, sigma] = get_post_para(X, mu_0, sigma_0)
    ep = 10e-20;
    var_0 = sigma_0 ^ 2;
    var_n = var(X(:)) + ep;
    N = numel(X);
    var_ = ( 1/var_0 + N/var_n )^-1;
    mu = ( mu_0/var_0 + sum(X(:))/var_n ) * var_;
    sigma = sqrt(var_);

function s = logsumexp(x, dim)
    if nargin == 1 
        dim = find(size(x)~=1,1);
        if isempty(dim), dim = 1; end
    end

    % subtract the largest in each column
    y = max(x,[],dim);
    x = bsxfun(@minus,x,y);
    % assert( isempty(find(x<-750)) )
    s = y + log(sum(exp(x),dim));
    i = find(~isfinite(y));
    if ~isempty(i)
        s(i) = y(i);
    end

function gen_movie(y, name)
     % Set up the movie.
    fn = sprintf('%s.avi',name);
    writerObj = VideoWriter(fn); 
    writerObj.FrameRate = 1; % frames per second.
    open(writerObj); 

    D = 7*2;
    % fid = figure;
    for i=1:size(y,3)
    %     pause(0.1);
        figure
        hold on
        y_tmp = y(:,:,i);
        if size(y_tmp,1)>1
            plot(y_tmp(1,1:4*24*D),'r','LineWidth',1.5);
            plot(y_tmp(2,1:4*24*D),'k','LineWidth',1.5);
        else
            plot(y_tmp(1:4*24*D),'r','LineWidth',1.5);
        end
        title(sprintf('iter %d',i))
    %     if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, frame);
    %     end

    end
    close(writerObj); % Saves the movie.
