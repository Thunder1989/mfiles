function samples = sensorMCMC(X,priors,ITERS,events,EQUIV)
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
end;
switch EQUIV(1) % h(t)
  case 1, fprintf('All days share time profile.\n');
  case 2, fprintf('Weekend/weekdays share time profile.\n');
  case 3, fprintf('Time profile unshared.\n');
end;
fprintf('Running for %d iterations, with %d for burn-in and plotting every %d.\n',Niter,Nburn,Nplot);

Z=zeros(size(X)); X_B=max(X,1); 
M=[.99,.5;.01,.5]; %[p00, p10; p01, p11]
Nd=7; Nh=size(X,1); Nw=size(X,2)/7;
samples.Z = zeros([size(Z),Niter]);
samples.M  = zeros([size(M),Niter]);
samples.X_B = zeros([size(X_B),Niter]);
samples.logp_NgLM = zeros(1,Niter); 
samples.logp_NgLZ = zeros(1,Niter);
samples.P_data = zeros(1,Niter);
samples.P_Z = zeros([2,numel(X),Niter]);
samples.P_E = zeros([2,numel(X),Niter]);

% MAIN LOOP: MCMC FOR INFERENCE
Bmu = repmat(priors.Hmu,1,Nw);
Bsigma = repmat(priors.Hsigma,1,Nw);
Mu0 = priors.mu0;
Sigma0 = priors.sigma0;
A0 = log( M^100 * [1;0] );
for iter=1:Niter+Nburn
    iter
    [Z,X_B,P_data,P_Z,P_E,A0] = draw_Z_Para(X,Bmu,Bsigma,Mu0,Sigma0,M,priors,A0); %E step
    [Bmu,Bsigma,Mu0,Sigma0] = draw_Para_SData(X,X_B,Mu0,Sigma0,priors,EQUIV); %M step
    M = draw_M_Z(Z,priors);
    samples.P_data(iter) = P_data;
    samples.P_Z(:,:,iter) = P_Z;
    samples.P_E(:,:,iter) = P_E;
    
  if (iter > Nburn)
    samples.Z(:,:,iter-Nburn) = Z;
    samples.M(:,:,iter-Nburn) = M;
    samples.X_B(:,:,iter-Nburn) = X_B;
  end
end

gen_movie(samples.P_Z,'event');
gen_movie(samples.P_E,'emission');

figure
plot(samples.P_data,'k','LineWidth',2);

%% E step
function [Z,X_B,P_data,P_Z,P_E,A0] = draw_Z_Para(X,Bmu,Bsigma,Mu0,Sigma0,M,prior,A_start)
    X_B=X; Z=0*X; ep=1e-50;

    a0 = A_start; %starting point for alpha
    M = log(M);
    pe=zeros(2,numel(X)); a=zeros(2,numel(X));

    %emission probability
    for t=1:numel(X)
        Bmu_t = Bmu(t);
        Bsigma_t = sqrt( Bsigma(t)^2 + prior.Bsigma^2 );
        pe(1,t) = log( normpdf(X(t), Bmu_t, Bsigma_t) ); %z_t=0
        
        sigma1 = Bsigma_t;
        sigma2 = sqrt( prior.Esigma^2 + Sigma0^2 );
        sigma12 = sqrt( sigma1^2*sigma2^2 / (sigma1^2+sigma2^2) );
        pe(2,t) = log( 1/(2*pi*sigma1*sigma2) * sqrt(2*pi*sigma12^2) ) ...
            + -(X(t)-Bmu_t-Mu0)^2 / 2*(sigma1^2+sigma2^2);
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
        if ( log(rand(1)) > p(1,t))	% if event at time t
            Z(t)=1;
            Bmu_t = Bmu(t);
            Bsigma_t = sqrt( Bsigma(t)^2 + prior.Bsigma^2);
            sigma1 = Bsigma_t; sigma2 = sqrt( prior.Esigma^2 + Sigma0^2);
            mu12 = ( (X(t)-Bmu_t)*sigma2^2 + Mu0*sigma1^2 ) / (sigma1^2 + sigma2^2);
            sigma12 = sigma1^2*sigma2^2 / (sigma1^2+sigma2^2);
            X_B(t) = normrnd(mu12, sigma12); % sampling X_B
        else
            Z(t)=0; X_B(t)=X(t); % no event
        end
    end
    P_Z = p;
    P_E = pe;

%% GIVEN Z, SAMPLE M
function [M] = draw_M_Z(Z,prior)
    n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0));
    n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
    z0 = betarnd(n01+prior.z01, n0-n01+prior.z00);
    z1 = betarnd(n10+prior.z10, n1-n10+prior.z11);
    M = [1-z0, z1; z0, 1-z1];

%% M step
function [Bmu,Bsigma,Mu0,Sigma0] = draw_Para_SData(X,X_B,Mu0_,Sigma0_,prior,EQUIV)
    Nd=7;	Nh=size(X_B,1);
    X_E = X - X_B;
    
    %compute posterior hyperparameters for mu_E
%     assert ( ~isempty(find(NE~=0,1)))
    if ~isempty( find(X_E~=0,1) )
        [mu, sigma] = get_post_para(X_E, prior.mu0, prior.sigma0);
        Mu0 = mu;
        Sigma0 = sigma;
    else
        Mu0 = Mu0_;
        Sigma0 = Sigma0_;
    end
    
    %compute posterior hyperparameters for mu_B(t)
    Hmu = zeros(Nh,Nd);
    Hsigma = zeros(Nh,Nd);
    for d=1:size(Hmu,2) %day
        for h=1:size(Hmu,1)   %interval
            hour_sample = X_B(h,d:7:end);
            [mu, sigma] = get_post_para(hour_sample, prior.Hmu(h,d), prior.Hsigma(h,d));
            Hmu(h,d) = mu;
            Hsigma(h,d) = sigma;
        end
    end
    Bmu = repmat(Hmu,1,size(X_B,2)/7);
    Bsigma = repmat(Hsigma,1,size(X_B,2)/7);

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
    
    %test block
%     figure
%     hold on
%     plot(reshape(X,1,[]),'r','LineWidth',1.5)
%     plot(reshape(Bmu,1,[]),'k','LineWidth',1.5)
%     plot(reshape(Bsigma,1,[]),'b','LineWidth',1)
%     pause
%     plot(reshape(X,1,[]) - reshape(Bmu,1,[]),'k','LineWidth',1.5)
%     legend('Original','Baseline','Event')
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

W = 7*2;
% fid = figure;
for i=1:size(y,3)
%     pause(0.1);
    figure
    hold on
    y_tmp = y(:,:,i);
    plot(y_tmp(1,1:4*24*W),'r','LineWidth',1.5);
    plot(y_tmp(2,1:4*24*W),'k','LineWidth',1.5);
    title(sprintf('iter %d',i))
%     if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
%     end
 
end
close(writerObj); % Saves the movie.
