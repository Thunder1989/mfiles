function samples = sensorMCMC(N,priors,ITERS,events,EQUIV)
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

if (nargin < 4), events=zeros(size(N)); end;
if (numel(events)==0), events = zeros(size(N)); end;

if (nargin < 5) EQUIV=[3,3]; end;

switch length(ITERS),
    case 0, Nplot=25;       Nburn=10;       Niter = 100;
    case 1, Nplot=25;       Nburn=10;       Niter = ITERS(1);
    case 2, Nplot=25;       Nburn=ITERS(2); Niter = ITERS(1);
    case 3, Nplot=ITERS(3); Nburn=ITERS(2); Niter = ITERS(1);
end;

fprintf('Running the Markov-modulated Poisson model...\n');
fprintf('Data is %d days long with %d measurements per day;\n',size(N,2),size(N,1));
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

% N=round(N);
% N=N-min(min(N));
% N = raw2count(N,1);
% events = zeros(size(N));
Z=zeros(size(N)); N0=max(N,1); NE=zeros(size(N)); L=(N+5)/2; %L is not used
M=[.99,.5;.01,.5]; 
% M=[0.75 0.25; 0.25 0.75];
% xs = 0:80;
Nd=7; Nh=size(N,1);
samples.L = zeros([size(L),Niter]);
samples.Z = zeros([size(Z),Niter]);
samples.PE = zeros([2,numel(N),Niter]);
samples.M  = zeros([size(M),Niter]);
samples.N0 = zeros([size(N0),Niter]);
samples.NE = zeros([size(NE),Niter]);
samples.logp_NgLM = zeros(1,Niter); 
samples.logp_NgLZ = zeros(1,Niter);

% MAIN LOOP: MCMC FOR INFERENCE
Bmu = repmat(priors.Dmu,Nh,1) + priors.Hmu;
Bmu = repmat(Bmu,1,4);
Bsigma = sqrt( repmat(priors.Dsigma,Nh,1).^2 + priors.Hsigma.^2 );
Bsigma = repmat(Bsigma,1,4);
Mu0 = priors.mu0;
Sigma0 = priors.sigma0;
for iter=1:Niter+Nburn,
    [Z,N0,NE,PE] = draw_Z_NLM(N,Bmu,Bsigma,Mu0,Sigma0,M,priors);
    [Bmu,Bsigma,Mu0,Sigma0] = draw_L_N0(N,N0,priors,EQUIV);
    M = draw_M_Z(Z,priors);
  
  if (iter > Nburn)
    samples.Z(:,:,iter-Nburn) = Z;
    samples.PE(:,:,iter-Nburn) = PE;    samples.M(:,:,iter-Nburn) = M;
    samples.N0(:,:,iter-Nburn) = N0; samples.NE(:,:,iter-Nburn) = NE;
    %samples.logp_NgLM(iter-Nburn) = eval_N_LM(N,L,M,priors);
    %samples.logp_NgLZ(iter-Nburn) = eval_N_LZ(N,L,Z,priors);
    %[logpC, logpGD, logpGDz] = logp(N,samples,priors,iter-Nburn,EQUIV);  
    %logpC=logpC/log(2); logpGD=logpGD/log(2); logpGDz=logpGDz/log(2); 
%     fprintf('\n Iter %d Est Marginal Likelihd: ln P(Data) = %.1f  (%.3f per time)\n',iter, logpC,logpC/numel(N));
%     mmppPlot(L,Z,N,NE,events,123); title(sprintf('iter %d',iter)); %pause(.5);
  end;
%   fprintf('.');         % DISPLAY / PLOT CURRENT SAMPLES & AVERAGES

%   if (mod(iter,Nplot)==0 && iter > Nburn)  
%     mmppPlot(mean(samples.L(:,:,1:iter-Nburn),3), ...
%       mean(samples.Z(:,:,1:iter-Nburn),3), N, events,101); figure(101); title('Posterior Averages');      
%     pause(.1);
%   end;
end;
% [logpC, logpGD, logpGDz] = logp(N,samples,priors,iter-Nburn,EQUIV); 
% logpC=logpC/log(2); logpGD=logpGD/log(2); logpGDz=logpGDz/log(2); 
% fprintf('\n Final Est Marginal Likelihd: ln P(Data) = %.1f  (%.3f per time)\n',logpC,logpC/numel(N));
%samples.logpC = logpC;
%samples.logpGD = logpGD;
%samples.logpGDz = logpGDz;

%Function for converting raw data to count data 
function count = raw2count(N, base)
hr = 2;
sample_freq = 4;
mean_day = mean(N);
Nd = size(N,2);
mean_day = repmat(mean_day, size(N,1)/(hr*sample_freq), 1);
mean_day = reshape(mean_day, 1, []); %comparison base
switch base
    case 1,
        N = reshape(N, sample_freq * hr, []);
        count = zeros(size(mean_day));
        for i=1:size(N,2)
            count(i) = sum( N(:,i)>=mean_day(i) );
        end
        count = reshape(count, [], Nd);
    case 2, 
end

%% EVALUATION FUNCTIONS
function [logpC, logpGD, logpGDz] = logp(N,samples,priors,iter,EQUIV);
% estimate the marginal likelihood of the data using our samples
% (Produces three estimates; Chib's is probably the best of the three)

  tmp = samples.logp_NgLZ(1:iter); tmpm=mean(tmp); tmp=tmp-tmpm;
  logpGDz = log(1./mean(1./exp(tmp)))+tmpm;    % Gelfand-Dey estimate

  tmp = samples.logp_NgLM(1:iter); tmpm=mean(tmp); tmp=tmp-tmpm;
  logpGD = log(1./mean(1./exp(tmp)))+tmpm;     % Gelfand-Dey estimate, marginalizing over Z

  Lstar = mean(samples.L(:,:,1:iter),3); Mstar = mean(samples.M(:,:,1:iter),3); 
  logp_LMgN=zeros(1,iter);
  logp_LM = eval_L_N0(Lstar,[],priors,EQUIV) + eval_M_Z(Mstar,[],priors);
  logp_NgLM= eval_N_LM(N,Lstar,Mstar,priors);
  for ii=1:iter
    logp_LMgN(ii) = eval_L_N0(Lstar,samples.N0(:,:,ii),priors,EQUIV)...
                  + eval_M_Z(Mstar,samples.Z(:,:,ii),priors);
  end;
  tmpm=mean(logp_LMgN); logp_LMgN=logp_LMgN-tmpm; logp_LMgN = log(mean(exp(logp_LMgN)))+tmpm;
  logpC = logp_NgLM + logp_LM - logp_LMgN;       % Chib estimate
  
function logp = eval_M_Z(M,Z,prior)		% evaluate p(M|Z)
  z1 = M(1,2); z0 = M(2,1);
  if (~isempty(Z))
    n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0));
    n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
  else n01=0; n0=0; n10=0; n1=0; 
  end;
  logp = log(betapdf(z0,n01+prior.z01,n0-n01+prior.z00)) + log(betapdf(z1,n10+prior.z10,n1-n10+prior.z11));

function logp = eval_L_N0(L,N0,prior,EQUIV)	% evaluate p(L | N0)
  L0 = mean(mean(L)); Nd = 7; Nh=size(L,1);
  for i=1:Nd, D(i) = mean(L(:,i)/L0); end;
  for i=1:Nd, for j=1:Nh, A(j,i) = L(j,i)/L0/D(i); end; end;
  logp = 0;

  % ENFORCE PARAMETER SHARING
  paD=prior.aD; aD=zeros(1,Nd); paH=prior.aH; aH=zeros(Nh,Nd);
  if (~isempty(N0))
    for i=1:Nd, aD(i) = sum(sum(N0(:,i:Nd:end))); end;
    for i=1:Nd, for j=1:Nh, aH(j,i) = sum(sum(N0(j,i:Nd:end))); end; end;
  end;
  switch EQUIV(1) % d(t)
    case 1, D = sum(D); paD = sum(paD); aD=sum(aD); 
    case 2, D = [D(1)+D(7),sum(D(2:6))]; paD=[paD(1)+paD(7),sum(paD(2:6))]; aD=[aD(1)+aD(7),sum(aD(2:6))];
    case 3, D = D; paD = paD; paH=paH;
  end;
  switch EQUIV(2) % tau(t)
    case 1, A=sum(A,2)/Nd; aH=sum(aH,2); paH=sum(paH,2);
    case 2, A = [(A(:,1)+A(:,7))/2,sum(A(:,2:6),2)/5]; aH = [aH(:,1)+aH(:,7),sum(aH(:,2:6),2)]; paH = [paH(:,1)+paH(:,7),sum(paH(:,2:6),2)]; 
    case 3, A=A; aH=aH; paH=paH;
  end;
  
  logp = logp + log(gampdf(L0,sum(sum(N0))+prior.aL,1/(numel(N0)+prior.bL)));
  logp = logp + dirlnpdf(D/Nd,aD + paD);
  for i=1:size(A,2),
    logp = logp + dirlnpdf(A(:,i)/Nh,aH(:,i)+paH(:,i));
  end;
  
function logp = eval_N_LZ(N,L,Z,prior)		% evaluate p(N|L,Z)
  logp = 0;
  for t=1:numel(N)
    if (N(t)~=-1),
      if (Z(t)==0), logp = logp + log(exppdf(N(t),L(t)));
      else          logp = logp + log(sum( exppdf(0:N(t),L(t)) .* nbinpdf(N(t):-1:0,prior.aE,prior.bE/(1+prior.bE)) ));
      end;
  end; end;

function logp = eval_N_LM(N,L,M,prior)		% evaluate p(N | L,M)
  PRIOR = M^100 * [1;0]; po=zeros(2,numel(N)); p=zeros(2,numel(N));
  for t=1:numel(N),
    if (N(t)~=-1)
      po(1,t) = exppdf(N(t),L(t));      
      po(2,t) = sum( exppdf(0:N(t),L(t)) .* nbinpdf(N(t):-1:0,prior.aE,prior.bE/(1+prior.bE)) );
    else po(1,t)=1; po(2,t)=1;
    end;
  end;
  p(:,1) = PRIOR .* po(:,1); sp=sum(p(:,1)); 
  logp = log(sp); p(:,1)=p(:,1)/sp;
  for t=2:numel(N), 
    p(:,t) = (M*p(:,t-1)).*po(:,t); sp=sum(p(:,t));
    logp = logp + log(sp); p(:,t)=p(:,t)/sp;
  end;

%% MISC PDF FUNCTIONS
function p = dirpdf(X,A)			% evaluate a dirichlet distribution
  k = length(X); if (k==1) p=1; return; end;
  logp = sum( (A-1).*log(X) ) - sum(gammaln(A)) + gammaln(sum(A));
  p = exp(logp);  
function logp = dirlnpdf(X,A)			% eval log(dirichlet)
  k = length(X); if (k==1) p=1; return; end;
  logp = sum( (A-1).*log(X) ) - sum(gammaln(A)) + gammaln(sum(A));  
function p = poisspdf(X,L)			% poisson distribution, and is faster than the one provided by matlab
  lnp = -L -gammaln(X+1) +log(L).*X;
  p = exp(lnp);
function lnp = poisslnpdf(X,L)			% log(poisson)
  lnp = -L -gammaln(X+1) +log(L).*X;
function p = nbinpdf(X,R,P)			% negative binomial distribution
  lnp = gammaln(X+R)-gammaln(R)-gammaln(X+1)+log(P).*R+log(1-P).*X;
  p = exp(lnp);
function lnp = nbinlnpdf(X,R,P)			% log(neg binomial)
  lnp = gammaln(X+R)-gammaln(R)-gammaln(X+1)+log(P).*R+log(1-P).*X;
function p = lomaxpdf(X,L,A) %L is scale, A is shape
    p = A/L * (1+X/L) .^ (-A-1);
function lnp = lomaxlnpdf(X,L,A) %L is scale, A is shape
    p = A/L * (1+X/L) .^ (-A-1);
    lnp = log(p);
function lnp = explnpdf(X,L)
    lnp = log(exppdf(X,L));

%% SAMPLING FUNCTIONS
function [Z,N0,NE,PE] = draw_Z_NLM(N,Bmu,Bsigma,Mu0,Sigma0,M,prior)
    N0=N; NE=0*N; Z=0*N; ep=1e-50;

    PRIOR = log( M^100 * [1;0] ); 
    M = log(M);
    pe=zeros(2,numel(N)); p=zeros(2,numel(N));
    %emission probability
    for t=1:numel(N)
        Bmu_t = Bmu(t);
        Bsigma_t = sqrt( Bsigma(t)^2 + prior.Bsigma^2 );
        pe(1,t) = log( normpdf(N(t), Bmu_t, Bsigma_t) ); %z_t=0
%         figure(333);hist(exppdf(0:N(t),L(t)),10,'replace');
%         figure(333);hist(nbinpdf(N(t):-1:0,prior.aE,prior.bE/(1+prior.bE)),20,'replace');
        
        sigma1 = Bsigma_t;
        sigma2 = sqrt( prior.Esigma^2 + Sigma0^2 );
        mu12 = ( (N(t)-Bmu_t)*sigma2^2 + Mu0*sigma1^2 ) / (sigma1^2 + sigma2^2);
        sigma12 = sqrt( sigma1^2*sigma2^2/ (sigma1^2+sigma2^2) );
%         pe(2,t) = 1/(2*pi*sigma1*sigma2) * sqrt(2*pi*sigma12^2) ...
%             * exp( -(N(t)-Bmu_t-prior.mu0)^2 / 2*(sigma1^2+sigma2^2) ) ...
%             * normcdf(N(t), mu12, sigma12);
        pe(2,t) = log( 1/(2*pi*sigma1*sigma2) * sqrt(2*pi*sigma12^2) ) ...
            + -(N(t)-Bmu_t-prior.mu0)^2 / 2*(sigma1^2+sigma2^2) ...
            + log( normcdf(N(t), mu12, sigma12) );
    end
%     figure
%     hold on
%     plot(pe(1,:),'r')
%     plot(pe(2,:),'k')
    
    % Compute forward posterior marginals
    p(:,1) = PRIOR + pe(:,1); 
    %p(:,1) = p(:,1)/sum(p(:,1));
    for t=2:numel(N)
        p(1,t) = logsumexp( M(1,:) + p(:,t-1)' ) + pe(1,t); 
        p(2,t) = logsumexp( M(2,:) + p(:,t-1)' ) + pe(2,t); 
    end
%     assert ( ~isempty(find(p(2,:)==1)) )

    % backward
    b=zeros(2,numel(N));
    b(:,end) = [0;0]; 
    for t=numel(N)-1:-1:1
        b(1,t) = logsumexp( M(:,1) + b(:,t+1) + pe(:,t+1) ); 
%         p(1,t) = p(1,t) + b(1,t);
        b(2,t) = logsumexp( M(:,2) + b(:,t+1) + pe(:,t+1) ); 
%         p(2,t) = p(2,t) + b(2,t);
    end
    figure
    hold on
    plot(b(1,:),'r')
    plot(b(2,:),'k')
    pause    
    
    % TBD: sampling
    for t=numel(N):-1:1
        if (rand(1) > p(1,t)),  % if event at time t
            if (N(t)~=-1)
                Z(t)=1;
                Bsigma_t = sqrt( Bsigma(t)^2 + prior.Bsigma^2);
                sigma1 = Bsigma_t; sigma2 = sqrt( prior.Esigma^2 + Sigma0^2);
                mu12 = ( (N(t)-Bmu_t)*sigma2^2 + Mu0*sigma1^2 ) / (sigma1^2 + sigma2^2);
                sigma12 = sigma1^2*sigma2^2/ (sigma1^2+sigma2^2);
                N0(t) = normrnd(mu12, sigma12); % draw sample of N0
                NE(t) = N(t) - N0(t);
            else
                Z(t)=1; N0(t)=poissrnd(L(t)); NE(t)=nbinrnd(prior.aE,prior.bE/(1+prior.bE));
            end
        else
            if (N(t)~=-1)
                Z(t)=0; N0(t)=N(t); NE(t)=0;              % no event at time t
            else
                Z(t)=0; N0(t)=poissrnd(L(t)); NE(t)=0;
            end
        end
        ptmp = zeros(2,1); ptmp(Z(t)+1) = 1;    % merge backward influence
        if (t>1), p(:,t-1) = p(:,t-1).*(M'*ptmp); p(:,t-1)=p(:,t-1)/sum(p(:,t-1)); end;
    end
    PE = p;
%     assert ( isequal(N, N0));

function [M] = draw_M_Z(Z,prior)
    % GIVEN Z, SAMPLE M
    n01 = length(find(Z(1:end-1)==0 & Z(2:end)==1)); n0=length(find(Z(1:end-1)==0));
    n10 = length(find(Z(1:end-1)==1 & Z(2:end)==0)); n1=length(find(Z(1:end-1)==1));
    z0 = betarnd(n01+prior.z01, n0-n01+prior.z00);
    z1 = betarnd(n10+prior.z10, n1-n10+prior.z11);
    M = [1-z0, z1; z0, 1-z1];

function [Bmu,Bsigma,Mu0,Sigma0] = draw_L_N0(N,N0,prior,EQUIV)
    Nd=7;   Nh=size(N0,1);
    NE = N - N0;
    Bmu = zeros(size(N0));
    Bsigma = zeros(size(N0));
    
    %update mu0, sigma0
%     assert ( ~isempty(find(NE~=0,1)))
    [mu, sigma] = get_post_para(NE, prior.mu0, prior.sigma0);
    Mu0 = mu;
    Sigma0 = sigma;
    
    %day effect
    Dmu = zeros(1,Nd);
    Dsigma = zeros(1,Nd);
    Hmu = zeros(Nh,Nd);
    Hsigma = zeros(Nh,Nd);
    
    for d=1:length(Dmu)
        daysum_d = sum(N0(:,d:7:end));
        daysum_d = daysum_d/Nh;
        [mu, sigma] = get_post_para(daysum_d, prior.Dmu(d), prior.Dsigma(d));
        Dmu(d) = mu;
        Dsigma(d) = sigma;
        %ToD effect
        for h=1:size(Hmu,1)   %interval
            hour_sample = N0(h,d:7:end);
            hour_sample = hour_sample - daysum_d;
            [mu, sigma] = get_post_para(hour_sample, prior.Hmu(h,d), prior.Hsigma(h,d));
            Hmu(h,d) = mu;
            Hsigma(h,d) = sigma;
        end  
    end
       
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
    
    % posterior hyperparameters for mu^B_t
    for d=1:size(Bmu,2)
        dd = mod(d-1,7)+1; 
        for h=1:size(Bmu,1)
            Bmu(h,d) = Dmu(dd) + Hmu(h,dd); 
            Bsigma(h,d) = sqrt( Dsigma(dd)^2 + Hsigma(h,dd)^2 ); 
        end
    end
    
    %test block
%     subplot(2,1,1)
%     hold on
%     plot(reshape(N,1,[]),'r')
%     plot(reshape(Bmu,1,[]),'k')
%     subplot(2,1,2)
%     plot(reshape(Bsigma,1,[]),'k')
%     plot(reshape(repmat(Dmu,Nh,1) + Hmu,1,[]),'k')
%     HeatMap(sqrt( repmat(Dsigma,Nh,1).^2 + Hsigma.^2 ))
%     Bmu_test = repmat(Dmu,Nh,1) + Hmu;
%     Bmu_test = repmat(Bmu_test,1,4);
%     Bsigma_test = sqrt( repmat(Dsigma,Nh,1).^2 + Hsigma.^2 );
%     Bsigma_test = repmat(Bsigma_test,1,4);
%     assert ( isequal(Bmu,Bmu_test) );
%     assert ( isequal(Bsigma,Bsigma_test) );
    
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
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end

function data = generate_base(D,H)
    data = repmat(D,96,1);
    data = data+ H;