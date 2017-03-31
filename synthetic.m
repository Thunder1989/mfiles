%% synthetic experiments for lassso-based clustering
clear
clc
close all
d = 100; %dimension
run = 1;
p = 1; %bernouli probability
%alpha = 0.85; % spectral norm of A
res1 = [];
res2 = [];
k = 15
%for k = 5:5:25 %# of clusters
for alpha = 0.8:0.2:0.8
% for alpha = 0.1:0.2:1.0 
    A_ = cell(k,1);
    count = [ceil(d/k)*ones(1,mod(d,k)), floor(d/k)*ones(1,k-mod(d,k))]; % # of TS in each cluster
%     tmp = (1:k)';
%     tmp = repmat(tmp,1,d/k);
%     gt = reshape(tmp',1,[]); %ground truth cluster id
    for T = 700:200:700 %sample size
%     for T = 100:200:1000 %sample size
        for beta = 0.45:0.05:0.45
%         for beta = 0.05:0.05:0.55
%         for beta = logspace(-3, 0.8, 10) 
            %generate A - block diagonal
            rand_sum = 0;
            vio_sum = 0;
            for c = 1:run
                gt = [];
                for i = 1:numel(A_)
                    num = count(i);
                    A_{i} = randi([1 10],1) * binornd(1,p,num,num); %randomly assign 1s to each block
                    gt = [gt i*ones(1,num)]; % cluster id
                end
                A = blkdiag(A_{:});
                A = eye(d) + A-diag(diag(A));
                %A = A-diag(diag(A));
                %rescale 2norm of A
                A = A/norm(A) * alpha;

                %generate Sigma
                tmp = zeros(d,d);
                tmp(:) = 0.05; %set all to 0.05, larger is better?
                Sigma = eye(d) + tmp-diag(diag(tmp));
                % Sigma = eye(d);
                Sigma = Sigma/norm(Sigma);

                %get Psi
                Psi = Sigma - A'*Sigma*A;
                %generate noise Z
                mu = zeros(T,d);
                Z = mvnrnd(mu, Psi);

                %generate X with VAR, each column is a sample, each row is a time point T
                mu = zeros(1,d);
                X1 = mvnrnd(mu, Sigma);
                X = zeros(T,d);
                X(1,:) = X1;
                for i = 2:T
                    X(i,:) = A*X(i-1,:)' + Z(i,:)';
                end

                %lasso based clustering
                W = zeros(d,d); %regression coefficient matrix
                X_S = X(1:T-1,:);
                X_T = X(2:T,:);

                b = beta * sqrt(log(d)/T);
%                 b = 0.03;
                for i = 1:d
                    dest = X_T(:,i); %i-th sample, T-1 by 1
                    src = X_S;
            %         src(:,i) = []; %all other samples, T-1 by d-1
                    % 0.015~0.03 for max-min normalization, 0.06 gives no zero rows
                    % 0.14~0.2 for u-std normalization, 0.14 gives no zero rows - 0.02-0.14 all reasonable
                    coef = lasso(src, dest, 'Lambda', b); %d-1 by 1
                    W(i,:) = coef;
                end

            %     loss = norm(W - A)
            %     figure
            %     subplot(1,2,1)
            %     imagesc(A);
            %     subplot(1,2,2)
            %     imagesc(W);

                %w = w + eye(d);
                [g_, idx] = sort(gt); %ground truth indexing
                W_ = W(idx,idx); %re-ordered by true type ID

                %compute violation rate
                B_ = cell(k,1);
                for i = 1:numel(B_)
                    num = count(i);
                    B_{i} = ones(num);
                end
                B = blkdiag(B_{:});
                blk_idx = logical(B);
                C = abs(W_);
                inner_c_sum = sum(sum(C(blk_idx)));
                outer_c_sum = sum(sum(C)) - inner_c_sum;
                vio = outer_c_sum / inner_c_sum;
                vio_sum = vio_sum + vio;
                spy(W_)
%                 continue;
                
%                 W_ = A;
                W_ = max(W_, W_'); %symmetrize w_, N by N
                D = diag(sum(W_,2));
                L = D - W_; %unormalized Laplacian
                [evc, evl] = eig(L); %each column of evc is an eigenvector
                idx = find(diag(evl)>=0);
                input = evc(:,idx(1:k));
%                 input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
                c_idx = kmeans(input,k);
                ari = adjrand(c_idx, g_)
                rand_sum = rand_sum + ari;
            %     figure
%                 spy(W_)
            %     res = [res adjrand(c_idx, gt)];
            end
            res1 = [res1; vio_sum/run];
            res2 = [res2; rand_sum/run];
        end
    end
end

save('synthetic.mat','X');
% save('vio_vs_alpha','res1');
% save('rand_vs_alpha','res2');

%% CC
% load('corr_syn.mat');
res = [];
for i=1:run
    corr = corrcoef(X);
    corr = corr - eye(size(corr));
    W_ = corr;
    run = 10;
    W_ = max(W_, W_'); %symmetrize w_, N by N
    D = diag(sum(W_,2));
    L = D - W_; %unormalized Laplacian
    [evc, evl] = eig(L); %each column of evc is an eigenvector
    idx = find(diag(evl)>=0);
    input = evc(:,idx(1:k));
    % input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
    c_idx = kmeans(input,k);
    ari = adjrand(c_idx, gt);
    res = [res ari];
end
%% Cosine
options = [];
options.NeighborMode = 'KNN';
options.WeightMode = 'Cosine';
options.bTrueKNN = 1;
% options.t = 1;

res_ = [];
for i = 7:2:7
    options.k = i;
    W_c = constructW(X', options);
    W_ = W_c(idx,idx);
%     figure
%     spy(W_)
    W_ = max(W_, W_'); %symmetric N by N
    D = diag(sum(W_,2));
    L = D - W_;
    [evc, evl] = eigs(L,4,'sa'); %N by N, each column is a principal component
    c_idx = kmeans(evc,4,'Distance','cosine');
    ari = adjrand(c_idx, gt)
    res_ = [res_ adjrand(c_idx, gt)];
end

%% ACF
load('acf_syn.mat');
W_ = acf;
W_ = W_ - min(W_(:));
W_ = W_ ./ max(W_(:)); %normalized to (0,1)
W_ = 1 - W_; %convert distance to correlation-like
D = diag(sum(W_,2));
L = D - W_; %unormalized Laplacian
[evc, evl] = eig(L); %each column of evc is an eigenvector
idx = find(diag(evl)>=0);
input = evc(:,idx(1:k));
% input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
c_idx = kmeans(input,k);
ari = adjrand(c_idx, g_)