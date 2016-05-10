%% synthetic experiments for lassso-based clustering
clear
clc
close all
k = 10; %# of clusters
d = 100; %dimension
T = 1000; %sample size

%generate A - block diagonal
p = 1; %bernouli probability
A_ = cell(k,1);
rand_sum = 0;
run = 5;
for c = 1:run
    for i = 1:numel(A_)
        A_{i} = binornd(1,p,d/k,d/k); %equal-sized cluster
    end
    A = blkdiag(A_{:});
    A = eye(d) + A-diag(diag(A));
    %A = A-diag(diag(A));
    %rescale by 2norm
    alpha = 0.85;
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
    tmp = (1:k)';
    tmp = repmat(tmp,1,d/k);
    gt = reshape(tmp',1,[]); %ground truth cluster id

    W = zeros(d,d); %regression coefficient matrix
    X_S = X(1:T-1,:);
    X_T = X(2:T,:);

    b  = 0.55*sqrt(log(d)/T);
    for i = 1:d
        cur = X_T(:,i); %i-th sample, T-1 by 1
        src = X_S;
%         src(:,i) = []; %all other samples, T-1 by d-1
        % 0.015~0.03 for max-min normalization, 0.06 gives no zero rows
        % 0.14~0.2 for u-std normalization, 0.14 gives no zero rows - 0.02-0.14 all reasonable
        coef = lasso(src, cur, 'Lambda', b); %d-1 by 1
%         coef = lasso(src, cur); %d-1 by 1
%         idx = 1;
%         for j = 1:length(coef)
%             if(idx==i) 
%                 idx = idx+1;
%             end
%             w(i,idx) = coef(j);
%             idx = idx+1;
%         end
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
    W_ = max(W_, W_'); %symmetrize w_, N by N
    D = diag(sum(W_,2));
    L = D - W_; %unormalized Laplacian
    [evc, evl] = eig(L); %each column of evc is an eigenvector
    idx = find(diag(evl)>=0);
    input = evc(:,idx(1:k));
%     input = evc(:,1:idx(k)); %trick: including extra negative and zero evls, slightly better
    c_idx = kmeans(input,k);
    rand_index = adjrand(c_idx, g_);
    rand_sum = rand_sum+rand_index;
%     figure
%     spy(w_)
%     res = [res adjrand(c_idx, gt)];
end
rand_sum/run
