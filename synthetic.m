%% synthetic experiments for lassso-based clustering
% clear
clc
k = 5; %# of clusters
d = 50; %dimension
T = 100; %sample size
X = zeros(T,d); %input data

%generate A - block diagonal
p = 0.5; %bernouli probability
A_ = cell(k,1); 
for i=1:numel(A_)
    A_{i} = binornd(1,p,d/k,d/k); %equal-sized cluster
end
A = blkdiag(A_{:});
%rescale by 2norm
alpha = 0.5;
A = A/norm(A) * alpha;

%generate Sigma
tmp = zeros(d,d);
tmp(:) = 0.1; %set to 0.1 for all
% Sigma = eye(d) + tmp-diag(diag(tmp));
Sigma = eye(d);
Sigma = Sigma/norm(Sigma);

%get Psi
Psi = Sigma - A'*Sigma*A;
%generate noise Z
mu = zeros(T,d);
Z = mvnrnd(mu, Psi);

%generate X with VAR, each column is a sample, each row is a time point T
mu = zeros(1,d);
X1 = mvnrnd(mu, Sigma);
X(1,:) = X1;
for i = 2:T
    X(i,:) = (A*X(i-1,:)' + Z(i,:)')';
end

%% lasso based clustering
tmp = (1:k)';
tmp = repmat(tmp,1,d/k);
gt = reshape(tmp',1,[]); %ground truth cluster id

% %normalization
% for i = 1:d
%     tmp = X(:,i);
% %     t_min = min(tmp);
% %     t_max = max(tmp);
% %     data(i,:) = (tmp-t_min) / (t_max-t_min);
%     mu = mean(tmp);
%     sigma = std(tmp);
%     X(:,i) = (tmp-mu) / sigma;
% end

w = zeros(d,d); %regression coefficient matrix
X_S = X(1:T-1,:);
X_T = X(2:T,:);
for b = 0.1:0.02:0.1
    b = 0.15;
    for i = 1:d
        cur = X_T(:,i); %i-th sample, T by 1
        src = X_S;
        src(:,i) = []; %all other samples, T by d-1
        % 0.015~0.03 for max-min normalization, 0.06 gives no zero rows
        % 0.14~0.2 for u-std normalization, 0.14 gives no zero rows - 0.02-0.14 all resonable
        coef = lasso(src, cur, 'Lambda', b); %N-1 by 1
        idx = 1;
        for j = 1:length(coef)
            if(idx==i) 
                idx = idx+1;
            end
            w(i,idx) = coef(j);
            idx = idx+1;
        end
    end

    [g_, idx] = sort(gt); %ground truth indexing
    w_ = w(idx,idx); %re-ordered by true type ID
    w_ = max(w_, w_'); %symmetrize w_, N by N
    D = diag(sum(w_,2));
%     find(diag(D)==0)
%     fprintf('all zero rows found!\n');
%     pause
    l = D - w_; %unormalized Laplacian
    [evc, evl] = eig(l); %N by N, each column in evc is an eigenvector
    idx = find(diag(evl)>0);
    input = evc(:,idx(1:k));
    c_idx = kmeans(input,k,'Distance','cosine');
    adjrand(c_idx, g_)
%     figure
    spy(w_)
%     res = [res adjrand(c_idx, gt)];
end
