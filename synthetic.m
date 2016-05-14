%% synthetic experiments for lassso-based clustering
clear
clc
close all
d = 100; %dimension
run = 20;
p = 1; %bernouli probability
alpha = 0.85; % spectral norm of A
res1 = [];
res2 = [];
for k = 20:5:25 %# of clusters
    A_ = cell(k,1);
    count = [ceil(d/k)*ones(1,mod(d,k)), floor(d/k)*ones(1,k-mod(d,k))];
%     tmp = (1:k)';
%     tmp = repmat(tmp,1,d/k);
%     gt = reshape(tmp',1,[]); %ground truth cluster id
    for T = 1000:200:1000 %sample size
%         for beta = 0.15:0.05:0.55
        for beta = logspace(-3, 0.8, 10) 
            %generate A - block diagonal
            rand_sum = 0;
            vio_sum = 0;
            for c = 1:run
                gt = [];
                for i = 1:numel(A_)
                    num = count(i);
                    A_{i} = binornd(1,p,num,num); %equal-sized cluster
                    gt = [gt i*ones(1,num)];
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

%                 b = beta * sqrt(log(d)/T);
                b = 1.0 * sqrt(log(d)/T);
%                 b = 0.03;
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
                vio = outer_c_sum / inner_c_sum
                vio_sum = vio_sum + vio;
%                 spy(W_)
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
                ari = adjrand(c_idx, g_);
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