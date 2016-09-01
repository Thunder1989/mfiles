clear
clc
% timescale = [1,3,5,7,14,21,28];
path = '/Users/hdz_1989/Documents/Dropbox/SDB/KETI/';%path to the data folers
folder = dir(path);
FN = 'coef_raw.txt';
FID = fopen(FN,'w');

% for k = 1:length(timescale)
%     T = timescale(k);
%     th1 = 8*T;%6h
%     th2 = 96*T;%0.5h
    data_length = 275770;
    %search for the shortest data length
%     for fn=5:size(folder,1)-1
%         path1 = strcat(path, folder(fn).name,'/');%path to the csv files
%         file = dir(strcat(path1,'*.csv'));
%         filename = [path1, file(1).name];
% %         rd = csvread(filename);
% %         ts = rd(:,1);
% %         data = rd(:,2);
% %         data = data(ts<(ts(1)+T*24*60*60));
%         data = csvread(filename);
%         if length(data)<=data_length
%             data_length = length(data);
%         end
%     end
    data_aggr = zeros(4*(size(folder,1)-5), data_length);

    %load data => direct xocrr as a baseline
    for fn=5:size(folder,1)-1
        path1 = strcat(path, folder(fn).name,'/');%path to the csv files
        file = dir(strcat(path1,'*.csv'));
        num = size(file,1);
        for n = 1:num
            filename = [path1, file(n).name];
%             rd = csvread(filename);
%             ts = rd(:,1);
%             data = rd(:,2);
%             data = data(ts<(ts(1)+T*24*60*60));
            data = csvread(filename);
            len = length(data);
            if len > data_length
                [p,q] = rat(data_length/len);
                data = resample(data, p, q);
                data = data(1:data_length);
            end
            if n~=3
                tmp=data(1);
                for k=2:length(data)
                    if abs(abs(data(k))-abs(tmp))>3*min(abs(tmp),abs(data(k))) && tmp~=0
                        data(k) = tmp;
                        continue;
                    end
                    tmp = data(k);
                end
            end
            data_aggr(((fn-5)*4+n),:) = data;
        end
    end
    
    %compute xcorr
    num = size(data_aggr,1);
    xcor_mat = eye(num);
    for i=1:num
      for j=1:(i-1)
        coef = corrcoef(data_aggr(i,:),data_aggr(j,:));
        xcor_mat(i, j) = coef(1,2);
        xcor_mat(j, i) = xcor_mat(i, j);
      end
    end
%     fprintf(FID, '==========xcormat obtained with %d days data==========\n', T);
    for i=1:num
        for j=1:num
            fprintf(FID, '%.2f\t', xcor_mat(i, j));
        end
        fprintf(FID,'\n');
    end
% end
fclose(FID);

%%
k = 51;
W_ = corr;
W_ = max(W_, W_'); %symmetrize w_, N by N
D = diag(sum(W_,2));
L = D - W_; %unormalized Laplacian
[evc, evl] = eig(L); %each column of evc is an eigenvector
idx = find(diag(evl)>=0);
input = evc(:,idx(1:k));
HeatMap(abs(input))

%% permutation test for significance
%permutations created by random sampling - one stream per type per room,
%without replacement
clc
N = 1000000; %number of trials
k = 51;
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
corr_ = abs(corr_);

%orginal score
score = 0;
for t = 0:k-1
    idx = 4*t+1;
    corr_cur = corr_(idx:idx+3, idx:idx+3);
    score = score + sum(sum(corr_cur));
end

% reorder by type for resampling
idx = 1:4;
idx = repmat(idx,1,51);
[Y,I] = sort(idx);
corr_ = corr_(I,I);

res = zeros(N,1);
ctr = 0;
for t = 1:N
    idx = zeros(k*4,1);
    for i = 0:3
        idx(k*i+1:k*(i+1)) = randperm(k) + k*i;
    end
    I = reshape(reshape(idx, k, [])', k*4, 1); %get the idx for resampling one type each per room
    corr_tmp = corr_(I,I);
    tmp = 0;
    for i = 0:k-1
        idx = 4*i+1;
        corr_cur = corr_tmp(idx:idx+3, idx:idx+3);
        tmp = tmp + sum(sum(corr_cur));
    end
    res(t) = tmp;
    if tmp >= score
        ctr = ctr+1;
    end
end

[n, x] = hist(res, 20);
figure
hold on
bar(x, n/N);
plot([score, score], [0,max(n)/N], 'k--')
significance = ctr / N

%% permutation test, with permutations created by swapping based on the ground truth arrangement and only permutate within type
% clear
clc
% load('keti_corr_typecleared.mat')
room_list = importdata('keti_room_list');
type_list = {'co2', 'hum', 'light', 'temp'};
k = 51;
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
corr_ = abs(corr_);

%orginal score
score = 0;
for t = 0:k-1
    idx = 4*t+1;
    corr_cur = corr_(idx:idx+3, idx:idx+3);
    score = score + sum(sum(corr_cur));
end

res = [];
ctr = 0;
times = 0;
idx = 1:k*4;
idx = reshape(idx, 4, []);
p = 3; %length of permutation subsequence
FN = sprintf('permutation_log_p%s.txt',num2str(p));
FID = fopen(FN,'w');
count = zeros(4,1);
for t = 1:4
    row = idx(t,:); %permute the t-th type, which is the t-th row
    for m = 1:length(row)-p+1
        combo = perms(row(m:m+p-1));
        for i = 1:size(combo,1)
            times = times + 1;
            if row(m:m+p-1) == combo(i,:)
                continue
            end
            row_ = row;
            row_(m:m+p-1) = combo(i,:); %set m to m+p-1 with the new permutation
            idx_ = idx;
            idx_(t,:) = row_;
            I = reshape(idx_, k*4, 1); %flatten the reshaped index array
            corr_tmp = corr_(I,I);
            
            tmp = 0;
            for j = 0:k-1
                idx_ = 4*j+1;
                corr_cur = corr_tmp(idx_:idx_+3, idx_:idx_+3);
                tmp = tmp + sum(sum(corr_cur));
            end
            res = [res tmp];
 
            if tmp >= score
                count(t) = count(t) + 1;
                ctr = ctr+1;
                room_id = floor((combo(i,:)-1)/4)+1;
                fprintf(FID, '%s--%s', type_list{t}, num2str(room_id));
                for n = 1:length(room_id)
                    fprintf(FID, ',%s', room_list{room_id(n)});
                end
                fprintf(FID, '\n');
            end
        end
    end
end
fclose('all');

% [n, x] = hist(res, 20);
% figure
% hold on
% bar(x, n/times);
% plot([score, score], [0,max(n)/times], 'r--', 'LineWidth', 2)
% significance = ctr / times;
% xlabel(['p-value=',num2str(significance),' on ',num2str(times),' with k=',num2str(p)]);
[count'/ctr ctr/ times], times

%% permutation test including across type permutations
% clear
clc
% load('keti_corr_typecleared.mat')
room_list = importdata('keti_room_list');
type_list = {'co2', 'hum', 'light', 'temp'};
k = 51;
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
corr_ = abs(corr_);

%orginal score
score = 0;
for t = 0:k-1
    idx = 4*t+1;
    corr_cur = corr_(idx:idx+3, idx:idx+3);
    score = score + sum(sum(corr_cur));
end

res = [];
ctr = 0;
times = 0;
idx = 1:k*4;
% idx = reshape(idx, 4, []);
p = 2; %length of permutation subsequence
FN = sprintf('permutation_log_p%s.txt',num2str(p));
FID = fopen(FN,'w');
count = zeros(4,1);
combo = nchoosek(1:4,1);
for r = 1:k-2
    for i = 1:length(combo)
        for j = 1:length(combo)
            for l = 1:length(combo)
                idx_1 = 4*(r-1) + i;
                idx_2 = 4*r + j;
                idx_3 = 4*(r+1) + l;
                pm = perms([idx_1,idx_2,idx_3]);
                for n = 1:size(pm,1)                    
                    times = times + 1;
                    I = idx;
                    pm_cur = pm(n,:);
                    if [idx_1,idx_2,idx_3] == pm_cur
                        continue
                    end
                    I(idx_1) = pm_cur(1);
                    I(idx_2) = pm_cur(2);
                    I(idx_3) = pm_cur(3);
                    corr_tmp = corr_(I,I);

                    tmp = 0;
                    for m = 0:k-1
                        idx_ = 4*m + 1;
                        corr_cur = corr_tmp(idx_:idx_+3, idx_:idx_+3);
                        tmp = tmp + sum(sum(corr_cur));
                    end
                    res = [res tmp];

                    if tmp >= score
        %                 count(r) = count(r) + 1;
                        ctr = ctr+1;
        %                 room_id = floor((combo(i,:)-1)/4)+1;
        %                 fprintf(FID, '%s--%s', type_list{r}, num2str(room_id));
        %                 for n = 1:length(room_id)
        %                     fprintf(FID, ',%s', room_list{room_id(n)});
        %                 end
        %                 fprintf(FID, '\n');
               
                    end
                end
            end
        end
    end
end
fclose('all');

% [n, x] = hist(res, 20);
% figure
% hold on
% bar(x, n/times);
% plot([score, score], [0,max(n)/times], 'r--', 'LineWidth', 2)
significance = ctr / times
% xlabel(['p-value=',num2str(significance),' on ',num2str(times),' with k=',num2str(p)]);


%% statistics
clc
ctr = 0;
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
for i=1:length(corr_)
    tmp = abs(corr_(i,:));
    id = find(tmp==max(tmp),1);
    if id <= 4*floor((i-1)/4)+4 && id >= 4*floor((i-1)/4)+1
            ctr = ctr+1;
    end

%     [v, id] = sort(tmp,'descend');
%     for j=1:3
%         if id(j) <= 4*floor((i-1)/4)+4 && id(j) >= 4*floor((i-1)/4)+1
%             ctr = ctr+1;
%             continue;
%         end
%     end
end
ctr

%% ILP
clear
load('keti_corr_typecleared.mat')

%get input sample weight matrix
ids = [1,15,25];
idx = [];
for i = 1:length(ids)
    id = ids(i);
    idx_ = 4*(id-1) + 1;
    idx = [idx idx_:idx_+3];
end
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
corr_ = abs(corr_);
corr_tmp = corr_(idx,idx);
% corr_tmp = corr_;

nEdges = length(corr_tmp)*3/2;
nNodes = length(corr_tmp);
idxs = nchoosek(1:nNodes, 2);
tmp = []
%remove edge connecting the same type
for i = 1:length(idxs)
    tmp_ = idxs(i,:);
    if mod(tmp_(1),4) ~= mod(tmp_(2),4)
        tmp = [tmp; tmp_];
    end
end
idxs = tmp;

% Get all the weights.
weight = zeros(length(idxs),1);
for i = 1:length(weight)
    weight(i) = corr_tmp(idxs(i,1),idxs(i,2));
end
% the sum of score is dist'*x_tsp, where x_tsp is the binary solution vector.

% Equality Constraints
% The first enforces that there must be N*3 edges total. 
% The second enforces that each node must have 3 edges attached to it.
% Specify the first constraint, Aeq * x_tsp = beq.
Aeq = spones(1:length(idxs)); % Add up the number of edges
beq = nEdges;
% To specify the second type of equality constraint, that there needs to be two trips attached to each stop, extend the Aeq matrix as sparse.
Aeq = [Aeq; spalloc(nNodes,length(idxs),nNodes*(nNodes-1))]; % allocate a sparse matrix
for ii = 1:nNodes
    whichIdxs = (idxs == ii); % find the edges that include node ii
    Aeq(ii+1,:) = sparse(sum(whichIdxs,2)); % include edges containing ii and add into the constraint matrix
end
beq = [beq; 3*ones(nNodes,1)];

% All decision variables are binary. 
% Now, set the intcon argument to the number of decision variables, 
% put a lower bound of 0 on each, and an upper bound of 1.
len_w = length(weight);
intcon = 1:len_w;
lb = zeros(len_w,1);
ub = ones(len_w,1);

% Optimize with intlinprog
opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(-weight,intcon,[],[],Aeq,beq,lb,ub,opts);
idxs(x_tsp==1,:)

%% new ILP formulation
clear
clc
load('keti_corr_typecleared.mat')

%get input sample weight matrix
ids = [1,10];
idx = [];
for i = 1:length(ids)
    id = ids(i);
    idx_ = 4*(id-1) + 1;
    idx = [idx idx_:idx_+3];
end
corr_ = corr - diag(diag(corr)); % remove 1s on the diagonal
corr_ = abs(corr_);
corr_tmp = corr_(idx,idx);
% corr_tmp = corr_;

M = 4; % # of sensors per cluster
nNodes = length(corr_tmp);
nClusters = length(ids);

% integer constraints
lbx = zeros(nClusters,nNodes); % x variables
lby = zeros(nClusters,nNodes,nNodes); % y variables
lb = [lbx(:); lby(:)]; % Column vector lower bound
ub = ones(size(lb)); % Binary variables have lower bound 0, upper bound 1

clearx = zeros(nClusters,nNodes); % 0 for the x variables
cleary = zeros(nClusters,nNodes,nNodes); % 0 for the y variables

% Equality Constraints
% 1) each sensor belongs to 1 cluster/clique
Aeq = spalloc(nNodes,length(lb),nClusters*nNodes); 
counter = 1;
for ii = 1:nNodes
    tmp = clearx;
    tmp(:,ii) = 1;
    addrow = [tmp(:); cleary(:)]';
    Aeq(counter,:) = sparse(addrow);
    counter = counter + 1;
end
beq = ones(nNodes,1);
% 2) each clique has M=4 sensors
Aeq_tmp = spalloc(nClusters,length(lb),nClusters*nNodes); 
counter = 1;
for ii = 1:nClusters
    tmp = clearx;
    tmp(ii,:) = 1;
    addrow = [tmp(:); cleary(:)]';
    Aeq_tmp(counter,:) = sparse(addrow);
    counter = counter + 1;
end
Aeq = [Aeq; Aeq_tmp];
beq = [beq; M*ones(nClusters,1)];
% % 3) each clique has M types of sensors
% Aeq_tmp = spalloc(nClusters*M,length(lb),nClusters*M*nClusters); 
% counter = 1;
% for ii = 1:nClusters
%     for jj = 1:M
%         tmp = clearx;
%         tmp(ii,jj:4:end) = 1;
%         addrow = [tmp(:); cleary(:)]';
%         Aeq_tmp(counter,:) = sparse(addrow);
%         counter = counter + 1;
%     end
% end
% Aeq = [Aeq; Aeq_tmp];
% beq = [beq; ones(nClusters*M,1)];

% Inequality Constraints
% 1) x_iu + x_iv - y_iuv <= 1
A = spalloc(nClusters*nNodes*nNodes,length(lb),nClusters*nNodes*nNodes*3); % each sensor belongs to 1 cluster/clique
counter = 1;
for ii = 1:nClusters
    for jj = 1:nNodes
        for kk = 1:nNodes
            if kk==jj
                continue
            end
            tmpx = clearx;
            tmpx(ii,jj) = 1;
            tmpx(ii,kk) = 1;
            tmpy = cleary;
            tmpy(ii,jj,kk) = -1;
            addrow = [tmpx(:); tmpy(:)]';
            A(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
b = ones(nClusters*nNodes*nNodes,1);
% 2) -x_iu + y_iuv <= 0
Atmp = spalloc(nClusters*nNodes*nNodes,length(lb),nClusters*nNodes*nNodes*2); % each sensor belongs to 1 cluster/clique
counter = 1;
for ii = 1:nClusters
    for jj = 1:nNodes
        tmpx = clearx;
        tmpx(ii,jj) = -1;
        for kk = 1:nNodes
            if kk==jj
                continue
            end
            tmpy = cleary;
            tmpy(ii,jj,kk) = 1;
            addrow = [tmpx(:); tmpy(:)]';
            Atmp(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
A = [A; Atmp];
b = [b; zeros(nClusters*nNodes*nNodes,1)];
% 3) -x_iv + y_iuv <= 0
Atmp = spalloc(nClusters*nNodes*nNodes,length(lb),nClusters*nNodes*nNodes*2); % each sensor belongs to 1 cluster/clique
counter = 1;
for ii = 1:nClusters
    for jj = 1:nNodes
        tmpx = clearx;
        for kk = 1:nNodes            
            if kk==jj
                continue
            end
            tmpx(ii,kk) = -1;
            tmpy = cleary;
            tmpy(ii,jj,kk) = 1;
            addrow = [tmpx(:); tmpy(:)]';
            Atmp(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
A = [A; Atmp];
b = [b; zeros(nClusters*nNodes*nNodes,1)];

% score to maxmize
% weight = cleary;
% for ii = 1:nNodes
%     for jj = 1:nNodes
%         if jj==ii
%             continue
%         end
%         weight(:,ii,jj) = corr_tmp(ii,jj);
%     end
% end
% weight = [clearx(:); weight(:)];

ty = zeros(nClusters,nNodes,nNodes);
weight = [[1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 1]'; ty(:)];

% Optimize with intlinprog
opts = optimoptions('intlinprog','Display','off');
[solution,costopt,exitflag,output] = intlinprog(-weight(1:16),1:length(weight(1:16)),[],[],Aeq(:,1:16),beq,lb(1:16),ub(1:16),opts);
exitflag