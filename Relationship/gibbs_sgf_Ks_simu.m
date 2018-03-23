clc
tic

load('320_output.mat') %load the Q R from previous experiments

%----generate data----
rep = 20;
len = repmat(10,1,rep);
vel = repmat([0 30 -30 0],1,rep/4);
state = repmat([1 2 2 1],1,rep/4);

F = [1 1; 0 1];
H = [1 0; 0 1];

Q_star = Q_ahu{2}';
R_star = R_ahu{2}';

Z_star = zeros(sum(len),1);
Z_star(11:10:end) = 1;
Z_star(1:40:end) = 0; %checked

y_prev = rand*100; %initial Y_0

Y = [];
X = [];

%Y ~ N(FY,Q), X ~ N(HX,R)
for i = 1:length(len)
    count = len(i);
    for t = 1:count
        tmp_y = mvnrnd( F*[y_prev vel(i)]',  Q_star{state(i)} );
        Y = [Y; tmp_y];
        y_cur = tmp_y(1);
        vy_cur = tmp_y(2);

        tmp_x = mvnrnd( H*[y_cur vy_cur]', R_star{state(i)} );
        X = [X; tmp_x];
        y_prev = y_cur;
    end

end


%------initialization------'
Ks = 2;
F_switch = 1;
Z = randi(Ks, size(X,1), 1); %initialize z seq
M = rand(Ks);
M = bsxfun(@rdivide, M, sum(M)); %transition matrix
F_ = { [1 1; 0 1], [1 0; 0 1] }; %1-steady, 2-event
F = F_{1}; %F will be estimated
H = [1 0; 0 1];

Q = mat2cell(rand(2*Ks,2), 2*ones(1,Ks), 2);
R = mat2cell(rand(2*Ks,2), 2*ones(1,Ks), 2);
non_event_idx = get_non_event_i(Q);
for i = 1:Ks
    if F_switch
        if i==non_event_idx
            F = F_{1};
        else
            F = F_{2};
        end
    end
    [Q{i}, R{i}] = get_Q_R(i,X,Y,Z,F,H);
end

Z_err = [];
Q_err = [];
R_err = [];
K = 10; % iters for EM
N = 100; % samples per iter
p_data = zeros(K,1);
Z_sample = repmat(Z,1,N);
p_tmp = zeros(Ks,1);

for iter = 1:10
    fprintf('----iter #%d----\n',iter);

    %------EM------
    for k = 1:K
        %fprintf('--------------iter# %d--------------\n',k);

        %---E step---
        non_event_idx = get_non_event_i(Q);

        %sample Z
        for n = 2:N
            for t = 2:length(Z)-1
                if F_switch
                    if Z_sample(t,n-1)==non_event_idx
                        F = F_{1};
                    else
                        F = F_{2};
                    end
                end

                for i = 1:length(p_tmp)
                    try
                        p_tmp(i) = M( i, Z_sample(t-1,n) ) * M( Z_sample(t+1,n-1), i ) * mvnpdf(X(t,:)', H*Y(t,:)', R{i}); %* mvnpdf(Y(t,:)', F*Y(t-1,:)', Q{i});
                    catch
                        fprintf('z sampling error!\n');
                        p_tmp(i) = randi(Ks);
                    end
                end
                p_tmp = p_tmp/sum(p_tmp);
                Z_sample(t,n) = find(cumsum(p_tmp) >= rand(1), 1);
            end
        end

        %sample Y
        Y_sample = repmat(Y,1,1,N);
%         p_tmp = zeros(N,1); %for data likelihood
        for n = 2:N
            for t = 2:size(Y,1)-1 %todo: shuffle the order

%                 Q_prev = Q{Z(t-1)+1};
                Q_next = Q{Z(t+1)}; %Question: use latest Z?
                Q_cur = Q{Z(t)};
                R_cur = R{Z(t)};

                if F_switch
                    if Z(t)==non_event_idx
                        F = F_{1};
                    else
                        F = F_{2};
                    end
                end

                if Z(t-1)==non_event_idx && Z(t)~=non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)~=non_event_idx && Z(t)==non_event_idx && Z(t+1)~=non_event_idx
                    %010 or 101: p(x_t|y_t)
                    Sigma_e = H^-1*R_cur*(H^-1)';
                    Mu_e = H^-1*X(t,:)';
                elseif Z(t-1)~=non_event_idx && Z(t)~=non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)==non_event_idx && Z(t)==non_event_idx && Z(t+1)~=non_event_idx
                    %110 or 001: p(x_t|y_t) p(y_t|y_t-1)
                    Sigma_e = ( Q_cur^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * Q_cur^-1 * (F*Y_sample(t-1,:,n)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');
                elseif Z(t-1)~=non_event_idx && Z(t)==non_event_idx && Z(t+1)==non_event_idx ...
                    || Z(t-1)==non_event_idx && Z(t)~=non_event_idx && Z(t+1)~=non_event_idx
                    %100 or 011: p(x_t|y_t) p(y_t|y_t+1)
                    Sigma_e = ( (F^-1*Q_next*(F^-1)')^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                    Mu_e = Sigma_e * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y_sample(t+1,:,n-1)') + Sigma_e * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');    
                end

                %normal case - product of 3 gaussians, always used for sampling position 
                Sigma_ = ( Q_cur^-1 + (F^-1*Q_next*(F^-1)')^-1 )^-1;
                Mu_ = Sigma_* Q_cur^-1 * (F*Y_sample(t-1,:,n)') + Sigma_ * (F^-1*Q_next*(F^-1)')^-1 * (F^-1*Y_sample(t+1,:,n-1)');
                Sigma = ( Sigma_^-1 + (H^-1*R_cur*(H^-1)')^-1 )^-1;
                Mu = Sigma * Sigma_^-1 * Mu_ + Sigma * (H^-1*R_cur*(H^-1)')^-1 * (H^-1*X(t,:)');

                if ( Z(t)~=Z(t-1) || Z(t)~=Z(t+1) ) && ...
                ( Z(t)==non_event_idx || Z(t-1)==non_event_idx || Z(t+1)==non_event_idx )
                    Y_sample(t,1,n) = normrnd(Mu(1), Sigma(1,1));
                    Y_sample(t,2,n) = normrnd(Mu_e(2), Sigma_e(2,2));
                else
                    Y_sample(t,:,n) = mvnrnd(Mu, Sigma);
                end

                %data likelihood
%                 p = 0;
%                 for t=2:size(Y,1)-1 
%                     p = p + log( mvnpdf(Y_sample(t,:,n), Y_sample(t-1,:,n)*F', Q) )...
%                         + log( mvnpdf(X(t,:), Y_sample(t,:,n)*H', R) );
%                 end
%                 p_tmp(n) = p;
            
            end          
        end
        p_data(k) = mean(p_tmp);

        %---M step---
        Y = mean(Y_sample, 3);
        Z = mode(Z_sample(:,end-10:end), 2);
        M = get_M(Z_sample, Ks);

        for i = 1:Ks
            if F_switch
                if i==non_event_idx
                    F = F_{1};
                else
                    F = F_{2};
                end
            end
            [Q{i}, R{i}] = get_Q_R(i,X,Y,Z_sample,F,H); %check: use Y_sample?
        end

    end
%         Z = mean(Z_sample(:,end-20:end),2);
    Z = Z_sample;      

    %calculate errors
    Z = mode(Z(:,11:2:end),2) - 1; %atn: Z

    assert( length(Z) == length(Z_star) );
    Z_ = remap_event(Z); %atn: Z
    Z_err = [Z_err sum(abs(Z_ - Z_star))/length(Z_)];

    tmp = cellfun(@minus, Q, Q_star, 'UniformOutput', false);
    tmp = cellfun(@abs, tmp, 'UniformOutput', false);
    tmp = cell2mat( cellfun(@(x) sum(x(:)), tmp, 'UniformOutput', false) );
    Q_err = [Q_err; tmp(:)'];

    tmp = cellfun(@minus, R, R_star, 'UniformOutput', false);
    tmp = cellfun(@abs, tmp, 'UniformOutput', false);
    tmp = cell2mat( cellfun(@(x) sum(x(:)), tmp, 'UniformOutput', false) );
    R_err = [R_err; tmp(:)'];

    flip = ~isequal(Z, Z_);
    for i = 1:length(Z)
        if rand>=0.2
            Z(i) = abs(1*flip - Z_star(i));
        end
    end

    Z = Z + 1; %atn: Z
    M = get_M(Z(:), Ks); %atn: Z
    for i = 1:Ks
        if F_switch
            if i==non_event_idx
                F = F_{1};
            else
                F = F_{2};
            end
        end
        [Q{i}, R{i}] = get_Q_R(i,X,Y,Z,F,H); %atn: Z
    end
    
end

data = X(:,1);
res = Y;
Z = remap_event(Z-1);
figure
hold on
yyaxis left
%events in shade
stem(1:length(Z), Z*max(data), 'Marker','None', 'LineWidth', 5, 'Color',[.8 .8 .8]);
%original data
plot(1:length(res), data, 'k--', 'LineWidth', 1.5)
%filtered data
%     plot(1:length(res), res(:,1),'g-')
yyaxis right
%velocity
plot(1:length(res), res(:,2), 'r-','LineWidth', 2)

figure
subplot(3,1,1)
plot(Z_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)
subplot(3,1,2)
plot(Q_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)
subplot(3,1,3)
plot(R_err, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 8)
    
toc