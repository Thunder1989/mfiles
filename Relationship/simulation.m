clear
clc
load('320_output');
len = [50 50 50 50];
vel = [0 30 -30 0];
state = [1 2 2 1];
F = [1 1; 0 1];
H = [1 0; 0 1];
Q = Q_ahu{2};
R = R_ahu{2};

y_prev = rand*100;
vy_prev = 0;

Y = [];
X = [];

%Y = N(FY,Q)
%X = N(HX,R)
for i = 1:length(len)
    count = len(i);
    for t = 1:count
        tmp_y = mvnrnd( F*[y_prev vel(i)]',  Q{state(i)} );
%         tmp = mvnrnd( F*[y_prev vy_prev]',  Q{state(i)} );
        Y = [Y; tmp_y];
        tmp_x = mvnrnd( H*[y_prev vy_prev]', R{state(i)} );
        X = [X; tmp_x];
        y_prev = tmp_y(1);
        vy_prev = tmp_y(2);
    end
    
end

data = X(:,1);
[res, Z, M, Q_, R_] = gibbs_sgf(data, 1, 0);
Z = mode(Z(:,11:2:end),2);

figure
hold on
yyaxis left
%events in shade
stem(1:length(Z), (1-Z)*max(data), 'Marker','None', 'LineWidth', 5, 'Color',[.8 .8 .8]);
%original data
plot(1:length(res), data, 'k--', 'LineWidth', 1.5)
%filtered data
%     plot(1:length(res), res(:,1),'g-')
yyaxis right
%velocity
plot(1:length(res), res(:,2), 'r-','LineWidth', 2)
