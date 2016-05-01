disp('Loading data');
%X = load('itmlSparse/data/irisSparse.mtx');
%y = load('itmlSparse/data/irisSparse.truth');

%rawX = csvread('20jsonX.csv');
%rawy = csvread('20jsonY.csv');
%X = spconvert(rawX);
%y = spconvert(rawy);

input = csvread('/Users/hdz/Downloads/code/sdb/rice_45min_forsdh');
X = input(:,1:end-1);
y = input(:,end);
disp('Running ITML to learn matrix A..');

num_folds = 3;
knn_neighbor_size = 4;
acc = CrossValidateKNN(y, X, @(y,X) MetricLearningAutotuneKnn(@ItmlAlg, y, X), num_folds, knn_neighbor_size);

disp(sprintf('kNN cross-validated accuracy = '));
disp(acc);


%% compute M-distance matrix
clear;
clc;
input = csvread('/Users/hdz/Downloads/code/sdb/rice_45min_forsdh');
X = input(:,1:end-1);
y = input(:,end);
A = load('matrixA.dat');
dist = zeros(size(X,1), size(X,1));
for i=1:size(X,1)
    for j=1:i
        D = X(i,:)-X(j,:);
        dist(i,j) = sqrt(D * A * D');
    end
end

%% generate CDF
same = [];
diff = [];
for i=1:size(X,1)
    for j=1:i-1
        if y(i) == y(j)
            same = [same; dist(i,j)];
        else
            diff = [diff; dist(i,j)];
        end
    end
end

figure
hold on
[p, v] = ecdf(same);
plot(v,p,'r--', 'LineWidth', 2)
[p, v] = ecdf(diff);
plot(v,p,'k--', 'LineWidth', 2)
