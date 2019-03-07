tic
n = 200;
A = 500;
a = zeros(n,1);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
toc
