function M = get_M(Z_sample,K) %K-class transition matrix, each col sums to 1, checked
    
    M = zeros(K); % # of transitions i-j
    N = zeros(K,1); % # of i
    for n = 1:size(Z_sample,2)
        Z = Z_sample(:,n);
        for i = 1:K
            N(i) = N(i) + length( find(Z(1:end-1)==i) );
            for j = 1:K
                N_ij = length(find(Z(1:end-1)==i & Z(2:end)==j)); 
                M(i,j) = M(i,j) + N_ij;
            end
        end
    end
    
    assert (sum(N) == numel(Z_sample) - size(Z_sample,2));

    M = bsxfun(@rdivide,M,N)';