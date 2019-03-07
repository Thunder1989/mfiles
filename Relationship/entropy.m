function e = entropy(X)
 X = X(:)/sum(X(:));
 e = sum( X .* log2(X) );