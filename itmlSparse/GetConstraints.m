function C = GetConstraints(X, y, num_constraints, l, u)
% C = GetConstraints(y, num_constraints, l, u)
%
% Get ITML constraint matrix from true labels.  See ItmlAlg.m for
% description of the constraint matrix format

consArray = [];
C = zeros(num_constraints, 4);
a = length(X);
for m=1:size(X,1)
    for n=m+1:size(X,1)
        if y(m)==y(n)
            consArray=[consArray ; [m n 1 l]];
        else
            consArray=[consArray ; [m n -1 u]];
        end
    end
end
          
for k=1:num_constraints,
    m = length(consArray);
    i = ceil(rand * m);
    C(k,:) = consArray(i,:);
end


