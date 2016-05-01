function list = kNN(a, k)
    %a is the array containing the distances of this to other vectors
    %and return the index of the k neighbors
    list = zeros(1,k);
    len = size(a,2);
    %create the heap
    for i = 1:k
        list(i) = i;
    end
    %compare and replace
    for i = k+1:len
        max = 0;
        index = 0;
        for j = 1:k
            if a(list(j)) > max
                max = a(list(j));
                index = j;
            end
        end
        
        if a(i) < max
            list(index) = i;  
        end
    end
end