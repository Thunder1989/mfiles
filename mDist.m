function distance = mDist(a, b)
    %compute the Manhattan distance of two vectors
    distance = 0;
    for i = 1:size(a,2)
        if a(i)~=b(i)
            distance = distance + 1;
        end
    end
end