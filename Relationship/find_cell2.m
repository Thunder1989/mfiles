function [x,y] = find_cell2(mat)
    metric = zeros(size(mat));
    curv = zeros(size(mat));
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            c = mat{i,j};
            [v,idx] = sort(c,'descend');
            metric(i,j) =  1 - normpdf( v(1), mean(v(2:end)), std(v(2:end)) );
            if isnan(metric(i,j))
                metric(i,j) = 0;
            end
            curv(i,j) = std(v);
        end
    end
    metric = metric .* curv;
    [x,y] = ind2sub( size(metric), find(metric==max(metric(:)), 1) );
