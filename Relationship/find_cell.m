function [x,y] = find_cell(mat)
    metric = zeros(size(mat));
    curv = zeros(size(mat));
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            c = mat{i,j};
%             c = c - mean(c);
%             c = c/std(c);
            [v,idx] = sort(c,'descend');
            metric(i,j) =  v(1)-v(2) + v(1)-v(3);
            curv(i,j) = std(c);
        end
    end
    metric  = metric .* curv;
    [x,y] = ind2sub( size(metric), find(metric==max(metric(:)), 1) );
