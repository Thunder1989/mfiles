%% re-arrange the cell
clear
clc
input = load('/Users/hdz_1989/Documents/Dropbox/scores.mat');
scores = input.scores;
score = {};
id = {};
k = 1;
% for i=1:20
%     for j=1:i-1
%         %remove mean+normalize
%         tmp = scores{i,j};
%         tmp = zscore(tmp);
% %         tmp = tmp-mean(tmp);
%         tmp = tmp/max(abs(tmp));
%         score{k} = tmp;
%         id{k} = [i,j];
%         k=k+1;
%     end
% end

%get pairs in the same room
for i=1:5
    sub = scores(4*i-3:4*i,4*i-3:4*i);
    for j=1:4
        for l=1:j-1
            %remove mean+normalize
            tmp = sub{j,l};
            tmp = zscore(tmp);
%             tmp = tmp-mean(tmp);
            tmp = tmp/max(abs(tmp));
            score{k} = tmp;
            id{k} = [i,j];
            k=k+1;
        end
    end
end

% % get all the H_T pairs
% for i=1:20
%     for j=1:i-1
%         %remove mean+normalize
%             if (mod(i,4)==0&&mod(j,4)==2) || (mod(i,4)==2&&mod(j,4)==0)
%             tmp = scores{i,j};
%             tmp = tmp-mean(tmp);
%             tmp = tmp/max(abs(tmp));
%             score{k} = tmp;
%             id{k} = [i,j];
%             k=k+1;
%         end
%     end
% end

% score = scores(3:end,2);

id = id';
mat = zeros(200,length(score));
for i=1:200
    mat(i,:) = cellfun(@(x) x(i), score);
end
% HeatMap((abs(diff(mat)))')
% HeatMap((diff(mat)))
residue = mat(2:end,:)-diff(mat); 
HeatMap(residue)

%% brute force loop to check the info in the data set
% TBD
    
    
    