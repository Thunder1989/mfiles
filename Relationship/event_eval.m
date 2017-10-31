close all
clear
clc

path = 'D:\HouseA\';
files = dir(strcat(path, '*.log'));
colors = containers.Map( 1:3, {[54/255,160/255,204/255], [211/255,142/255,194/255], [.8 .8 .455]} );

%% gibbs sampling based
% close all
num = length(files);
for n = 1:num

    fn = [path, files(n).name];
    cur = importdata(fn); %skip the 1st row, which is the headers
    cur = cur(1:3000,:);
    cur = cellfun(@strsplit, cur, 'UniformOutput', false);
	type = regexp(files(n).name,'[A-Z]+','match', 'once');
    if strcmp(type, 'PIR')
        continue;
%         cur = cellfun(@(x) x(end-1), cur, 'UniformOutput', false);
%         cur = cell2mat( cellfun(@(x) strcmp(x, 'True'), cur, 'UniformOutput', false) );
    else
        cur = cell2mat( cellfun(@(x) str2double( x(end-1) ), cur, 'UniformOutput', false) );
    end
    [~, Z, M, ~, ~] = gibbs_sgf_K(cur, 2, 1, 0);
    Z_ = mode(Z(:,11:3:end),2);

    figure
    hold on
    y_lim = max(cur);
    for i=1:length(Z_)
        stem(i, y_lim, 'Marker','None', 'LineWidth', 4, 'Color', colors(Z_(i)) );
    end
    plot(cur,'k')
    pause;

end