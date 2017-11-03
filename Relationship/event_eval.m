function event_eval(debug)
    
%     close all
    clear
    clc

    path = 'D:\HouseA\';
    files = dir(strcat(path, '*.log'));
    colors = containers.Map( 1:3, {[54/255,160/255,204/255], [211/255,142/255,194/255], [.8 .8 .455]} );

    num = length(files);
    step = 3000;
    for n = 2:2

        fn = [path, files(n).name];
        type = regexp(files(n).name,'[A-Z]+','match', 'once');
        if strcmp(type, 'PIR')
            continue;
        %         cur = cellfun(@(x) x(end-1), cur, 'UniformOutput', false);
        %         cur = cell2mat( cellfun(@(x) strcmp(x, 'True'), cur, 'UniformOutput', false) );
        end
        
        cur = importdata(fn); %skip the 1st row, which is the headers
        FN = sprintf('%s_output.txt', files(n).name);
        FID = fopen(FN,'w');
        
        i = 1;
        len = length(cur);
        while i<len
            ending = min(i+step-1, len);
            tmp = cur(i:ending, :);
            tmp = cellfun(@strsplit, tmp, 'UniformOutput', false);
            tmp = cell2mat( cellfun(@(x) str2double( x(end-1) ), tmp, 'UniformOutput', false) );

            [~, Z, ~, ~, ~] = gibbs_sgf_K(tmp, 2, 1, 0);
            Z_ = mode(Z(:,11:3:end),2);
            fprintf(FID, '%d', Z_);
            i = i + 2999;
            fprintf('done for peroid starting from %d ...\n',i);
        end

        fclose(FID);
        
        if debug==1
            figure
            hold on
            y_lim = max(cur);
            for i=1:length(Z_)
                stem(i, y_lim, 'Marker','None', 'LineWidth', 4, 'Color', colors(Z_(i)) );
            end
            plot(cur,'k');
        end
        
    end