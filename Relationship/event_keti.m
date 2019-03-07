function event_eval(debug)
    
    clear
    clc
   
    p = gcp('nocreate'); % If no pool, do not create new one.
    if ~isempty(p)
        delete(gcp('nocreate'))
    end
    num = feature('numcores'); %check #of cores
    c = parcluster('local');
    %c.NumWorkers = num*2;
    c.NumWorkers = 5;
    parpool(c, c.NumWorkers);
    
    direc = '../../co-location/KETI_oneweek/'
    rooms = dir(direc);
   
    spmd
        idx = labindex + 48;
        %idx = labindex;
        sensors = dir(strcat(direc, rooms(idx).name, '/*.csv'));
        
        %res = [];
        for i = 1:length(sensors)
        
            fn = [direc, rooms(idx).name, '/', sensors(i).name] 
            cur = csvread(fn); 
            cur = cur(:,2);
            cur = cur(1:50:end);
            
            if length(unique(cur)) ~= 1 
                try
                    [~, Z, ~, ~, ~] = gibbs_sgf_K(cur, 2, 1, 0);
                catch e
                    fprintf(1, 'err was:\n%s', e.message)
                    Z = ones(size(cur,1),11);
                end
            else
                fprintf('constant values only for %s ...\n',fn);
                continue
            end

            Z_ = mode(Z(:,11:3:end),2)-1;
            Z_ = remap_event(Z_);
            %res = [res; Z];
            fprintf('done for %s ...\n',fn);
            FN = sprintf('./events/keti/%s', [rooms(idx).name,'_', sensors(i).name]);
            FID = fopen(FN,'w');
            fprintf(FID, '%d', Z_);
            fclose(FID);
       end 
    end
    
    %out = cell(length(res),1);
    %for i=1:length(res)
        %out{i} = res{i};
    %    cur = res{i};
        
    %    FN = sprintf('./events/keti/%s.txt', [rooms(idx).name, sensors(i).name] );
    %    FID = fopen(FN,'w');
    %    fprintf(FID, '%d', Z_);
    %    fclose(FID);
    %end
    %celldisp(out)
