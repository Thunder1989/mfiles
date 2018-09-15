function get_keti_acc(event, canny)

    path = 'D:\TraneData\KETI_oneweek\';
    rooms = dir(path);

%     event = cellfun(@transpose,event,'UniformOutput',false);
%     max_len = max(cell2mat(cellfun(@length,event,'UniformOutput',false)));
%     empty_idx = find( cellfun(@isempty,event) );
%     for i = 1:length(empty_idx)
%         event{empty_idx(i)} = zeros(1,max_len);
%     end
%     event{175} = zeros(1,max_len);
%     min_len = min(cell2mat(cellfun(@length,event,'UniformOutput',false)));
%     event = cellfun(@(x) x(1:min_len),event,'UniformOutput',false);
%     event = cellfun(@(x) resample(x, min_len, length(x)),event,'UniformOutput',false);
%     event = tfidf(cell2mat(event));
%     event = mat2cell(event,ones(1,size(event,1)),[size(event,2)]);
%     disp(event);
    
    ref = 2;
    ctr = 0;
    ctr1 = 0;
    for m = 1:length(event)
        if mod(m,4)==ref
            continue
        end
        event_cur = event{m};
        edge_cur = canny{m};

        sims = zeros(length(rooms)-2,1);
        scores = zeros(length(rooms)-2,1);
        for n = ref:4:length(event)
            event_ref = event{n};
            edge_ref = canny{n};

            if sum(edge_ref)==0 || sum(edge_cur)==0
                sims(ceil(n/4)) = 0;
            else
                [edge_ref, edge_cur] = align_vector(edge_ref, edge_cur);
                cur_sim = dot(edge_ref, edge_cur)/(norm(edge_ref)*norm(edge_cur));
                sims(ceil(n/4)) = cur_sim;
            end
            
            if sum(event_ref)==0 || sum(event_cur)==0
                scores(n) = 0;
            else
                [event_ref, event_cur] = align_vector(event_ref, event_cur);
                scores(n) = dot(event_ref, event_cur) / ( norm(event_ref) * norm(event_cur) );
            end
        end

        true = ceil(m/4);
        if ismember(true, find(sims==max(sims))) && length( find(sims==max(sims)) ) < length(sims);
            ctr = ctr + 1;
        else
            flag = ismember( true, find(scores==max(scores)) ) & length( find(scores==max(scores)) ) < length(scores);
            if flag == 1
                ctr1 = ctr1 + 1;
            end
        end
    end

    fprintf('acc on canny edge is %.4f\n', ctr/((length(rooms)-2)*3));
    fprintf('acc on combination is %.4f\n', (ctr+ctr1)/((length(rooms)-2)*3));

function [a_, b_] = align_vector(a,b)
    len = min(length(a), length(b));
    a_ = a(1:len);
    b_ = b(1:len);
    