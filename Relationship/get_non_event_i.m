function idx = get_non_event_i(R)

    R_min = min(cell2mat(cellfun(@(x) x(1), R, 'UniformOutput', false))); 
   
    for i=1:length(R)
        if R{i}(1) == R_min 
            idx = i;
        end
    end
