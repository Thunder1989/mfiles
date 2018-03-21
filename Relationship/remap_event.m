function reading = remap_event(input)
    %remap event sequence to make 1 for events, 0 for non-events
    input = input(:);
    
    num_zero = sum(input == 0);
    num_one = length(input) - num_zero;
    
    if num_one > num_zero
        reading = abs(1 - input)';
    else
        reading = input';
    end
    