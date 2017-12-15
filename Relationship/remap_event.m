function reading = remap_event(input)

    input = input(:);
    
    num_zero = sum(input == 0);
    num_one = length(input) - num_zero;
    
    if num_one > num_zero
        reading = abs(1 - input)';
    else
        reading = input';
    end
    