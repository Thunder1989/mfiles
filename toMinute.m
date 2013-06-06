function minute = toMinute(a)
    minute = floor(a)*60 + mod(a,1)*100;
end