prompt = 'What is the original value? ';
res = input(prompt,'s');
if strcmp(res, 'yes') || strcmp(res, 'y')
    fprintf('We got a yes\n')
elseif strcmp(res, 'no') || strcmp(res, 'n')
    fprintf('We got a no\n')
else
    fprintf('We got a %s\n', res)
end