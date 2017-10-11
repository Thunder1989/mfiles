function [event, canny] = get_keti_events()
    close all
    clc

    path = 'D:\TraneData\KETI_oneweek\';
    rooms = dir(path);

    event = cell((length(rooms)-2)*4, 1);
    canny = cell((length(rooms)-2)*4, 1);
    for n = 3:length(rooms)
        fprintf('processing %s ...\n',rooms(n).name);
        files = dir(strcat(path, rooms(n).name, '\*.csv'));
        ctr = 1;
        for m = 1:length(files)
            if m==4
                continue %skip pir for now
            end
            fn = strcat(path, rooms(n).name, '\', files(m).name);
            cur_data = csvread(fn,1); %skip the 1st row, which is the headers
            cur_data = cur_data(:,2);
            cur_data = cur_data(1:30:end); %downsample
%             plot(cur_data)
%             pause
            [res,Z,M] = gibbs_sgf(cur_data,1,0);
%             event{(n-3)*4+ctr} = get_KF_residue(cur_data');
            event{(n-3)*4+ctr} = mean(Z(:,21:3:end),2);
            canny{(n-3)*4+ctr} = get_canny_edges(cur_data');
            ctr = ctr + 1;
        end
    end

function res = get_KF_residue(X) %X is a row input vector
    N = 5*2; %Kalman Filter lookback window size
    M = 8; %EWMA window size
    kf = StandardKalmanFilter(X,M,N,'EWMA'); 
    res = abs(X - kf);
    res(isnan(res)) = 0;
    
function event = get_canny_edges(X) %X is a row input vector
    th = 1.25;    
    e = edge(repmat(X,3,1), th);
    e = e(1,:);
    e = e | [false e(1:end-1)] | [e(2:end) false];
    event = double(e);
