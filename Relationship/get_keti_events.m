function [event, canny] = get_keti_events()
    close all
    clc

%     path = 'D:\TraneData\KETI_oneweek\';
%     rooms = dir(path);
% 
%     event = cell((length(rooms)-2)*4, 1);
%     canny = cell((length(rooms)-2)*4, 1);
%     for n = 3:length(rooms)
%         fprintf('processing %s ...\n',rooms(n).name);
%         files = dir(strcat(path, rooms(n).name, '\*.csv'));
%         ctr = 1;
%         for m = 1:length(files)
%             if m==4
%                 continue %skip pir for now
%             end
%             fn = strcat(path, rooms(n).name, '\', files(m).name);
%             cur_data = csvread(fn,1); %skip the 1st row, which is the headers
%             cur_data = cur_data(:,2);
%             cur_data = cur_data(1:30:end); %downsample
% %             plot(cur_data)
% %             pause
%             try
%                 [res,Z,M] = gibbs_sgf(cur_data,1,0);
%             catch
%                 warning('GloME exception on %s ...\n',rooms(n).name);
%                 Z = zeros(size(cur_data,1),11);
%             end
% %             event{(n-3)*4+ctr} = get_KF_residual(cur_data');
%             event{(n-3)*4+ctr} = mean(Z(:,11:3:end),2);
%             canny{(n-3)*4+ctr} = get_canny_edges(cur_data');
%             ctr = ctr + 1;
%         end
%     end

    %use pre-processed matrix
    load('keti_processed_ssc.mat');
    data = input_data;
    data(4:5:end,:) = []; %remove pir since mostly are zeros
    data = data(:, 1:2:end); %downsample

    event = cell( size(data,1), 1 );
    canny = cell( size(data,1), 1 );
    for i = 1:size(data,1)
        fprintf('processing %d-th sensor...\n',i);
        cur_data = data(i,:);
        try
            [res,Z,M] = gibbs_sgf(cur_data,1,0);
        catch
            warning('GloME exception on %s-th sensor ...\n',i);
            Z = zeros(size(cur_data,1),11);
        end
%             event{(n-3)*4+ctr} = get_KF_residual(cur_data');
        event{i} = mean(Z(:,11:3:end),2);
        canny{i} = get_canny_edges(cur_data);
    end  
    
function res = get_KF_residual(X) %X is a row input vector
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
