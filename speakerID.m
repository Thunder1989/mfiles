%% pitch experiment
close all
clear
clc
f = dir('./voice/');
for i = 3:3
    [y,Fs] = audioread(f(i).name);
    N = length(y);
    cur = 1;
    nsample = round(Fs*40/1000);
    noverlap = nsample/2;
    pitch =  [];
    ts = [];
    while cur + nsample < N
        ts = [ts cur];
        y_tmp = y(cur:cur+nsample-1);
        y_f = fft(y_tmp);
        ms1 = 1;
        ms20 = 20;
% ceptrum based pitch calculation
        Y = fft(y_tmp .* hamming(length(y_tmp)));
        C = ifft(log(abs(Y)+eps));
        q = (ms1:ms20)/Fs;
%         maxlag = 20;
%         r = xcorr(y_tmp, maxlag, 'coeff');
%         r = r(floor(length(r)/2):end);
        ms2 = 2; % 2ms
        [max_amp, idx] = max(abs(C(ms2:ms20)));
%         [max_amp, idx] = max(r(ms2:ms20));

        pitch_cur = Fs / (ms2+idx-1);
        pitch = [pitch pitch_cur];
        cur = cur + nsample - noverlap;
    end
    figure
    hold on
    yyaxis left
    plot(y)
    yyaxis right
    plot(ts, pitch, 'ro', 'MarkerSize', 4);
    
%     b = fir1(48,[400/Fs 4000/Fs]);
%     y = filter(b,1,y);
%     figure
%     spectrogram(y,256,128,1024,Fs,'yaxis')
%     ylim([0.4, 4])
end

%% mfcc works
close all
clear
clc
files = dir('./voice/*christian_ondesk*');
X = [];
for i = 1:length(files)
    fname = strcat('./voice/', files(i).name)
    str = strsplit(files(i).name,'_');
    if strcmp(str{1},'christian')
        class = 1;
    elseif strcmp(str{1},'dezhi')
        class = 2;
    else
        class = 3;
    end
    [y,Fs] = audioread(fname);
% b = fir1(48,[200/Fs 4000/Fs]); %filtering doesn't make much diff
% y = filter(b,1,y);
    f = get_acousic_features(y,Fs,1);
    f = [f class * ones(length(f),1)];
% calculate ZCR
% cur = 1;
% nsample = ceil(Fs*25/1000);
% nshift = ceil(Fs*10/1000);
% zcr =  [];
% ts = [];
% N = length(y);
% while cur + nsample < N
%     ts = [ts cur];
%     y_tmp = y(cur:cur+nsample-1);
%     zcr = [zcr ZCR(y_tmp)];
%     cur = cur + nshift-1;
% end

% idx = 1:ceil(10/1000*Fs):length(y);
% figure
% hold on
% yyaxis left
% plot(y)
% yyaxis right
% plot(idx(1:length(mfccs)), mfccs(:,1),'ro', 'MarkerSize', 4);
% plot(idx(1:length(mfccs)), tao,'k--', 'LineWidth', 2);
% plot(ts, zcr, 'ro', 'MarkerSize', 4);
    X = [X; f];
end

csvwrite('v_f.csv', X);

% options = statset('MaxIter',1000);
% gmfit = fitgmdist(X(:,2:end),3,'CovarianceType','full',...
%             'SharedCovariance',false,'Options',options);
% clusterX = cluster(gmfit,X(:,2:end));
% 
% coef = pca(X(:,2:end));
% X_ = X(:,2:end) * coef(:,1:2);
% cluster_km = kmeans(X_,3);

%% plot
res = [0.512396694215 0.549311294766 0.695456814778 0.912131802297 0.456783056669 0.511734401832;
    0.463483146067 0.714606741573 0.594381468704 0.9122720552 0.431246655966 0.477260567148;
    0.481142241379 0.623922413793 0.665474974464 0.917773237998 0.617009750813 0.732394366197];
h = bar(res,'grouped');
set(h,'BarWidth',0.9); % The bars will now touch each other
set(gca,'YGrid','on')
%set(gca,'GridLineStyle','+')
%set(gca,'XTicklabel','LW>=sD&&LH>=sD|LW>=sD&&LH<=sD|LW<=sD&&LH>=sD|LW<=sD&&LH<=sD')
set(get(gca,'YLabel'),'String','SpeakerID Acc')
lh = legend('Christian-before','Christian-after','Dezhi-before','Dezhi-after','Vijay-before','Vijay-after' );
set(lh,'Location','SouthEast','Orientation','vertical')
colormap(summer) 