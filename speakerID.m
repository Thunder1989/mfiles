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

%%
f = dir('./voice/');
[y,Fs] = audioread(strcat('./voice/', f(5).name));
[f0_time,f0_value,SHR,f0_candidates]=shrp(y,Fs);
figure
hold on
yyaxis left
plot(y)
yyaxis right
plot(f0_time/1000'*Fs,f0_value', 'ro', 'MarkerSize', 4);
