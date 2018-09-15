%% get features function
function results = get_acousic_features(y,Fs,delta)
mfccs = msf_mfcc(y,Fs,'winlen',ceil(25/1000*Fs)/Fs,'nfilt',26,'ncep',20,'appendenergy',true);
mfccs(isnan(mfccs)| isinf(mfccs)) = -10;
if delta == 1
    mfccs_ = [zeros(2, size(mfccs,2)); mfccs; zeros(2, size(mfccs,2))];
    d_mfcc = mfccs_(4:end-1,:) - mfccs_(2:end-3,:);
    d_mfcc = d_mfcc + 2*(mfccs_(5:end,:) - mfccs_(1:end-4,:));
    d_mfcc = d_mfcc/10;
    mfccs = [mfccs d_mfcc];
end
power = mfccs(:,1);
win_len = 80; %reasonably longer window is better
p_smoothed = EWMA(power, win_len);

tao = zeros(size(p_smoothed));%threshold needs finer tuning, TBD
for i=win_len:length(p_smoothed)%also, tao per point or same tao per window?
    m = min(p_smoothed(i-win_len+1:i));
    a = mean(p_smoothed(i-win_len+1:i));
    tao(i-win_len+1:i) = repmat(m+0.7*(a-m),1,win_len);
end

% idx = 1:ceil(10/1000*Fs):length(y);
% figure
% hold on
% yyaxis left
% plot(y)
% yyaxis right
% % plot(idx(1:length(mfccs)), mfccs(:,1),'ro', 'MarkerSize', 4);
% % plot(idx(1:length(mfccs)), tao,'k--', 'LineWidth', 2);
% plot(idx(1:length(zcr)), zcr,'ro', 'MarkerSize', 4);

zcr = msf_framezcr(y,ceil(25/1000*Fs),10/1000*Fs);
mfccs = [mfccs zcr];
results = mfccs(p_smoothed>tao,:);