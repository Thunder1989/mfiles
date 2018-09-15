function frame_zcr = msf_framezcr(signal, frame_len, frame_step)  
    %signal is the input signal
    %frame_len is the # of samples per frame
    %frame_step is the # of samples shifted between two frames

    if size(signal,1) ~= 1
        signal = signal';
    end
    
    signal_len = length(signal);
    if signal_len <= frame_len  % if very short frame, pad it to frame_len
        num_frames = 1;
    else
        num_frames = 1 + ceil((signal_len - frame_len)/frame_step);
    end
    padded_len = (num_frames-1)*frame_step + frame_len;
    % make sure signal is exactly divisible into N frames
    pad_signal = [signal, zeros(1,padded_len - signal_len)];
    
    % build array of indices
    indices = repmat(1:frame_len, num_frames, 1) + ...
        repmat((0: frame_step: num_frames*frame_step-1)', 1, frame_len);
    frames = pad_signal(indices);
    
    % calculate ZCR
    frame_zcr = zeros(size(frames,1),1);
    for i = 1:size(frames,1)
        frame_zcr(i) = ZCR(frames(i,:));
    end
end