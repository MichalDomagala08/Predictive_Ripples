
function artifact_timestamps = ngo_spike(dat,fsample)
    fsample = 500;
    dat = dd(5,:);
    % x: 1 x nSamples, jeden kanał z całego skonkatenowanego sygnału
    % Fs: sampling rate
        
    % Optional but consistent with Ngo et al.: 0.3–150 Hz signal
    xBP = ft_preproc_bandpassfilter(dat,fsample, [0.3 150], ...
        [], 'fir', 'twopass');
    
    % Gradient / first derivative
    dx = diff(xBP);
    
    % Robust, channel-specific threshold
    gradThr = median(dx, 'omitnan') + 4 * iqr(dx);
    gradMask = abs(dx - median(dx, 'omitnan')) > gradThr;
    
    % diff has one fewer sample: map gradient hits to original-signal indices
    gradSamples = find(gradMask) + 1;
    
    % Build a time mask and pad each detection by ±250 ms
    padSamples = round(0.250 *fsample);
    
    iedMask = false(1, numel(dat));
    for k = 1:numel(gradSamples)
        idx = max(1, gradSamples(k) - padSamples) : ...
              min(numel(dat), gradSamples(k) + padSamples);
        iedMask(idx) = true;
    end
    
    artifact_timestamps{chan} = find(iedMask);


end