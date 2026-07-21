

function [artifacts_byChan,data_bha_artf,artifact_timestamps,currentSpikes] = func_artifact_rejection_Alt3(data_norm,data_bha,goodChannels,params)
    
    
        
    fs = data_norm.fsample;
    dat = data_norm.trial{1};
    Nchan = size(dat,1);
    
    data_artf = data_norm;
    rejection_timestamps = cell(Nchan,1);
    
    pad = round(0.25 * fs);     % ±250 ms
    min_clean = round(3 * fs);  % 3 s
    rms_win = round(0.1 * fs);  % 100 ms
    allPaddedIdx = {};
    parfor chan = 1:Nchan
    
        x = dat(chan,:);
    
        % 1) 0.3–150 Hz filter
        x_filt = ft_preproc_bandpassfilter(x, fs, [0.3 150], [], 'fir','twopass');
    
        art_idx = [];
    
        % 2) Amplitude artifacts ±750 µV
        amp_idx = find(abs(x_filt) > 750);
    
        % 3) Gradient artifacts (median ± 6*IQR)
        dx = diff(x_filt);
        med_dx = median(dx,'omitnan');
        iqr_dx = iqr(dx);
    
        thr_up = med_dx + 6*iqr_dx;
        thr_low = med_dx - 6*iqr_dx;
    
        grad_idx = find(dx > thr_up | dx < thr_low);
    
        % 4) High-frequency bursts (>150 Hz)
        x_hp = ft_preproc_highpassfilter(x, fs, 150, [], 'fir','twopass');
    
        rms_sig = sqrt(movmean(x_hp.^2, rms_win,'omitnan'));
    
        med_rms = median(rms_sig,'omitnan');
        iqr_rms = iqr(rms_sig);
    
        thr_rms = med_rms + 4*iqr_rms;
    
        hf_idx = find(rms_sig > thr_rms);
    
        % Merge all artifact indices
        art_idx = unique([amp_idx grad_idx hf_idx]);
    
        % Padding ±250 ms
        padded_idx = [];    
        for i = 1:length(art_idx)
            s = max(1, art_idx(i)-pad);
            e = min(length(x), art_idx(i)+pad);
            padded_idx = [padded_idx s:e];
        end
    
        padded_idx = unique(padded_idx);
        allPaddedIdx{chan} = padded_idx
    
    end
    
    for chan = 1:Nchan
         data_artf.trial{1}(chan,allPaddedIdx{chan}) = NaN;
    
        % Remove short clean intervals (<3 s)
        clean = setdiff(1:length(dat(chan,:)), allPaddedIdx{chan});
    
        if ~isempty(clean)
            diff_clean = diff(clean);
            breaks = [0 find(diff_clean>1) length(clean)];
            
            for b = 1:length(breaks)-1
                segment = clean(breaks(b)+1 : breaks(b+1));
                if length(segment) < min_clean
                    data_artf.trial{1}(chan,segment) = NaN;
                    allPaddedIdx{chan} = unique([allPaddedIdx{chan} segment]);
                end
            end
        end
    
        rejection_timestamps{chan} = allPaddedIdx{chan};
    end

    
    

    
    %%% Spike Detection
    fprintf("Spike Detection ...")
    clear ts
    
    ts.rawTs= zscore(data_norm.trial{1}',[],1);
    ts.Fs = data_bha.fsample;
    ts.chanNames = data_bha.label;
    
    %%% FindSpikes - Function for our SPIKE detection:
    mySingleton = findSpikeTimes_reviewed(ts,params.spikePeakWin,...
        params.spikeZThresh,params.spikeCtsThresh, params.spikeWindow, ...
        'amp_scale',params.spikeAmpScale,'maxNegPeakWidth',...
        params.spikeMNegPeakW,'trackPeaks',params.spikeTrackPeaks);
    
    currentSpikes = cell(1,length(mySingleton.spikeTime));
    for channel = 1:length(mySingleton.spikeTime)
        currentSpikes{channel} = mySingleton.spikeTime{channel};
    end
    
    clear mySingleton ts segm trials spikeTimes tr_count
    fprintf("done! \n")
    
    for channel = 1:length(goodChannels)
        if ~isempty(currentSpikes{channel})
            artifact_timestamps{channel} =[artifact_timestamps{channel} currentSpikes{channel}(1):currentSpikes{channel}(2)];
            artifact_timestamps{channel} = unique(artifact_timestamps{channel});
        end
    end
    
    
    %%% Cluster Together Arrtifacts
    artifacts_byChan = cell(1,size(data_bha.trial{1},1));
    data_bha_artifactless = data_bha;
    for chan = 1:size(data_bha.trial{1},1)
        clust_count = 1;
        artifs = [];
        artifact_clusters = {};
        data_bha_artifactless.trial{1}(chan,artifact_timestamps{chan}(artifact_timestamps{chan}>0)) = nan;
        for i = 1:length(artifact_timestamps{chan})-1
            if (artifact_timestamps{chan}(i+1) - artifact_timestamps{chan}(i) < 3)% Clustering when 3 timestamps long
                
            else
                clust_count =clust_count +1;
                artifs = [];
            end
            artifs = [ artifs artifact_timestamps{chan}(i+1)];
            artifact_clusters{clust_count} = artifs;
    
        end
        artifacts_byChan{chan} = artifact_clusters;
    end
    
    clear artifact_clusters artifs chan clust_count
    fprintf("done! \n")
end