%------------------------%
%%% Artifact Rejection %%%
%------------------------%



function [artifacts_byChan,data_bha_artf,artifact_timestamps,currentSpikes] = func_artifact_rejection_Alt2(data_norm,data_bha,goodChannels,params)
fprintf("Artifact Rejection: \n")

    
    %%% Ripple Detection Scheme: 
    data_bha_artf = data_norm;
    artifact_timestamps = cell(1,size(data_norm.trial{1},1));
    
    parfor chan = 1:size(data_norm.trial{1},1)
    
        fprintf('Processing channel %d/%d ...\n', chan, size(data_norm.trial{1},1));
    
        dat = data_norm.trial{1}(chan,:);
    
        %%% Bandpass filtering 
        dat_filt = ft_preproc_bandpassfilter(dat, 500, [20 80], [], 'fir','twopass');
    
        %%%Normalisation and power
        dat_z = (dat_filt - mean(dat_filt)) / std(dat_filt);
        nss = dat_z.^2;
        mask = nss > 3;
    
        %%% Joining gaps together (if they are less than 7)
        mask2 = imclose(mask, ones(1,7));
        CC = bwconncomp(mask2);
    
        IED_idx = [];
        IED_nums = 0;
        for c = 1:CC.NumObjects
    
            cluster = CC.PixelIdxList{c};
            L = numel(cluster);
    
            if L >=  round(50  /1000 * 500) && L <= round(250 /1000 * 500)
    
                % Findpeaks inside the cluster
                [pks, locs] = findpeaks(nss(cluster), 'MinPeakHeight', 10);
    
                if ~isempty(locs)
                    %%% Max Peak in clusters 
                    [~, maxpk] = max(pks);
                    peak_idx = cluster(locs(maxpk));
    
                    pad = params.spikeWindow;
    
                    start_idx = max(1, peak_idx - pad);
                    end_idx   = min(length(dat), peak_idx + pad);
                    IED_nums = IED_nums+1;
                    IED_idx = [IED_idx start_idx:end_idx];
                end
            end
        end
        artifact_timestamps{chan} = unique(IED_idx);
    end
    
    for chan = 1:size(data_norm.trial{1},1)
        data_bha_artf.trial{1}(chan, artifact_timestamps{chan}) = NaN;
        fprintf('Channel %d → %d samples removed\n', chan, length(artifact_timestamps{chan}));
    end
    fprintf('IED detection complete!\n');
    
    
    artifact_timestamps = artifact_timestamps';
    
    for channel = 1:length(goodChannels)
        currentSpikes{channel} = [];
    end
    
    
    
    %%
    %%% Cluster Together Arrtifacts
    artifacts_byChan = cell(1,size(data_bha.trial{1},1));
    data_bha_artifactless = data_bha;
    for chan = 1:size(data_bha.trial{1},1)
        clust_count = 1;
        artifs = [];
        artifact_clusters = {};
        data_bha_artifactless.trial{1}(chan,artifact_timestamps{chan}(artifact_timestamps{chan}>0)) = nan;
        for i = 1:length(artifact_timestamps{chan})-1
            if ~(artifact_timestamps{chan}(i+1) - artifact_timestamps{chan}(i) < 3)% Clustering when 3 timestamps long
                clust_count =clust_count +1;
                artifs = [];
            else
                continue; % XXX
            end
            artifs = [ artifs artifact_timestamps{chan}(i+1)];
            artifact_clusters{clust_count} = artifs;
    
        end
        artifacts_byChan{chan} = artifact_clusters;
    end
    
    clear artifact_clusters artifs chan clust_count
    fprintf("done! \n")

end