

function [artifacts_byChan,data_bha_artifactless,artifact_timestamps,currentSpikes] = MW_3_1_ArtifactRejection_Alt1(data_norm,data_bha,goodChannels,params)
    
    %  Rejects artifacts on the basis of Range and IQR: 
    %    1) IQR:   Segments and then rejects given segment if a mean of a segment
    %              is more than Q3 and 2.3 of IQR of a given Data
    %    2) RANGE: Filtering data with High Pass 225, rejecting data points
    %              exceeding 5 stds of a whole data; Complememmntary computing
    %              differentaition, taking absolute, and rejecting data points
    %              exceedin 5 stds
    %
    %  Later, each rejected point is clustered and +/- 100 samples are added to Padd each artifact 
    %
    % Returns:
    %       - artifacts:     raw data from range artifacts, 
    %                               detected as artifactual, without padding
    %       - artifact_timestamps: a list of rejected timestamps including
    %                               IQR
    %       - data_bha_artf:        Data in which every padded artifact 
    %                               have been replaces by NaNs 

    
    
    clear Seltrls
    clear Seltrls_all
    clear artifacts
    
    rejMode = "precise"; % Two modes of Rejection:
                         %  + precise: 10 timepoint around artefacts
                         %  + approx: 200 ms around peaks Work in progr!
    
    clear dat_bha dat_norm
    
    dat_norm = data_norm.trial{1};
    Nchan = size(dat_norm,1);
    
    w = 2.3; % Iqr thresholds
    data_bha_artf = data_norm;
    
    
    
    for chan = 1:Nchan
        artif_estim = [];
        temp_artif = [];
        dat = dat_norm(chan,:);
        clear dat_trl_mean dat_trl_std dat_trl_raw trl2large trl2small
        
        
        %------------------%
        % 1: IQR rejection %
        %------------------%
        fprintf('    Removing IQR Artifacts ... ')
    
        % Reject those trials
        iqr_artifacts = [];
        data_iqr = iqr(dat);
        q3 = quantile(dat,0.75);
        %sham segmentation for
        for begin = 1:250:length(dat)
            if begin +250 < length(dat)
                current_segment_time = begin:begin+250;
            else 
                current_segment_time = begin:length(dat);
            end
            current_segment = dat(current_segment_time);
            curr_mean = nanmean(current_segment);
            if abs(curr_mean) < abs(q3) + w*data_iqr
                continue
            else
                iqr_artifacts = [iqr_artifacts current_segment_time];
            end
        end
        fprintf("done! \n")
    
        %------------------%
        % 2:  Range %
        %------------------%
        fprintf('    Removing Range Artifacts ... ')
    
    
        %note - probably it would be better to resample and then filter to
        %250 HighPass Filter
        filtered_range = ft_preproc_highpassfilter(dat, data_norm.fsample,225,[],'fir','twopass');
    
        dat_diff= zscore(diff(squeeze(dat)));
        %dat_amp = zscore(squeeze(dat)); % For Staresina Purposes
        dat_filtered = zscore(squeeze(filtered_range));
        artif_estim = sort(unique([find(abs(dat_diff) >5) find(abs(dat_filtered) >5)]));
        temp_artif = [artif_estim(1)];
        artifacts_env = {};
        count = 1;  
    
        %Clustering Artifacts 
        for i = 1:length(artif_estim)-1
        
            if artif_estim(i+1) - artif_estim(i) < 2
                continue
            else
                temp_artif = [temp_artif artif_estim(i)];
            end
            artifacts_env{count} = temp_artif; 
            temp_artif = [artif_estim(i+1)];
            count = count +1;
        end
        fprintf("done! \n")
        fprintf('    Rejecting Artifacts ... ')
        temp_rejection = [];
    
        artifact_timestamps{chan} = [];
    
        %%% Get artifact Timestamps by different means:
        %    + approx - for every Artifacr Envelope get -200 and + 200 around Mean timestamp of artifact 
        %    + precise - for every artefat timestsamp get envelope of 100 timestamps around 
    
        switch rejMode
            case "approx"
                for i = artifacts_env
    
                   for j = mean([i{1}(1),i{1}(2)])-200: mean([i{1}(1),i{1}(2)])+200
                        artifact_timestamps{chan} = [artifact_timestamps{chan} j];
                   end
                end
            case "precise" % Get Artifact 
                 parfor i = 1:length(artifacts_env)
                    if artifacts_env{i}(1) < 100 % Fringe situation of artefact being in the beginning
                        addit1 = 0;
                    else
                        addit1 = 100;
                    end
    
                    if artifacts_env{i}(2) > length(data_bha_artf.trial{1}) -100 % Fringe situation of artefact being at an end
                        addit2 = 0;
                    else
                        addit2 = 100;
                    end
    
                    %Creating Padding for Artifact 
                    for j = artifacts_env{i}(1)-addit1:artifacts_env{i}(2)+addit2
                        temp_rejection = [temp_rejection j];
                    end
                end
        end
        artifact_timestamps{chan} = temp_rejection;
   
        artifact_timestamps{chan} = unique([artifact_timestamps{chan} iqr_artifacts]);
        nonArtifacts = setdiff([1:size(data_bha_artf.trial{1},2)],artifact_timestamps{chan});
    
        artifacts = artif_estim;
        data_bha_artf.trial{1}(chan,artifact_timestamps{chan}) = NaN;
        fprintf("done! \n")
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
        for channel = 1:length(goodChannels)
            if ~isempty(currentSpikes{channel})
                
                allSpikePoints = []; 
                for i = 1:size(currentSpikes{channel},1)
                    allSpikePoints = [allSpikePoints, currentSpikes{channel}(i,1):currentSpikes{channel}(i,2)];
                end
                allSpikePoints = allSpikePoints(allSpikePoints >= 1 & allSpikePoints <= length(data_norm.trial{1}(channel,:))); %%%% CHANGE                         !!!!!! CHANGE!!!!
        
        
                artifact_timestamps{channel} =[artifact_timestamps{channel} unique(allSpikePoints)];
                artifact_timestamps{channel} = unique(artifact_timestamps{channel});
            end
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
            if (artifact_timestamps{chan}(i+1) - artifact_timestamps{chan}(i) < 3)% Clustering when 3 timestamps long period
                
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
