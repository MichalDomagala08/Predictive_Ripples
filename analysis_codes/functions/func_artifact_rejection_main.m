

function [artifacts_byChan,ieds_byChan,iqr_byChan,range_byChan,data_artifact,ied_timestamps] = func_artifact_rejection_main(data_norm,params)
    
    %  Rejects artifacts on the basis of Range and IQR: 
    %    0) SPIKE Detection: 
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



    %%% -- PREPARATION -- %%%
    clear Seltrls Seltrls_all artifacts dat_bha dat_norm

    % Check if signal is a single trial or more (if so, concat them together)
    if size(data_norm.trial,1) == 1
        singleTrialFlag = 1; else; singleTrialFlag = 0; end

    dat_norm = data_norm.trial;
    norm_concat = cell2mat(dat_norm);

    pars = params; pars.lowpassfreq = 70; pars.highpassfreq = 150;
    bha_sign = func_rippleband_filtering(data_norm,pars);
    bha_concat = cell2mat(bha_sign.trial);
    data_artf = norm_concat;

    Nchan = size(dat_norm{1},1);
    trialBoundaries = 0:7001:size(norm_concat,2);
    trialBoundaries(trialBoundaries == 0) = []; % Usuwamy zero, jeśli nie chcemy maskować początku
    boundary_samples = unique(bsxfun(@plus, trialBoundaries', -5:5));  % Tworzymy wektor relatywny [-5:5] i dodajemy go do każdego punktu granicznego

    % Parameters
    artPadding = params.artPadding; % how many samples should we pad the artifacts with
    cluster_tolerance = params.cluster_tolerance;  % How many samples is it "close" 
    range_threshold = params.range_threshold; % How many stds from mean is considered an artifact 
    w  = params.iqr_w;

    %---------------------%
    % 1:  Spike Detection %
    %---------------------%

    fprintf("(1) Spike Detection ...\n")
    clear ts
    
    ts.rawTs= zscore(norm_concat',[],1);
    ts.Fs = data_norm.fsample;
    ts.chanNames = data_norm.label;
    
    
    mySingleton = findSpikeTimes_reviewed(ts,params.spikePeakWin,...
        params.spikeZThresh,params.spikeCtsThresh, params.spikeWindow, ...
        'amp_scale',params.spikeAmpScale,'maxNegPeakWidth',...
        params.spikeMNegPeakW,'trackPeaks',params.spikeTrackPeaks);
    
    currentSpikes = cell(1,length(mySingleton.spikeTime));
    for channel = 1:length(mySingleton.spikeTime)
        currentSpikes{channel} = mySingleton.spikeTime{channel};
        tmp = arrayfun(@(s,e) max(1,s):min(length( ts.rawTs),e), ...
            currentSpikes{channel}(:,1), currentSpikes{channel}(:,2), 'UniformOutput', false);
        ied_timestamps{channel} = unique([tmp{:}]);
        fprintf('      %s: %d spikes\n', data_norm.label{channel}, length(currentSpikes{channel}));
    end
    clear mySingleton ts segm trials spikeTimes tr_count




    %----------------------------%
    % 1: Channel (IQR) rejection %
    %----------------------------%
    fprintf("\n(2) Bad Channel IQR detection ...")

    trialIQR = cell2mat(cellfun(@(x) iqr(x,2), dat_norm, 'UniformOutput',false)); % channels × trials
    chanIQR  = median(trialIQR,2);                                                        % 1 IQR / kanał
    badChans  = find(chanIQR > median(chanIQR) + 6*mad(chanIQR,1)); % 5 times median deviation
    

    fprintf('  Noisy Channels (IQR):\n')
    fprintf('      %s : %\n',string(data_norm.label(badChans)))

    fprintf("\n(3) Per-trial Artifact Rejection ...\n")

    artifacts_byChan =  cell(Nchan,1); ieds_byChan = cell(Nchan,1); iqr_byChan = cell(Nchan,1); range_byChan= cell(Nchan,1);
    for chan = 1:Nchan
        clear dat_trl_mean dat_trl_std dat_trl_raw trl2large trl2small


        %------------------%
        % 1: IQR rejection %
        %------------------%
        fprintf('    For Channel %s \n',data_norm.label{chan})
    
       iqr_artifacts = [];

        currTrial = norm_concat(chan,:);
        data_iqr = iqr(currTrial);
        channelBaseline = median(currTrial, 'omitnan');

        for begin = 1:250:length(currTrial)

            if begin + 250 < length(currTrial)
                current_segment_time = begin:begin + 250;
            else
                current_segment_time = begin:length(currTrial);
            end

            current_segment = currTrial(current_segment_time);
            segmentCenter = median(current_segment, 'omitnan');

            if abs(segmentCenter - channelBaseline) >= w * data_iqr
                iqr_artifacts = [iqr_artifacts current_segment_time];
            end
        end
        
        fprintf('        Channel IQR artifs (n.samples) : %d \n',length(iqr_artifacts))


        %------------------%
        % 2:  Range %
        %------------------%
    
        % (1) 225 Hz FIltered signal
        filtered_range = ft_preproc_highpassfilter(norm_concat(chan,:), data_norm.fsample,225,[],'fir','twopass'); 
        dat_filtered = zscore(squeeze(filtered_range));

        % (2) differential Z-score 
        dat_diff= zscore(diff(squeeze(norm_concat(chan,:)))); % Plain Differential Signal 

        % (3) BHA zscore
        bha_diff = zscore(diff(squeeze(bha_concat(chan,:))));

        % Finding artefacts that exceed >5 
        artif_estim = sort(unique([find(abs(dat_diff) >range_threshold) find(abs(dat_filtered) >range_threshold)   find(abs(bha_diff) >range_threshold)]));
        artif_estim = artif_estim(~ismember(artif_estim, boundary_samples)); % Removing those artifacts that interject on boundary samples
        artifacts_env = {};

        % Making clusters by concatenating 5 consecutive sampes 
        i = 1;
        while i <= length(artif_estim)
            cluster_start = i;
            while i < length(artif_estim) && artif_estim(i+1) - artif_estim(i) <= cluster_tolerance
                i = i + 1;
            end
            artifacts_env{end+1} = artif_estim(cluster_start:i); %
            i = i + 1;
        end

        fprintf('        Channel Range artifs (n.events): %d \n',length(artifacts_env))

 
        temp_rejection = [];
        %%% Get artifact Timestamps by different means: 
        for i = 1:numel(artifacts_env)
        
            thisArtifact = artifacts_env{i};
        
            if isempty(thisArtifact); continue; end
        
            % Działa dla pojedynczego sampla i wektora sampli
            startSample = max(1, min(thisArtifact) - artPadding);
            endSample   = min(length(data_artf), max(thisArtifact) + artPadding);
        
            temp_rejection = [temp_rejection, startSample:endSample];
        end


        %%% Join Together artefact Timestamps! (And filter out Trial Boundaries
        temp_ar2tif_timestamps = unique(sort([ied_timestamps{channel}, temp_rejection, iqr_artifacts]));
        temp_ar2tif_timestamps = temp_ar2tif_timestamps(~ismember(temp_ar2tif_timestamps, boundary_samples)); % Removing those artifacts that interject on boundary samples


        %%% Filter out data artifacts for IEDs, Range, IQR and JOINT
        data_ied_nan = data_artf; data_ied_nan(chan,ied_timestamps{channel}) = NaN;
        data_range_nan = data_artf; data_range_nan(chan,temp_rejection) = NaN;
        data_iqr_nan = data_artf; data_iqr_nan(chan,iqr_artifacts) = NaN;

        data_artf(chan,temp_ar2tif_timestamps) = NaN;
  

        %%% Divide into trials once again (and create Per-trial artifact timestamps 
        trials_artf = mat2cell(data_artf(chan,:)', repmat(length(data_norm.trial{1}), 1, length(data_norm.trial)))';
        trials_ied  = mat2cell(data_ied_nan(chan,:)', repmat(length(data_norm.trial{1}), 1, length(data_norm.trial)))';
        trials_iqr  = mat2cell(data_iqr_nan(chan,:)', repmat(length(data_norm.trial{1}), 1, length(data_norm.trial)))';
        trials_rage = mat2cell(data_range_nan(chan,:)', repmat(length(data_norm.trial{1}), 1, length(data_norm.trial)))';
        
        artifacts_byTrial = cell(1,80); ieds_byTrial = cell(1,80); iqr_byTrial = cell(1,80); range_byTrial = cell(1,80); 
        for i = 1:length(trials_artf)
            artifacts_byTrial{i} = find(isnan(trials_artf{i}));
            ieds_byTrial{i} = find(isnan(trials_ied{i}));
            iqr_byTrial{i}  = find(isnan(trials_iqr{i}));
            range_byTrial{i} = find(isnan(trials_rage{i}));
        end

        artifacts_byChan{chan} = artifacts_byTrial; ieds_byChan{chan} = ieds_byTrial; iqr_byChan{chan} = iqr_byTrial; range_byChan{chan} = range_byTrial;
        fprintf('        Trial n. with artifacts        : %d / %d \n',sum(cellfun(@length, artifacts_byChan{chan}) > 0),length(artifacts_byChan{chan}))
        fprintf('        Trials with long artifs (0.1)  : %d / %d \n',sum(cellfun(@length, artifacts_byChan{chan}) > 700),length(artifacts_byChan{chan}))
        fprintf('        Average artif Length (samples) : %f \n',mean(cellfun(@length, artifacts_byChan{chan})))

    end
    data_artifact = data_norm;
    data_artifact.trial = mat2cell(data_artf, size(data_norm.trial{1},1), repmat(length(data_norm.trial{1}), 1, length(data_norm.trial)));
  
    
    clear artifact_clusters artifs chan clust_count
end