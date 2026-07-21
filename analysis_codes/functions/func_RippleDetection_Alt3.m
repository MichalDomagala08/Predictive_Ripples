


function [hipothetical_ripples_cluster,varargout] = func_RippleDetection_Alt3(currentRippleData,data_viz,fileID,ObservationsCount,channel,params,unfilteredData)
    hipothetical_ripples = zeros(1,length(currentRippleData));   % first Hypothetical Ripple Struct

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 1-2 Detection: %%%
    % RMS Cacluation (Fell Only)  %%%%

    fprintf("1 Detection (2.5 ~ 9 RMS) ...")

    fs = 500;
    win = round(0.02*fs); % 20 ms
            
    rms_signal = sqrt(movmean(currentRippleData.^2, win,'omitnan'));
    
    mu = mean(rms_signal,"omitnan");
    sigma = std(rms_signal,"omitnan");
    
    lower_thr = mu + 2.5*sigma;
    upper_thr = mu + 9*sigma;
    hipothetical_ripples = find(rms_signal > lower_thr & rms_signal < upper_thr);    
    ripp_only = hipothetical_ripples;

    % Clustering of all Hippotethical Ripples For Clarity Purposes
    hipothetical_ripples_cluster = {}; % cell array with Ripple Clusters for current channel
    cluster_count = 1;
    cluster_contain = (ripp_only(1));
    for i = 1:length(ripp_only)-1
        if abs(ripp_only(i) - ripp_only(i+1)) <=2
            cluster_contain = [cluster_contain ripp_only(i+1)];
        else
            cluster_count = cluster_count+1;
            cluster_contain =  (ripp_only(i+1));
        end
        hipothetical_ripples_cluster{cluster_count} = cluster_contain;
    end
    
    clear cluster_contain cluster_count ripp_only hipothetical_ripples
    fprintf("done! \n")
    fprintf(fileID,'1 Detection (2 Std): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster)); % Remove empty Cells
    %%%%%%%%%%%%%%%%%%%%%%%
    
    

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 3 Detection:  %%%%
    fprintf("3 Detection (max) ....")
    
    
    % % Checking Whether "ripple" is not Too long
    parfor clust = 1:length(hipothetical_ripples_cluster)
        curr_clust = hipothetical_ripples_cluster{clust};
        %%% Here WE just check whethr it is not too long, and whether
        %%% it can Sustain 3 Cycles --> Later, as Ripple could be
        %%% Shorter than 38 ms (3 cycles and 80Hz) We just checks
        %%% Number of Peaks!!!!
    
        if length(curr_clust) <params.maxRippleLength && length(curr_clust) > 3 %floor(1/80*3*500)
            continue
        else
            hipothetical_ripples_cluster{clust} = {};
        end
    end
    
    % Remove Empty Clusters
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    
    fprintf(fileID,'3 Detection (Max 200): ripple candidates: %d \n',length(hipothetical_ripples_cluster));
    %%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 4 Detection:  %%%%
    % Detecting at least 3 cycles
    fprintf("5 Detection (thrpk) ...")
    parfor clust = 1:length(hipothetical_ripples_cluster)

        curr_clust = data_viz.trial{1}(channel,hipothetical_ripples_cluster{clust}); % Currenr Ripple
        [peaks ,locs,w,p] = findpeaks(curr_clust,500);  % Peaks Computation

       if length(peaks) >= params.num_peaks  && length(curr_clust) >= params.sampleThresh
            continue
        else
            hipothetical_ripples_cluster{clust} = {};
        end
    end
    %%
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")

    fprintf(fileID,'5 Detection (thr_pk): ripple candidates:  %d \n',length(hipothetical_ripples_cluster));
    
    clear VizTimeCourse margins curr_clust addition1 addition2 ripple_time ripple_timecourse temp_cmplx pov freq timestamp curr_clust_data peaks locs w p std_chan mean_chan
    %%%%%%%%%%%%%%%%%%%%%%%



   
    %%%%%%%%%%%%%%%%%%%%%%
    %%% 5 Detection: %%%%%
    % % Find Spectral Peaks 
    for clust = 1:length(hipothetical_ripples_cluster)
        current_cluster = (hipothetical_ripples_cluster{clust}(1):hipothetical_ripples_cluster{clust}(end));
        currentData = currentRippleData(hipothetical_ripples_cluster{clust});
        [~, ripplePeakIdx] =  max(abs(currentData));
        ripple_peak = current_cluster(ripplePeakIdx);

        if ripple_peak-params.SpectralPeak_segmL < 1 || ripple_peak+params.SpectralPeak_segmL > size(currentRippleData,2)
            hipothetical_ripples_cluster{clust} = {};
            continue
        end
        seg = unfilteredData(ripple_peak-params.SpectralPeak_segmL : ripple_peak+params.SpectralPeak_segmL);

        %%% NAN remover: not to confound the signal:
        if any(isnan(seg))
            seg(isnan(seg)) = mean(seg, 'omitnan');
        end

        [pow,freqs,times ] = ft_specest_wavelet(seg, (1:length(seg))/500, 'width', 7, 'foi', 65:2:135,  'verbose', 0,  'pad', 1);              % padding dla FFT
        pow = abs(squeeze(pow)); %%% Getting Amplitude Values
        if params.SpectralPeak_PowerScaling
            pow = pow.*freqs(:); %%% Getting Power Scaling 
        end

        %%% Collapsing Temporal Indexes
        center_idx =ceil(length(times)/2);      % peak jest w czasie 0 w segmencie
        t_idx = center_idx- round(params.SpectralPeak_tfrEnv * 500) : center_idx+ round(params.SpectralPeak_tfrEnv * 500);
        freq_profile = mean(pow( :, t_idx), 2,'omitnan');  % freq x 1

        %%% Finding Peaks 
        [~, locs] = findpeaks(freq_profile(freqs >= 65 & freqs <= 135), freqs(freqs >= 65 & freqs <= 135), 'MinPeakProminence', params.SpectralPeak_Prominence*max(freq_profile(freqs >= 65 & freqs <= 135)));  
        
        % Manual Algorithm Check first
        % figure(); imagesc(abs(pow(freqs >= 65 & freqs <= 135,:))); axis xy 
        % figure(); plot(freq_profile(freqs >= 65 & freqs <= 135))

        % If peak location isoutside of Ripple Range - Do not save a ripple
        if sum(locs >= params.SpectralPeak_freqMin & locs <= params.SpectralPeak_freqMax) 
            continue
        else
            hipothetical_ripples_cluster{clust} = {};
        end
        close all;

    end
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf(fileID,'8 Spectral Peak Present: ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
    %%%%%%%%%%%%%%%%%
    

    
    
    % viz plot
    if params.vizualization == 2
        %%% Preprocess Spike for Visualizations
        ourSpikesClusters = cell(1,length(ourSpikes));
        for i = 1:length(ourSpikes)
            ourSpikesClusters{i} = [ourSpikes(i,1):ourSpikes(i,2)];
        end
    
        pathSpikes = fullfile(saveFolder,params.save_modifier,"Ripples",current_ds_name,strcat(goodChannels{channel},"_Spikes",".ps"));
        %%% Visualize Spikes Without saving any in memory
        [~,~] =  getSubjRipples(params,ourSpikesClusters,data_viz,data_norm,channel,goodChannels,current_ds_name,pathSpikes);
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pathSpikes,".ps", ".pdf"),pathSpikes))
        delete(pathSpikes)
    elseif params.vizualization ==3
        %%% Fisualize Faulty Channels
        pathSpikesUndetected = fullfile(saveFolder,params.save_modifier,"Ripples",current_ds_name,strcat(goodChannels{channel},"_SpikesUndet",".ps"));
        faultyRipples = {hipothetical_ripples_cluster{[artefacts{ObservationsCount}]}};
        [~,~] = getSubjRipples(params,faultyRipples,data_viz,data_norm,channel,goodChannels,current_ds_name,pathSpikesUndetected);
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pathSpikesUndetected,".ps", ".pdf "),pathSpikesUndetected))
        delete(pathSpikesUndetected)
    end
        
    
    
    %%
    clear ourSpikes numCov coverage spike clust ourSpikesClusters
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    


end