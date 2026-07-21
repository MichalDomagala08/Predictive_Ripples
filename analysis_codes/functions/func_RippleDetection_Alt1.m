


function [hipothetical_ripples_cluster,varargout] = func_RippleDetection_Alt1(currentRippleData,data_viz,std_chan,mean_chan,artifacts_byChan,currentSpikes,fileID,ObservationsCount,channel,params,unfilteredData)
    hipothetical_ripples = zeros(1,length(currentRippleData));   % first Hypothetical Ripple Struct

   
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 1 Detection:  %%%%
    % % Checking for all moments when the amplitude was higher than 2 SDs

    fprintf("1 Detection (2 SD) ...")
    
    parfor i = 1:length(currentRippleData(:))
        if abs(mean_chan-currentRippleData(i)) > 2*std_chan
            hipothetical_ripples(i) = 1; %% One only when current time point can be a ripple candidate
        end
    end
    ripp_only = find(hipothetical_ripples); % get only times when there is a ripple candidate

  

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
    %%%% 2 Detection:  %%%%
    % Checking for all moments when the select high amplitude clusters
    % contains at least three 4 std peak!


    fprintf("2 Detection (4 SD) ...")
    
    for clust = 1:length(hipothetical_ripples_cluster)
        curr_clust = hipothetical_ripples_cluster{clust};
        if sum(abs(currentRippleData(curr_clust)-mean_chan) >= params.maximum_peak*std_chan) >1
            continue
        else
            hipothetical_ripples_cluster{clust} = {}; %If no >4Std, make an empty Cluster
        end
    end
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    fprintf(fileID,'2 Detection (4 Std): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
    


    %%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%
    %%% 2.5 Detection: %%%
    % Get Rid of Too Big Ripples

    fprintf("2.5 Detection (9std) ..")
    parfor clust = 1:length(hipothetical_ripples_cluster)
        curr_clust = hipothetical_ripples_cluster{clust};
        if sum(abs(currentRippleData(curr_clust)-mean_chan) > 9*std_chan)
            hipothetical_ripples_cluster{clust} = {}; % If there  is a Large peak - get it out! 
        else
            continue
        end
    end

    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    fprintf(fileID,'2.5 Detection (9 Std): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));



    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 3 Detection:  %%%%
    fprintf("3 Detection (max) ....")
    
    % % Checking Whether "ripple" is not Too long
    parfor clust = 1:length(hipothetical_ripples_cluster)
        curr_clust = hipothetical_ripples_cluster{clust};
     
        if length(curr_clust) <params.maxRippleLength && length(curr_clust) > 3
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
    fprintf("4 Detection (join) ...")
    % join clusters that are 10 ms from each other
    for clust = 1:length(hipothetical_ripples_cluster)-1
        curr_clust = hipothetical_ripples_cluster{clust};
        if abs(curr_clust(end) - hipothetical_ripples_cluster{clust+1}(1)) <= params.cluster_join
            hipothetical_ripples_cluster{clust+1} = (curr_clust(1): hipothetical_ripples_cluster{clust+1}(end));
            hipothetical_ripples_cluster{clust} = {};
        else
            continue
        end
    
    end
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    
    fprintf(fileID,'4 Detection (join): ripple candidates:    %d \n',length(hipothetical_ripples_cluster));
    %%%%%%%%%%%%%%%%%%%%%%%
    


    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 5 Detection:  %%%%
    % Check whether Ripple has at least 3 cycles
    fprintf("5 Detection (thrpk) ...")
    parfor clust = 1:length(hipothetical_ripples_cluster)

        curr_clust = data_viz.trial{1}(channel,hipothetical_ripples_cluster{clust}); % Currenr Ripple
        [peaks ,locs,w,p] = findpeaks(curr_clust,500);  % Peaks Computation
        [throughs ,locs,w,p] = findpeaks(-curr_clust,500); % Through Computation
       
       if length(peaks) >= params.num_peaks && length(throughs) >= params.num_peaks && length(curr_clust) >= params.sampleThresh
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


    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 6 Detection: %%%%%
    % % Remove Artifacts

    fprintf("6 Detection (artif) ..")
    
    numCov = 0;
    currentArtifacts = artifacts_byChan{channel};
    parfor clust = 1:length(hipothetical_ripples_cluster)
        for artf =  1:length(currentArtifacts)
    
            coverage = setdiff(hipothetical_ripples_cluster{clust},currentArtifacts{artf});
            if length(coverage) == length(hipothetical_ripples_cluster{clust})
                if sum(coverage == hipothetical_ripples_cluster{clust}) ==length(hipothetical_ripples_cluster{clust})
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                    numCov = numCov+1;
                    break
                end
            else
                hipothetical_ripples_cluster{clust} = {};
                numCov = numCov+1;
    
                break
            end
        end
    end
        
    clear currentArtifacts clust coverage numCov artf clust
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")
    
    fprintf(fileID,'6 Detection (Artif): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% 7 Detection: %%%%%
    % % Remove Ripples that coincide with Inter-ictal Spike events
    fprintf("7 Detection (Spike) ..")
    
    numCov = 0;
    ourSpikes = currentSpikes{channel};
    parfor clust = 1:length(hipothetical_ripples_cluster)
        current_cluster = (hipothetical_ripples_cluster{clust}(1):hipothetical_ripples_cluster{clust}(end));
        for spike =  1:size(ourSpikes,1)
    
            coverage = setdiff(current_cluster,(ourSpikes(spike,1):ourSpikes(spike,2)));
            %disp(coverage == hipothetical_ripples_cluster{clust})
            if length(coverage) == length(current_cluster)
                if sum(coverage == current_cluster) ==length(current_cluster)
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                    numCov = numCov+1;
                    break
                end
            else
                hipothetical_ripples_cluster{clust} = {};
                numCov = numCov+1;
    
                break
            end
        end
    end
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf(fileID,'7 Detection (Spikes): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));    
    %%%%%%%%%%%%%%%%%%%%%




    %%%%%%%%%%%%%%%%%%%%%%
    %%% 8 Detection: %%%%%
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

    

    
    %%% viz code
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