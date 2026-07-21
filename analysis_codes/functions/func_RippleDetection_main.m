
function [hipothetical_ripples_all,varargout] = func_RippleDetection_main(RippleData,current_data_viz,std_chan,mean_chan,artifacts_byChan,channel,params,logFilePath,saveFolder,varargin)
    

    if ~isempty(varargin)
        unfilteredData = varargin{1};
        currentChannelName = varargin{2};
    end

    fid = fopen(logFilePath, 'w');
    hipothetical_ripples_all = cell(length(RippleData.trial),1);
    total_ripple_number = 0;

    detection_log = zeros(length(RippleData.trial), 7);
    for it = 1:length(RippleData.trial) %    Iterate through trials

        currentRippleData = RippleData.trial{it}(channel,:);
        hipothetical_ripples = zeros(1,length(currentRippleData));   % first Hypothetical Ripple Struct
        


       
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 1 Detection:  %%%%
        % % Checking for all moments when the amplitude was higher than 2 SDs
    
        fprintf("1 Detection (2 SD) ...")
        
        for i = 1:length(currentRippleData)
            if abs(mean_chan(channel)-currentRippleData(i)) > 2*std_chan(channel)
                hipothetical_ripples(i) = 1; %% One only when current time point can be a ripple candidate
            end
        end
        ripp_only = find(hipothetical_ripples); % get only times when there is a ripple candidate

        if isempty(ripp_only)
            warning("No Ripples found in 1st Stage!")
            continue;
        end
    
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
        detection_log(it,1) = length(hipothetical_ripples_cluster);
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster)); % Remove empty Cells
        
        %%%%%%%%%%%%%%%%%%%%%%%
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 2 Detection:  %%%%
        % Checking for all moments when the select high amplitude clusters
        % contains at least three 4 std peak!
    
        fprintf("2 Detection (4 SD) ...")
        
        for clust = 1:length(hipothetical_ripples_cluster)
            curr_clust = hipothetical_ripples_cluster{clust};
            if sum(abs(currentRippleData(curr_clust)-mean_chan(channel)) >= params.maximum_peak*std_chan(channel)) >1
                continue
            else
                hipothetical_ripples_cluster{clust} = {}; %If no >4Std, make an empty Cluster
            end
        end
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
        fprintf("done! \n")
       detection_log(it,2) = length(hipothetical_ripples_cluster);

        %%%%%%%%%%%%%%%%%%%%%%%
    
     
        
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 3 Detection:  %%%%
        % Check whether ripple is not above maximum length
        fprintf("3 Detection (max) ....")
         
        % % Checking Whether "ripple" is not Too long
        for clust = 1:length(hipothetical_ripples_cluster)
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
        detection_log(it,3) = length(hipothetical_ripples_cluster);

        %%%%%%%%%%%%%%%%%%%%%%%
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 4 Detection:  %%%%
        % Joining together ripple clusters that are close
        fprintf("4 Detection (join) ...")
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
        detection_log(it,4) = length(hipothetical_ripples_cluster);

        %%%%%%%%%%%%%%%%%%%%%%%
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 5 Detection:  %%%%
        % Check whether Ripple has at least 3 cycles
        fprintf("5 Detection (%s) ..",params.detection_peaks)

        if strcmp(params.detection_peaks, "peaks") % by peaks
         
            for clust = 1:length(hipothetical_ripples_cluster)
    
               curr_clust = current_data_viz.trial{it}(channel,hipothetical_ripples_cluster{clust}); % Currenr Ripple
               [peaks ,locs,w,p] = findpeaks(curr_clust,500);  % Peaks Computation
    
               if length(peaks) >= params.num_peaks
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                end
            end
        elseif strcmp(params.detection_peaks, "thrsh") % by threshold
        
    
            for clust = 1:length(hipothetical_ripples_cluster)
                curr_clust = hipothetical_ripples_cluster{clust};
                if length(curr_clust) >= params.sampleThresh
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                end
            end
        end
        %%
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
        fprintf("done! \n")
        detection_log(it,5) = length(hipothetical_ripples_cluster);

    
        %%%%%%%%%%%%%%%%%%%%%%%

    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 6 Detection: %%%%%
        % % Remove Artifacts
        fprintf("6 Detection (artif) ..")
        currentArtifacts = artifacts_byChan{channel}{it};
        %%
        for clust = 1:length(hipothetical_ripples_cluster)
        
            coverage = setdiff(hipothetical_ripples_cluster{clust},currentArtifacts);
            if length(coverage) == length(hipothetical_ripples_cluster{clust})
                if sum(coverage == hipothetical_ripples_cluster{clust}) ==length(hipothetical_ripples_cluster{clust})
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                    break
                end
            else
                hipothetical_ripples_cluster{clust} = {};    
                break
            end
        end
            %%
        clear currentArtifacts clust coverage numCov artf clust
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
        fprintf("done! \n")
        detection_log(it,6) = length(hipothetical_ripples_cluster);

        %%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%
        %%% 7 Detection: %%%%%
        % % Find Spectral Peaks 

        if params.spectralPeak 
            for clust = 1:length(hipothetical_ripples_cluster)
                current_cluster = (hipothetical_ripples_cluster{clust}(1):hipothetical_ripples_cluster{clust}(end));
                currentData = currentRippleData(hipothetical_ripples_cluster{clust});
                [~, ripplePeakIdx] =  max(abs(currentData));
                ripple_peak = current_cluster(ripplePeakIdx);
        
                if ripple_peak-params.SpectralPeak_segmL < 1 || ripple_peak+params.SpectralPeak_segmL > size(currentRippleData,2)
                    hipothetical_ripples_cluster{clust} = {};
                    continue
                end
                seg = unfilteredData{it}(channel,ripple_peak-params.SpectralPeak_segmL : ripple_peak+params.SpectralPeak_segmL);
        
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
                %[~, locs] = findpeaks(freq_profile(freqs >= 65 & freqs <= 135), freqs(freqs >= 65 & freqs <= 135), 'MinPeakProminence', params.SpectralPeak_Prominence*max(freq_profile(freqs >= 65 & freqs <= 135)));  
                 [~, locs] = findpeaks(freq_profile(freqs >= params.SpectralPeak_freqMin & freqs <= params.SpectralPeak_freqMax),...
                                             freqs(freqs >= params.SpectralPeak_freqMin & freqs <= params.SpectralPeak_freqMax), ...
                                             'MinPeakProminence', params.SpectralPeak_Prominence*max(freq_profile(freqs >= params.SpectralPeak_freqMin & freqs <= params.SpectralPeak_freqMax)));  
                
                % % Manual Algorithm Check first
                % freqogr = freqs(freqs >= 65 & freqs <= 135);
                % figure("Visible","Off"); 
                % try
                %     sgtitle(sprintf("Trial: %s",string(it)))
                %     subplot(2,2,[2 4]);imagesc(abs(pow(freqs >= 65 & freqs <= 135,:))); yticklabels(freqogr(1:10:length(freqogr))); title("Power"); axis xy 
                %     subplot(2,2,1); plot(freq_profile(freqs >= 65 & freqs <= 135)); xticklabels(freqogr(1:10:length(freqogr))); title("Frequency Power Profile")
                %     subplot(2,2,3); plot(unfilteredData{it}(channel,hipothetical_ripples_cluster{clust}(1)-50:hipothetical_ripples_cluster{clust}(end)+50));  title("Ripple Timecourse")
                %     xline(50)
                %     xline(hipothetical_ripples_cluster{clust}(end)-hipothetical_ripples_cluster{clust}(1)+50)
                % 
                %     xticks([1,round(median(hipothetical_ripples_cluster{clust}))-hipothetical_ripples_cluster{clust}(1)+50,hipothetical_ripples_cluster{clust}(end)-hipothetical_ripples_cluster{clust}(1)+100])
                %     xticklabels([-100,0,100])
                % catch
                %     warning("There was a mistake while generating figure")
                % end
                % print(gcf,fullfile(saveFolder,sprintf("%s_ripplePowerProfile.ps",currentChannelName)),'-dpsc','-append','-bestfit');
                % If peak location isoutside of Ripple Range - Do not save a ripple
                if sum(locs >= params.SpectralPeak_freqMin & locs <= params.SpectralPeak_freqMax) 
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                end
                close(gcf);

            end
            close all;
            hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));

        end
        detection_log(it,7) = length(hipothetical_ripples_cluster);

        %%%%%%%%%%%%%%%%%
            
        
    
        %%
        clear ourSpikes numCov coverage spike clust ourSpikesClusters
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));

        hipothetical_ripples_all{it} = hipothetical_ripples_cluster;
        fprintf("done! \n")

        total_ripple_number = total_ripple_number + length(hipothetical_ripples_cluster);
        
    end

    %%% Getting Trial-wise ripple LOG
    T = array2table(detection_log, ...
        'VariableNames', {'Det1_2SD','Det2_4SD','Det3_MaxLen','Det4_Join','Det5_Cycles','Det6_Artif_Rej','Det7_Spctr_Peak'});
    T.Trial = (1:length(RippleData.trial))';
    T = movevars(T, 'Trial', 'Before', 1);
    disp(T)
    fprintf(fid, '%-6s %-10s %-10s %-13s %-11s %-13s %-16d %-19d\n\n', ...
        'Trial','Det1_2SD','Det2_4SD','Det3_MaxLen','Det4_Join','Det5_Cycles','Det6_Artif_Rej','Det7_Spctr_Peak');
    for i = 1:height(T)    % Rows
        fprintf(fid, '%-6d %-10d %-10d %-13d %-11d %-13d %-16d %-19d\n', ...
            T.Trial(i), T.Det1_2SD(i), T.Det2_4SD(i), T.Det3_MaxLen(i), ...
            T.Det4_Join(i), T.Det5_Cycles(i), T.Det6_Artif_Rej(i), T.Det7_Spctr_Peak(i));
    end
    fprintf(fid, '%s\n', repmat('-',1,82));    % sum

    fprintf(fid, '%-6s %-10d %-10d %-13d %-11d %-13d %-16d %-19d\n\n', ...
        'SUM', sum(T.Det1_2SD), sum(T.Det2_4SD), sum(T.Det3_MaxLen), ...
        sum(T.Det4_Join), sum(T.Det5_Cycles), sum(T.Det6_Artif_Rej), sum(T.Det7_Spctr_Peak));

    fprintf(fid, '%s\n', repmat('-',1,82));    fprintf(fid, '\n Total Ripple Number: %d\n', total_ripple_number);
    fprintf(fid, '\n Detection parameters:\n');
    fprintf(fid, '     maximum_peak:    %s\n', string(params.maximum_peak));
    fprintf(fid, '     maxRippleLength: %s\n', string(params.maxRippleLength));
    fprintf(fid, '     detection_peaks: %s\n', string(params.detection_peaks));
    fprintf(fid, '     num_peaks:       %s\n', string(params.num_peaks));
    if isfield(params,"samplThresh")
        fprintf(fid, '     sampleThresh:    %s\n', string(params.sampleThresh));
    end
    fclose(fid);


end


