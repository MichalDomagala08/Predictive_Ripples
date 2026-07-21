


function [hipothetical_ripples_cluster,varargout] = func_RippleDetection_Alt2(currentRippleData,data_viz,std_chan,mean_chan,fileID,ObservationsCount,channel,params)
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

 
    

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 3 Detection:  %%%%
    % Check whether ripple is not above maximum length
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
    
    fprintf(fileID,'4 Detection (join): ripple candidates:    %d \n',length(hipothetical_ripples_cluster));
    %%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% 5 Detection:  %%%%
    
    % Check whether Ripple has at least 3 cycles
    
    fprintf("5 Detection (thrsh) ..")

    parfor clust = 1:length(hipothetical_ripples_cluster)
        curr_clust = hipothetical_ripples_cluster{clust};
        if length(curr_clust) >= params.sampleThresh
            continue
        else
            hipothetical_ripples_cluster{clust} = {};
        end
    end
    hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));
    fprintf("done! \n")

    fprintf(fileID,'5 Detection (thresh): ripple candidates:  %d \n',length(hipothetical_ripples_cluster));

    clear VizTimeCourse margins curr_clust addition1 addition2 ripple_time ripple_timecourse temp_cmplx pov freq timestamp curr_clust_data peaks locs w p std_chan mean_chan
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
 
    
    % Viz
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