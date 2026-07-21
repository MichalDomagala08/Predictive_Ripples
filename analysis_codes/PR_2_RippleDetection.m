
ft_defaults()
%%% Ripple Detection Parameters %%%
params.lowpassfreq     = 80;  % cutoffs for filtering ripple band
params.highpassfreq    = 120;
params.cluster_join    = 7;   % timepoint distance on which clusters will be joined %15ms                                            
params.maxRippleLength = 100; % Maximum length (in samples of a ripple - Current 100ms)
params.maximum_peak    = 4;   % SDs of maximum peak from mean       
params.num_peaks       = 3;   % minimum peaks that are detected in ripple                                                             
params.detection_peaks = "peaks"; % thrsh - lowerbound threshold; 
                                  % peaks - minimum peak attained;
params.sampleThresh= 18;



params.vis = 3; %  
                          
params.vizbefT = 250; %How many timestamps before max Ripple for Viz! 
params.vizaftT = 250; 
params.FFT_Length = 250; %
params.fftmaxFreq = 150; % Viz FFT parameters
params.fftminFreq = 60;  % Viz FFT parameters

fieldtripPath = 'C:\Users\barak\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118';  % <- zmień na rzeczywistą ścieżkę
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions')
                                  
% Paths
ripple_detection_scheme_name = "RippleDetection";
dataFolder = "D:\Documents_Dell\Predictive_Ripples_2025\preprocessed";
saveFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",ripple_detection_scheme_name);
mkdir(saveFolder);

files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractBefore(x, '_artif_rej'), '.mat'), files, 'UniformOutput', false);

% Single subject ( for now) 

% load and unpack the data structures
iter = 1;
data = load(fullfile(dataFolder,files{iter}));
currSubjName = subjNames{iter};

data.badChans = 5; % I need to load it not manualy but by: badChans = 

data_eeg = data.short_data;
data_artifact = data.data_artifact;
artifacts_byChan = data.artifacts_byChan(setdiff(1:length(data.artifacts_byChan),data.badChans));
iqr_byChan = data.iqr_byChan(setdiff(1:length(data.iqr_byChan),data.badChans));
range_byChan = data.range_byChan(setdiff(1:length(data.range_byChan),data.badChans));
ieds_byChan = data.ieds_byChan(setdiff(1:length(data.ieds_byChan),data.badChans));

good_trials = data.good_trials(setdiff(1:length(data.good_trials),data.badChans));


%%% get Ripple Envelope: 
[data_ripples,data_viz] = func_rippleband_filtering(data_eeg,params);

%%% Get Analytical Signals used for establishing all
data_ripples_zscore = data_ripples;
data_artif = data_ripples;

all_rippledata_viz = zscore(cell2mat(data_ripples.trial),0,2); % Zscoring the data 
%re-join trials: 
for tr = 1:length(data_artif.trial)
    data_ripples_zscore.trial{tr} = all_rippledata_viz(:,(tr-1)*length(data_artif.trial{1}) +1:tr*length(data_artif.trial{1}));
    for chan = 1:length(data_artif.label)
        data_artif.trial{tr}(chan,isnan(data_artifact.trial{tr}(chan,:))) = nan;
    end
end


joint_artif_data = cell2mat(data_artif.trial);
std_chan = std(joint_artif_data',1,"omitnan");
mean_chan = mean(joint_artif_data',1,"omitnan");

path_ripple_mean_save = fullfile(saveFolder, "RippleChars_XXX_.ps");
path_ripple_mean_indiv = fullfile(saveFolder, currSubjName, "RippleChars_XXX_.png");
path_artif_mean_save = fullfile(saveFolder, "ArtifChars_XXX_.ps");
path_artif_mean_indiv = fullfile(saveFolder, currSubjName, "ArtifChars_XXX_.png");



addpath(fieldtripPath);
ft_defaults;


% Konfiguracja każdego workera
%%
for channel = 1:length(data_eeg.label)
    fid = fopen(fullfile(saveFolder,currSubjName,sprintf('%s_ripple_detection_log.txt',data_eeg.label{channel})), 'w');
    [hipothetical_ripples_cluster] = func_RippleDetection_main(data_ripples,data_viz,std_chan,mean_chan,artifacts_byChan,channel,params,fid,fullfile(saveFolder,currSubjName));
            
    fprintf("     Channel: %s",data_eeg.label{channel})


   %%% VISUALIZE DATA


    psDir   = fullfile(saveFolder, currSubjName,sprintf('%s_trial_ripple.ps',data_artif.label{channel}));
    

    powerSWR_all      = [];timecourseSWR_all = []; FFT_SWR_all       = [];
    powerARF_all      = [];timecourseARF_all = []; FFT_ARF_all       = [];

    path_ripple_save = fullfile(saveFolder,currSubjName, "single_ripples_"+data_artif.label{channel}+".ps");
    path_artif_save = fullfile(saveFolder,currSubjName, "single_artifs_"+data_artif.label{channel}+".ps");

    artifs_list = {};
    for it = 1:length(data_ripples.trial)

        fprintf("     Trial N: %d",it)

        idx = artifacts_byChan{channel}{it}; mask = false(1, max(idx)); mask(idx) = true;
        stats = regionprops(logical(mask), 'PixelIdxList', 'Area');
        artifactsClustered = {stats([stats.Area] >= 5).PixelIdxList};
        artifs_list = {artifs_list artifactsClustered};
        %%% Get Artifact Characteristic Data
        [powerARF,timecourseARF,ripplerangeARF,FFT_ARF] = func_RippleChar(params,artifactsClustered,data_viz.trial{it},data_eeg.trial{it},channel);
        powerARF_all      = cat(1, powerARF_all,      powerARF);
        timecourseARF_all = cat(1, timecourseARF_all, timecourseARF);
        FFT_ARF_all       = cat(1, FFT_ARF_all,       FFT_ARF);

        %%% Get Ripple Characteistic Data 
        [powerSWR,timecourseSWR,ripplerangeSWR,FFT_SWR] = func_RippleChar(params,hipothetical_ripples_cluster{it},data_viz.trial{it},data_eeg.trial{it},channel);
        powerSWR_all      = cat(1, powerSWR_all,      powerSWR);
        timecourseSWR_all = cat(1, timecourseSWR_all, timecourseSWR);
        FFT_SWR_all       = cat(1, FFT_SWR_all,       FFT_SWR);
        
        %%% (1) Detection Timecourses with marked ripples 
        if params.vis
            plot_ripple_trials(it, channel, hipothetical_ripples_cluster, data_viz, data_eeg, data_ripples_zscore, psDir);
        end
        

        %%% (2) Optional - All Detected Ripple Characteristics (with 500 ms boundaries) 
        if params.vis > 1
            for rip = 1:length(hipothetical_ripples_cluster{it})
               plot_single_ripple(powerSWR(rip,:,:),timecourseSWR(rip,:),ripplerangeSWR(rip,:), currSubjName,data_eeg.label{channel},params,path_ripple_save)
            end
        end

        %%% (3) Optional - All Artifact Characteristics (with 500 ms boundaries) 
        if  params.vis == 3
            for rip = 1:length(artifactsClustered)
               plot_single_ripple(powerARF(rip,:,:),timecourseARF(rip,:),ripplerangeARF(rip,:), currSubjName,data_eeg.label{channel},params,path_artif_save)
            end
        end

        %%
    end

    %%% (4) Mean ripple characteristics     result = [hipothetical_ripples_cluster{:}];

    plot_rippleSignal(timecourseSWR_all,mean(powerSWR_all,1,"omitnan"),FFT_SWR_all,diff(cellfun(@median,[hipothetical_ripples_cluster{:}])),...
                    size(timecourseSWR_all,1),path_ripple_mean_save,params,currSubjName+" // "+data_artif.label{channel},path_ripple_mean_indiv)
    

    %%% (5) Mean Artif characteristics
    plot_rippleSignal(mean(timecourseARF_all,1,"omitnan"),mean(powerARF_all,1,"omitnan"),mean(FFT_ARF_all,1,"omitnan"),diff(cellfun(@median,artifs_list)),...
                    size(timecourseARF_all,1),path_artif_mean_save,params,currSubjName+" // "+data_artif.label{channel},path_artif_mean_indiv)
    

    % Transform PS into PDFs
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir)
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_ripple_save,".ps", ".pdf "),path_ripple_save));  delete(path_ripple_save)
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_artif_save,".ps", ".pdf "),path_artif_save));  delete(path_artif_save)


end
%%










function [hipothetical_ripples_all,varargout] = func_RippleDetection_main(RippleData,current_data_viz,std_chan,mean_chan,artifacts_byChan,channel,params,fid,saveFolder)
    
    hipothetical_ripples_all = cell(length(RippleData.trial),1);
    total_ripple_number = 0;

    detection_log = zeros(length(RippleData.trial), 6);
    for it = 1:length(RippleData.trial) %    Iterate through trials

        currentRippleData = RippleData.trial{it}(channel,:);
        hipothetical_ripples = zeros(1,length(currentRippleData));   % first Hypothetical Ripple Struct
        


       
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 1 Detection:  %%%%
        % % Checking for all moments when the amplitude was higher than 2 SDs
    
        fprintf("1 Detection (2 SD) ...")
        
        parfor i = 1:length(currentRippleData)
            if abs(mean_chan(channel)-currentRippleData(i)) > 2*std_chan(channel)
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
        detection_log(it,1) = length(hipothetical_ripples_cluster);
        %fprintf(fileID,'1 Detection (2 Std): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
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
       % fprintf(fileID,'2 Detection (4 Std): ripple candidates:   %d \n',length(hipothetical_ripples_cluster));
       detection_log(it,2) = length(hipothetical_ripples_cluster);

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
        detection_log(it,3) = length(hipothetical_ripples_cluster);

      %  fprintf(fileID,'3 Detection (Max 200): ripple candidates: %d \n',length(hipothetical_ripples_cluster));
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

      %  fprintf(fileID,'4 Detection (join): ripple candidates:    %d \n',length(hipothetical_ripples_cluster));
        %%%%%%%%%%%%%%%%%%%%%%%
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%% 5 Detection:  %%%%
        % Check whether Ripple has at least 3 cycles
        fprintf("5 Detection (%s) ..",params.detection_peaks)

        if strcmp(params.detection_peaks, "peaks") % by peaks
         
            parfor clust = 1:length(hipothetical_ripples_cluster)
    
               curr_clust = current_data_viz.trial{it}(channel,hipothetical_ripples_cluster{clust}); % Currenr Ripple
               [peaks ,locs,w,p] = findpeaks(curr_clust,500);  % Peaks Computation
    
               if length(peaks) >= params.num_peaks
                    continue
                else
                    hipothetical_ripples_cluster{clust} = {};
                end
            end
        elseif strcmp(params.detection_peaks, "thrsh") % by threshold
        
    
            parfor clust = 1:length(hipothetical_ripples_cluster)
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

      %  fprintf(fileID,'5 Detection: ripple candidates:  %d \n',length(hipothetical_ripples_cluster));
    
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
        
        
   

        % labeled = bwlabel(artifacts_byChan{channel}{it});  % numeruje ciągłe regiony
        % stats = regionprops(labeled, 'PixelIdxList', 'Area');
        % clusters = {stats([stats.Area] >= 5).PixelIdxList};
        

        
    
        %%
        clear ourSpikes numCov coverage spike clust ourSpikesClusters
        hipothetical_ripples_cluster =  hipothetical_ripples_cluster(~cellfun('isempty',hipothetical_ripples_cluster));

        hipothetical_ripples_all{it} = hipothetical_ripples_cluster;
        fprintf("done! \n")

        total_ripple_number = total_ripple_number + length(hipothetical_ripples_cluster);
        
    end

    %%% Getting Trial-wise ripple LOG
    T = array2table(detection_log, ...
        'VariableNames', {'Det1_2SD','Det2_4SD','Det3_MaxLen','Det4_Join','Det5_Cycles','Det6_Artif_Rej'});
    T.Trial = (1:length(RippleData.trial))';
    T = movevars(T, 'Trial', 'Before', 1);
    disp(T)
    fprintf(fid, '%-6s %-10s %-10s %-13s %-11s %-13s %-16s\n', ...
        'Trial','Det1_2SD','Det2_4SD','Det3_MaxLen','Det4_Join','Det5_Cycles','Det6_Artif_Rej');
    for i = 1:height(T)    % Rows
        fprintf(fid, '%-6d %-10d %-10d %-13d %-11d %-13d %-16d\n', ...
            T.Trial(i), T.Det1_2SD(i), T.Det2_4SD(i), T.Det3_MaxLen(i), ...
            T.Det4_Join(i), T.Det5_Cycles(i), T.Det6_Artif_Rej(i));
    end
    fprintf(fid, '%s\n', repmat('-',1,82));    % sum

    fprintf(fid, '%-6s %-10d %-10d %-13d %-11d %-13d %-16d\n', ...
        'SUM', sum(T.Det1_2SD), sum(T.Det2_4SD), sum(T.Det3_MaxLen), ...
        sum(T.Det4_Join), sum(T.Det5_Cycles), sum(T.Det6_Artif_Rej));

    fprintf(fid, '%s\n', repmat('-',1,82));    fprintf(fid, '\n Total Ripple Number: %s\n', total_ripple_number);
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


