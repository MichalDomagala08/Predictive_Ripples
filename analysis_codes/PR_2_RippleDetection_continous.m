


ft_defaults()
%%% Ripple Detection Parameters %%%
params.lowpassfreq     = 80;  % cutoffs for filtering ripple band
params.highpassfreq    = 120;
params.cluster_join    = 7;   % timepoint distance on which clusters will be joined %15ms                                            
params.maxRippleLength = 100; % Maximum length (in samples of a ripple - Current 100ms)
params.maximum_peak    = 3;   % SDs of maximum peak from mean       
params.num_peaks       = 3;   % minimum peaks that are detected in ripple                                                             
params.detection_peaks = "thrsh"; % thrsh - lowerbound threshold; 
                                  % peaks - minimum peak attained;
params.sampleThresh= 13;

params.spectralPeak = true;
params.SpectralPeak_Prominence = 0.1;    % How much % of a signal are we considering a Peak 
params.SpectralPeak_freqMax = 125;       % Upper Cutoff in mean spectral peak: Properly done should be to Ripple Range
params.SpectralPeak_freqMin = 75;        % Lower Cutoff in mean spectral peak: Properly done should be  to Ripple Range 
params.SpectralPeak_PowerScaling = true; % Whether we are Scaling with Frequencies 
params.SpectralPeak_segmL  = 50;         % time around SWR in miliseconds to get a 
params.SpectralPeak_tfrEnv = 0.05;       % time around SWR in seconds to Average for


params.vis = 3; %  
                          
params.vizbefT = 250; %How many timestamps before max Ripple for Viz! 
params.vizaftT = 250; 
params.FFT_Length = 250; %
params.fftmaxFreq = 150; % Viz FFT parameters
params.fftminFreq = 60;  % Viz FFT parameters

fieldtripPath = 'C:\Users\barak\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118';  % <- zmień na rzeczywistą ścieżkę
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions')
   

% ----- ANALYSIS NAME -----
params.ripple_detection_scheme_name = "RippleDetection_frank_et_al_RipplePeak";
params.artif_rejection_scheme_name = "ArtifactRejection_alternative";

%    -----------------


% Paths

dataFolder = "D:\Documents_Dell\Predictive_Ripples_2025\preprocessed";
saveFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.ripple_detection_scheme_name);
saveDataFolder = fullfile(dataFolder,params.ripple_detection_scheme_name);

mkdir(saveFolder); mkdir(saveDataFolder);


files = {dir(fullfile(dataFolder,params.artif_rejection_scheme_name, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractBefore(x, '_artif_rej'), '.mat'), files, 'UniformOutput', false);

% Single subject ( for now) 

countInstances = 1;

all_hipothetical_ripples_cluster = cell(1,length(subjNames)); 
all_timecourseSWR_trial= cell(1,length(subjNames)); 
all_ripplePeaks = cell(1,length(subjNames)); 

for s = 1:length(subjNames)

    
    % load and unpack the data structures
    data = load(fullfile(dataFolder,params.artif_rejection_scheme_name,files{s}));
    currSubjName = subjNames{s};
        
    data_eeg = data.short_data;
    data_artifact = data.data_artifact;
    artifacts_byChan = data.artifacts_byChan;
    iqr_byChan = data.iqr_byChan;
    range_byChan = data.range_byChan; %(setdiff(1:length(data.range_byChan),data.badChans))
    ieds_byChan = data.ieds_byChan;
    
    good_trials = data.good_trials;
    
    all_hipothetical_ripples_cluster{s} = cell(1,length(data_eeg.label));
    all_timecourseSWR_trial{s} = cell(1,length(data_eeg.label));
    all_ripplePeaks{s} = cell(1,length(data_eeg.label));

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
    path_artif_mean_save = fullfile(saveFolder, "ArtifChars_XXX_.ps");
    
    
    
    
    addpath(fieldtripPath);
    ft_defaults;
    
    
    % Konfiguracja każdego workera
    %%
    mkdir(fullfile(saveFolder, currSubjName));
    
    for channel = 1:length(data_eeg.label)
        logFilePath = fullfile(saveFolder,currSubjName,sprintf('%s_ripple_detection_log.txt',data_eeg.label{channel}));
        [hipothetical_ripples_cluster] = func_RippleDetection_main(data_ripples,data_viz,std_chan,mean_chan,artifacts_byChan,channel,params,logFilePath,fullfile(saveFolder,currSubjName),data_eeg.trial,data_eeg.label{channel});
                
        fprintf("     Channel: %s",data_eeg.label{channel})
    
        path_artif_mean_indiv = fullfile(saveFolder, currSubjName, sprintf("%s_ArtifChars_XXX_.png",data_eeg.label{channel}));
        path_ripple_mean_indiv = fullfile(saveFolder, currSubjName, sprintf("%s_RippleChars_XXX_.png",data_eeg.label{channel}));
       %%% VISUALIZE DATA
    
    
        psDir   = fullfile(saveFolder, currSubjName,sprintf('%s_trial_ripple.ps',data_artif.label{channel}));
        
    
        powerSWR_all      = [];timecourseSWR_all = []; FFT_SWR_all       = [];
        powerARF_all      = [];timecourseARF_all = []; FFT_ARF_all       = [];
    
        path_ripple_save = fullfile(saveFolder,currSubjName, "single_ripples_"+data_artif.label{channel}+".ps");
        path_artif_save = fullfile(saveFolder,currSubjName, "single_artifs_"+data_artif.label{channel}+".ps");
    
        artifs_list = [];
        timecourseSWR_trial = cell(1,length(data_ripples.trial));
        ripplePeaks = cell(1,length(data_ripples.trial));
        for it = 1:length(data_ripples.trial)
    
            fprintf("     Trial N: %d",it)
    
            idx = artifacts_byChan{channel}{it}; mask = false(1, max(idx)); mask(idx) = true;
            stats = regionprops(logical(mask), 'PixelIdxList', 'Area');
            artifactsClustered = {stats([stats.Area] >= 5).PixelIdxList};
            artifs_list = [artifs_list artifactsClustered];
            %%% Get Artifact Characteristic Data
    
            if ~isempty(idx)
                [powerARF,timecourseARF,ripplerangeARF,FFT_ARF] = func_RippleChar(params,artifactsClustered,data_viz.trial{it},data_eeg.trial{it},channel);
                powerARF_all      = cat(1, powerARF_all,      powerARF);
                timecourseARF_all = cat(1, timecourseARF_all, timecourseARF);
                FFT_ARF_all       = cat(1, FFT_ARF_all,       FFT_ARF);
            end
            %%% Get Ripple Characteistic Data 
            [powerSWR,timecourseSWR,ripplerangeSWR,FFT_SWR] = func_RippleChar(params,hipothetical_ripples_cluster{it},data_viz.trial{it},data_eeg.trial{it},channel);
            powerSWR_all      = cat(1, powerSWR_all,      powerSWR);
            timecourseSWR_all = cat(1, timecourseSWR_all, timecourseSWR);
            FFT_SWR_all       = cat(1, FFT_SWR_all,       FFT_SWR);
            timecourseSWR_trial{it} = timecourseSWR;

            %%% Compute Ripple peak: 
            
            I = []; 
            %%
            for clust  = 1:length(hipothetical_ripples_cluster{it})
                [~,It] = max(data_viz.trial{it}(channel, hipothetical_ripples_cluster{it}{clust}));
                I(clust) = It + hipothetical_ripples_cluster{it}{clust}(1);
            end
            %%
            ripplePeaks{it} = I(clust);

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
            if  params.vis == 3 & ~isempty(idx)
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
        plot_rippleSignal(timecourseARF_all,mean(powerARF_all,1,"omitnan"),FFT_SWR_all,diff(cellfun(@median,artifs_list)),...
                        size(timecourseARF_all,1),path_artif_mean_save,params,currSubjName+" // "+data_artif.label{channel},path_artif_mean_indiv)
        
    
        % Transform PS into PDFs
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir)
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_ripple_save,".ps", ".pdf "),path_ripple_save));  delete(path_ripple_save)
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_artif_save,".ps", ".pdf "),path_artif_save));  delete(path_artif_save)

        rpower_path = fullfile(saveFolder,currSubjName,sprntf("%s_ripplePowerProfile.ps",data_eeg.label{channel}))
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(rpower_path,".ps", ".pdf "),rpower_path));  delete(rpower_path)

        all_hipothetical_ripples_cluster{s}{channel} = hipothetical_ripples_cluster;
        all_timecourseSWR_trial{s}{channel} = timecourseSWR_trial;
        all_ripplePeaks{s}{channel}  =  ripplePeaks{it};

    end
end

psPath = fullfile(saveFolder,"ArtifChars_FFT.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Power.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Ripple.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_FFT.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_FFT.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Power.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Ripple.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)


save(fullfile(saveDataFolder,"rippleData.mat"),"all_hipothetical_ripples_cluster","all_timecourseSWR_trial","all_ripplePeaks")
%%









