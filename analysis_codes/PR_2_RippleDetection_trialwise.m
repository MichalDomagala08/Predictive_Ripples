
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -- RIPPLE DETECTION SCRIPT -- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script detects Ripple using algorithm adapted from Vaz et al (std
% threshold and peak detection) and Ngo et al (spectral peak detection)
%
% It saves the  underlying ripple timestamps per trial, ripple peaks,
% timcourses with names.
%
% Additionally, when specified, this script can visualise ripples,
% their chars, per channel, trial or individual


ft_defaults()


%%% ----- PARAMETERS ---- %%%

%%% Analysis and Artif Rej Names
params.ripple_detection_scheme_name = "RippleDetection_vaz_et_al_RipplePeak";
params.artifact_reject_scheme_name = "ArtifactRejection_alternative2";
params.manual_rejection = true;

%%% Ripple Detection Parameters %%%
params.lowpassfreq     = 80;  % cutoffs for filtering ripple band
params.highpassfreq    = 120;
params.cluster_join    = 7;   % timepoint distance on which clusters will be joined %15ms                                            
params.maxRippleLength = 100; % Maximum length (in samples of a ripple - Current 100ms)
params.maximum_peak    = 4;   % SDs of maximum peak from mean       
params.num_peaks       = 3;   % minimum peaks that are detected in ripple                                                             
params.detection_peaks = "peaks"; % thrsh - lowerbound threshold; 
                                  % peaks - minimum peak attained;
params.sampleThresh= 13;

%%% Peak detection parameters %%%
params.spectralPeak = true;
params.SpectralPeak_Prominence = 0.1;    % How much % of a signal are we considering a Peak 
params.SpectralPeak_freqMax = 125;       % Upper Cutoff in mean spectral peak: Properly done should be to Ripple Range
params.SpectralPeak_freqMin = 75;        % Lower Cutoff in mean spectral peak: Properly done should be  to Ripple Range 
params.SpectralPeak_PowerScaling = true; % Whether we are Scaling with Frequencies 
params.SpectralPeak_segmL  = 50;         % time around SWR in miliseconds to get a 
params.SpectralPeak_tfrEnv = 0.05;       % time around SWR in seconds to Average for


%%% Visualisation Parameters %%%
params.vis = 3; %  
params.vizbefT = 250; %How many timestamps before max Ripple for Viz! 
params.vizaftT = 250; 
params.FFT_Length = 250; %
params.fftmaxFreq = 150; % Viz FFT parameters
params.fftminFreq = 60;  % Viz FFT parameters


%%% bad channels and trial
% artfi_subj_chan = cell(1,19);
% artfi_subj_chan{3} = {'RDa2-RDa3','RDa3-RDa4','RDa4-RDa5'};  % Qutie allrght
% 
% artfi_subj_chan_tr = cell(1,19);
% artfi_subj_chan_tr{1}{2} = [35,50,54,55,56,82,83,91,124,125];
% artfi_subj_chan_tr{1}{3} = [5,19,42,79,114,124,152,];
% artfi_subj_chan_tr{1}{4} = [5,19,42,79,114,124,152,];
% artfi_subj_chan_tr{4}{1} = [2,4,16,17,26,30,37,44,84,86,95,134,135,151,188];
% artfi_subj_chan_tr{4}{2} = [34,45,50,93,107,108,117,118,124,136,137138,150];
% artfi_subj_chan_tr{4}{3} = [28,30,101,103,118,130,133];
% artfi_subj_chan_tr{5}{1} =  [35,36, 44, 47, 49, 60, 61, 69, 76];
% artfi_subj_chan_tr{5}{2} =  [47,76, 77, 79, 123, 131, 132, 167];
% artfi_subj_chan_tr{5}{5} =  [ 2, 18, 17, 27, 37, 41, 43, 95, 96, 103, 104, 107,117 129, 123, 125, 141];
% artfi_subj_chan_tr{5}{7} = [41,42, 43, 56, 55, 57, 58, 63, 66,69, 79, 84,86,89, 90, 93, 94, 125, 127, 130, 156, 157];
% artfi_subj_chan_tr{5}{8} = [9,79,80,86,106,108,133];
% artfi_subj_chan_tr{5}{9} = [9, 20, 21, 79, 80, 118];
% artfi_subj_chan_tr{6}{1} = [32,34,35,62,69,70,71,72,73, 74,75,77,78,81,83,84,86,95,99,101,105,110,111,113,114,137,139,160,162,165,171,174,198,205,210,...
%                             212,213,214,216,218,223,224,233,234,235,236,237,243,245,246,]; 
% artfi_subj_chan_tr{6}{2} = [19,20,21,23,26,39,55,56,57,58,59,60,63,68,74,77,78,85,87,88,89,91,92,93,95,96,111,112,114,115,118,119,126,127,170,174,175,...
%                             1776,177,184,188,193,204,206,210,213,214,218,219,223,224,226,233,234,236,237]; 
% artfi_subj_chan_tr{6}{3} = [2,19,20,21,25,26,29,33,34,35,36,37,38,39,51,52,53,59,63,72,77,78,80,94,96,103,106,108,133,132,134,135,137,138,139,144,146,...
%                             153,155,157,163,175,200,203,205,212,215,223,228,230,231,232,233,235,236,237,251,282,288,302,303,304,310,311,312]; 
% artfi_subj_chan_tr{6}{4} = [17,19,21,22,23,24,30,31,34,46,54,68,76,90,94,111,114,117,118,122,128,130,131,133,135,136,138,141,142,143,146,147,173,175,...
%                             178,188,191,195,204,214,220,223,229,232,237,238,240,243,244,255,269,274,275,279,281,284,287,290,297,300,302,303,305,307,308,326,330,]; 
% artfi_subj_chan_tr{7}{1} = [59,61,90,93,102,103,109,131,134]; 
% artfi_subj_chan_tr{8}{2} = [40,65]; 
% artfi_subj_chan_tr{10}{1} = [9,10,14,15,17,18,40,41];
% artfi_subj_chan_tr{10}{2} = [12,13,20,22];
% artfi_subj_chan_tr{10}{3} = [2,3,4,16,18,19,24];
% artfi_subj_chan_tr{10}{4} = [5,6,7,16,17,24,26,29,31,34,37,38,39,40,41,42,45];
% artfi_subj_chan_tr{10}{5} = [11,13,27,29,33,36,40,42,56,57,59,64,66];
% artfi_subj_chan_tr{10}{6} = [1,2,3,4,6,7,8,9,10,11,13,15,22,23,24,26,27,29,31,35,36,38,39,40,41,42,43]; % Remove
% artfi_subj_chan_tr{10}{7} = [3,4,8,9,14,25,27,28,31,32,35];
% artfi_subj_chan_tr{10}{8} = [8,12,13,16,20,21,23,25,26,28,29,31,37,38,39,41,42,43,44,45,46,47,48,49];
% artfi_subj_chan_tr{11}{2} = [11,17,27,28,42,48,49,51,52,80,92,93,94,95,96,97,98,137,141,142,153,154,155,156,157,158,186,187,188];
% artfi_subj_chan_tr{11}{1} = [15,10,16,67,70,96,112,114,116,118,134,135,142,1143,147,148,164,168,173174,192,193];
% artfi_subj_chan_tr{11}{3} = [27,33,54,62,82,93,106,107,108,110,111,142,143,144,153,163,164,165,166,192];
% artfi_subj_chan_tr{12}{1} = [58,123,129,141,143,150];
% artfi_subj_chan_tr{12}{2} = [6,13,69,97,113,];
% artfi_subj_chan_tr{12}{3} = [25,26,37];
% artfi_subj_chan_tr{12}{4} = [28,29,64,73];
% artfi_subj_chan_tr{12}{5} = [15,16,19,100,];
% artfi_subj_chan_tr{13}{3} = [16,20,23,24,90,91,151];
% artfi_subj_chan_tr{14}{2} = [11,14,39,90];
% artfi_subj_chan_tr{14}{4} = [48,131];
% artfi_subj_chan_tr{14}{6} = [64];
% artfi_subj_chan_tr{17}{1} = [4,11,12,13,14,16,17,19,20,24,27,33,39,41,50,60,63,64,67,69,71];
% artfi_subj_chan_tr{17}{2} = [3,4,5,12,14,22,24,25,27,29,31,33];
% artfi_subj_chan_tr{18}{2} = [74,76,104,109];
% artfi_subj_chan_tr{18}{3} = [40,42,81,103,105];
% artfi_subj_chan_tr{19}{2} = [57,58,59,61,6263,64,67];
% 

% Load bad Channels: 
ut_bad_channels

%%% ---- PATHS ----- %%%

fieldtripPath = 'C:\Users\barak\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118';  % <- zmień na rzeczywistą ścieżkę
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions')
dataFolder = "D:\Documents_Dell\Predictive_Ripples_2025\preprocessed";
saveFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.ripple_detection_scheme_name);
saveDataFolder = fullfile(dataFolder,params.ripple_detection_scheme_name);
path_ripple_mean_save = fullfile(saveFolder, "RippleChars_XXX_.ps");
path_artif_mean_save = fullfile(saveFolder, "ArtifChars_XXX_.ps");
mkdir(saveFolder); mkdir(saveDataFolder);

% Get file and subject names
files = {dir(fullfile(dataFolder,params.artifact_reject_scheme_name, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractBefore(x, '_artif_rej'), '.mat'), files, 'UniformOutput', false);


%%% Create a result data structs:
all_hipothetical_ripples_cluster = cell(1,length(subjNames)); 
all_timecourseSWR_trial= cell(1,length(subjNames)); 
all_ripplePeaks = cell(1,length(subjNames)); 
all_subjNames = cell(1, length(subjNames));
all_chanNames = cell(1, length(subjNames));


%%% Save parameters in a log file:
T_params = struct2table(params, 'AsArray', true);
logFile = fullfile(saveFolder, 'processing_parameters.txt');
writetable(stack(T_params, 1:width(T_params)), logFile, ...
    'WriteVariableNames', false, 'Delimiter', '\t');


%%% Pre-make PS files
psPaths = {fullfile(saveFolder,"ArtifChars_Hist_.ps"),fullfile(saveFolder,"ArtifChars_FFT_.ps"),fullfile(saveFolder,"ArtifChars_Power_.ps"),...
           fullfile(saveFolder,"ArtifChars_Ripple_.ps"),fullfile(saveFolder,"RippleChars_FFT_.ps"),fullfile(saveFolder,"RippleChars_FFT_.ps"),...
           fullfile(saveFolder,"RippleChars_FFT_.ps"),fullfile(saveFolder,"RippleChars_Hist_.ps"),fullfile(saveFolder,"RippleChars_Power_.ps"),fullfile(saveFolder,"RippleChars_Ripple_.ps")};
parfor l = 1:length(psPaths)
    figure(); print(gcf,psPaths{l}); % Generating PS pre-files for averages
end



%%% ---- ANALYSIS ----- %%%
for s = 1:length(all_subjNames)

    %%% load and unpack the data structures
    data = load(fullfile(dataFolder,params.artifact_reject_scheme_name,files{s}));
    currSubjName = subjNames{s};
        
    data_eeg = data.short_data;
    data_artifact = data.data_artifact;
    artifacts_byChan = data.artifacts_byChan;
    iqr_byChan = data.iqr_byChan;
    range_byChan = data.range_byChan; %(setdiff(1:length(data.range_byChan),data.badChans))
    ieds_byChan = data.ieds_byChan;
    good_trials = data.good_trials;
    
    % Get the size of Channel based data for a given subject
    all_hipothetical_ripples_cluster{s} = cell(1,length(data_eeg.label));
    all_timecourseSWR_trial{s} = cell(1,length(data_eeg.label));
    all_ripplePeaks{s} = cell(1,length(data_eeg.label));
    all_chanNames{s}  = cell(1,length(data_eeg.label)); 
    all_subjNames{s} = currSubjName;

    %%% get Ripple Envelope: 
    [data_ripples,data_viz] = func_rippleband_filtering(data_eeg,params);
    
    %%% Get Analytical Signals used for establishing thresholds
    data_ripples_zscore = data_ripples;
    data_artif = data_ripples;
    
    all_rippledata_viz = zscore(cell2mat(data_ripples.trial),0,2); % Zscoring the data 
    for tr = 1:length(data_artif.trial)
        data_ripples_zscore.trial{tr} = all_rippledata_viz(:,(tr-1)*length(data_artif.trial{1}) +1:tr*length(data_artif.trial{1}));
        for chan = 1:length(data_artif.label)
            data_artif.trial{tr}(chan,isnan(data_artifact.trial{tr}(chan,:))) = nan;
        end
    end
    

    %%% Compute mean and std for thresholding
    joint_artif_data = cell2mat(data_artif.trial);
    std_chan = std(joint_artif_data',1,"omitnan");
    mean_chan = mean(joint_artif_data',1,"omitnan");


    %%% Creating Folder and file structure for saving information
    mkdir(fullfile(saveFolder, currSubjName));
    for channel = 1:length(data_eeg.label)
        psDir   = fullfile(saveFolder, currSubjName,sprintf('%s_trial_ripple.ps',data_eeg.label{channel}));                    % Trial Ripple Detection Viz PS
        path_ripple_save = fullfile(saveFolder,currSubjName, "single_ripples_"+data_eeg.label{channel}+".ps");                 % Individual ripples per channel PS
        path_artif_save = fullfile(saveFolder,currSubjName, "single_artifs_"+data_eeg.label{channel}+".ps");                   % Individual Artifacts per channel PS

        figure(); print(gcf,psDir); % Generating PS pre-files for averages
        figure(); print(gcf,path_ripple_save); % Generating PS pre-files for averages
        figure(); print(gcf,path_artif_save); % Generating PS pre-files for averages
    end


    
    for channel = 1:length(data_eeg.label)


        if ismember(data_eeg.label{channel},artfi_subj_chan{s})
            continue
        end

        %%% --- RIPPLE DETECTION --- %%%
        logFilePath = fullfile(saveFolder,currSubjName,sprintf('%s_ripple_detection_log.txt',data_eeg.label{channel}));
        [hipothetical_ripples_cluster] = func_RippleDetection_main(data_ripples,data_viz,std_chan,mean_chan,artifacts_byChan,channel,params,logFilePath,fullfile(saveFolder,currSubjName),data_eeg.trial,data_eeg.label{channel});
        fprintf("     Channel: %s",data_eeg.label{channel})
    
        % Selecte channel paths
        path_artif_mean_indiv = fullfile(saveFolder, currSubjName, sprintf("%s_ArtifChars_XXX_.png",data_eeg.label{channel}));   % mean indvidual artifact PNGs
        path_ripple_mean_indiv = fullfile(saveFolder, currSubjName, sprintf("%s_RippleChars_XXX_.png",data_eeg.label{channel})); % mean indvidual Ripple PNGs
        psDir   = fullfile(saveFolder, currSubjName,sprintf('%s_trial_ripple.ps',data_artif.label{channel}));                    % Trial Ripple Detection Viz PS
        path_ripple_save = fullfile(saveFolder,currSubjName, "single_ripples_"+data_artif.label{channel}+".ps");                 % Individual ripples per channel PS
        path_artif_save = fullfile(saveFolder,currSubjName, "single_artifs_"+data_artif.label{channel}+".ps");                   % Individual Artifacts per channel PS
        rpower_path = fullfile(saveFolder,currSubjName,sprintf("%s_ripplePowerProfile.ps",data_eeg.label{channel}));              % power Pea visualisation per channel PS

        %%% Pre-alocating all infornation
        powerSWR_all      = [];timecourseSWR_all = []; FFT_SWR_all       = [];
        powerARF_all      = [];timecourseARF_all = []; FFT_ARF_all       = [];
        timecourseSWR_trial = cell(1,length(data_ripples.trial));
        ripplePeaks = cell(1,length(data_ripples.trial));
        artifs_list = [];

        %%% Manual Correction by checking ripple candidate numbers:
        ripcount = 0;


        for it = 1:length(data_ripples.trial)

            %%% Manual Correction
            cluste = 1;
            while cluste <= length(hipothetical_ripples_cluster{it})
                ripcount = ripcount + 1;

                if ~isempty(artfi_subj_chan_tr{s}); if length(artfi_subj_chan_tr{s}) >=channel
                    if ismember(ripcount, artfi_subj_chan_tr{s}{channel})
                        hypothetical_temp = hipothetical_ripples_cluster{it};
                        hypothetical_temp(cluste) = [];
                        hipothetical_ripples_cluster{it} = hypothetical_temp;
                        continue
                    end; end
                end
                cluste = cluste + 1;
            end


            %%% Compute Ripple peak: 
            I = []; 
            for clust  = 1:length(hipothetical_ripples_cluster{it})

                [~,It] = max(data_viz.trial{it}(channel, hipothetical_ripples_cluster{it}{clust}));
                
                I(clust) = It + hipothetical_ripples_cluster{it}{clust}(1);

            end
            ripplePeaks{it} = I;

         
    
            fprintf("\n     Trial N: %d",it)
            
            % Getting per-trial artifact
            idx = artifacts_byChan{channel}{it}; mask = false(1, max(idx)); mask(idx) = true;
            stats = regionprops(logical(mask), 'PixelIdxList', 'Area');
            artifactsClustered = {stats([stats.Area] >= 5).PixelIdxList};
            artifs_list = [artifs_list artifactsClustered];

            %%% --- RIPPLE CHARACTERISTIC COMPUTATION --- %%%

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

       


            %%% --- VISUALISATION --- %%%

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
    
            
        end
    
        %%% (4) Mean ripple characteristics     result = [hipothetical_ripples_cluster{:}];

        plot_rippleSignal(timecourseSWR_all,mean(powerSWR_all,1,"omitnan"),FFT_SWR_all,diff(cellfun(@median,[hipothetical_ripples_cluster{:}])),...
                        size(timecourseSWR_all,1),path_ripple_mean_save,params,currSubjName+" // "+data_artif.label{channel},path_ripple_mean_indiv)

        %%% (5) Mean Artif characteristics
        plot_rippleSignal(timecourseARF_all,mean(powerARF_all,1,"omitnan"),FFT_SWR_all,diff(cellfun(@median,artifs_list)),...
                        size(timecourseARF_all,1),path_artif_mean_save,params,currSubjName+" // "+data_artif.label{channel},path_artif_mean_indiv)

    
        % Transform PS into PDFs
        if params.vis; system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); end
        if params.vis > 1; system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_ripple_save,".ps", ".pdf "),path_ripple_save));  delete(path_ripple_save); end
        if params.vis == 3; system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(path_artif_save,".ps", ".pdf "),path_artif_save));  delete(path_artif_save); end
        if params.vis > 1; system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(rpower_path,".ps", ".pdf "),rpower_path));  delete(rpower_path); end


        %%% Save current channel data: 
        all_hipothetical_ripples_cluster{s}{channel} = hipothetical_ripples_cluster;
        all_timecourseSWR_trial{s}{channel} = timecourseSWR_trial;
        all_ripplePeaks{s}{channel}  =  ripplePeaks{it};
        all_chanNames{s}{channel} = data_eeg.label{channel};

    end
end

% Transform ALL Per-channel Ripple and Artif characteristics: 
psPath = fullfile(saveFolder,"ArtifChars_FFT_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Hist_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Power_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"ArtifChars_Ripple_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_FFT_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_FFT_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Hist_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Power_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

psPath = fullfile(saveFolder,"RippleChars_Ripple_.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath));  %delete(psPath)

%%% Final Saving of the data 
save(fullfile(saveDataFolder,"rippleData.mat"),"all_hipothetical_ripples_cluster","all_timecourseSWR_trial","all_ripplePeaks","all_chanNames","all_subjNames")










