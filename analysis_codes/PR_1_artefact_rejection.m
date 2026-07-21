
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -- ARTIFACT REJECTION SCRIPT -- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script Rejects artfiacts n the data usng several modules:
% - Spike Artifacts: Using Diamond et al., algorithm to detect prominent IEDs
% - Range Artifacts: Using diff on filtered, bha and absolute signal to 
%                    account for abnormal and rapidpeaks and troughs
% - IQR with MAD - to account for long and rapid 
% - Absolute 9sstd threshold: to account for large deflections: 
%
% This script also visualise each artif type and saves the filtered data as
% well as artifact timestamp per subject - channel - trial



% params
params.hipOnly = true;

%%% Spike Detection Parameters %%%
params.spikeCtsThresh  = 2;      % Zawsze 1 na czas debugowania!
params.spikePeakWin    = 0.3;   % Zwiększ okno dopasowania (hp) do 150ms
params.spikeZThresh    = 4;    % Obniż próg (z-score na dużych danych rzadko dobija do 4 dla rozlazłych fal)
params.spikeAmpScale   = 2.5;    % Obniż sumaryczny wymóg Peak-to-Trough
params.spikeMNegPeakW  = 150;    % KLUCZ: Pozwól negatywnej fazie trwać do 300ms
params.spikeTrackPeaks = false;   % Szukaj od pozytywnego (tak jak na Twoim obrazku)
params.spikeWindow     = 250;    % Zwiększ margines wycinania artefaktu wokół IED

% Artifact Rejection parameters:
params.artPadding = 70; % how many samples should we pad the artifacts with
params.cluster_tolerance = 20;  % How many samples is it "close" 
params.range_threshold = 6; % How many stds from mean is considered an artifact 
params.iqr_w  = 3; %How many IQR do we have to surpass to have an artefact
params.artif_mad = 8;
analysisName = 'ArtifactRejection_alternative2';

% Paths
codeFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions';
dataPreprocessedFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\preprocessed';
dataFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\data\reref'; %Folder where preprocessed Electrophysiological data are bieng housed
saveFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\results'; %Folder where preprocessed Electrophysiological data are bieng housed
imageFolder =    'D:\Documents_Dell\Predictive_Ripples_2025\data\electrodes_viz';
addpath("functions")

warning('off', 'MATLAB:printerOutput:PostScriptNotSupported');
mkdir(fullfile(saveFolder,analysisName))
mkdir(fullfile(dataPreprocessedFolder,analysisName))


%%% Save parameters in a log file:

T_params = struct2table(params, 'AsArray', true);
logFile = fullfile(saveFolder, analysisName, 'processing_parameters.txt');
writetable(stack(T_params, 1:width(T_params)), logFile, ...
    'WriteVariableNames', false, 'Delimiter', '\t');


%%% Preprocessing and Preparing 


% Electrode info
T = readtable('D:\Documents_Dell\Predictive_Ripples_2025\data\Electrodes.xlsx','ReadRowNames',true, 'Sheet', 'responsive');
if params.hipOnly
    electrodeFields = {'antHP',	'postHP'};
else
    electrodeFields = {'ITC',	'Vx', 'antHP',	'postHP'};
end



% Data 
files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractAfter(x, 'sub'), '.mat'), files, 'UniformOutput', false);

goodChannels = {};

fprintf("\n------------------------------------------\n")
disp("   ==== Artifact Rejection ====")
fprintf("------------------------------------------\n\n\n")

for iter = 1:length(subjNames)


    %%%% --- READ single subjeft data ---- %%%%
    raw_data = load(fullfile(dataFolder,files{iter}),'data').data.data_eeg;
    currSubjName = subjNames{iter};
    

    diary(fullfile(saveFolder,analysisName, sprintf('%s_ar_log.txt', currSubjName)));
    diary on
    disp("--------------------------------------------------------")
    fprintf("CURRENT SUBJECT: %s \n\n",currSubjName)
    

    
    % Select Rippleband channels only:   
    channel_strings = cellfun(@(f) T.(f){currSubjName}, electrodeFields, 'UniformOutput', false);
    channel_strings = channel_strings(~cellfun('isempty', channel_strings));  % usuń puste
    splits          = cellfun(@(s) split(s,','), channel_strings, 'UniformOutput', false);
    goodChannels{iter} = vertcat(splits{:});

      % For Each Subject Get the Candidate Hippocampal Channel: Verify Later
      % (THis is just a placeholder! When I'll get channel pairs as good
      % channels THEN i will patch it!)
    if isempty(goodChannels{iter})
    %     final_goodChannels{iter} = raw_data.label(cellfun(@(x) any(ismember(strsplit(x, '-'), strtrim(goodChannels))), raw_data.label));
    % else
        continue % for situations when there are NO hippocampal channels ready! 
    end


    
    
    % Select only Hip Channels:
    chanHipIdx = ismember(raw_data.label,goodChannels{iter});

    % CHECK (We have to have the same number of raw_data.label as in
    % goodChans)

    if sum(chanHipIdx) ~= length(goodChannels{iter})
        error(sprintf("A channel got lost in translation: %s",strjoin(setdiff(goodChannels{iter},raw_data.label), ' ')))
    end
    raw_data.label = raw_data.label(chanHipIdx);
    for tr = 1:length(raw_data.trial)
        raw_data.trial{tr} = raw_data.trial{tr}(chanHipIdx,:);
    end
    
   %[artifacts_byChan,ieds_byChan,iqr_byChan,range_byChan,data_artifact,ied_timestamps] =  func_artifact_rejection_main(raw_data,params);
    [artifacts_byChan,ieds_byChan,iqr_byChan,range_byChan,data_artifact,ied_timestamps] =  func_artifact_rejection_alternative(raw_data,params);

    short_data = raw_data;

    diary oFf

    %%% Cut first and last 1s from trials (to avoid edge artefacts) & Cut out Bad channels 
    % data_artifact.label = data_artifact.label(setdiff(1:size(data_artifact.trial{1},1),badChans));
    % short_data.label = short_data.label(setdiff(1:size(short_data.trial{1},1),badChans));
    % 
    % for i = 1:length(data_artifact.trial)
    %     data_artifact.trial{i} = data_artifact.trial{i}(:,1001:length(data_artifact.trial{i})-1000);
    %     short_data.trial{i} = short_data.trial{i}(:,1001:length(short_data.trial{i})-1000);
    % 
    %     data_artifact.time{i} = data_artifact.time{i}(1001:length(data_artifact.time{i})-1000);
    %     short_data.time{i} = short_data.time{i}(1001:length(short_data.time{i})-1000);
    % end
    % 
    
    %%% Find trials that had more than 10% of time artefactual in ANY OF THE
    %%% CHANNELS (for later rejection)
    good_trials = cell(1,length(artifacts_byChan));
    for cells = 1:length(artifacts_byChan)
        currChan = artifacts_byChan{cells};
        good_trials{cells} =  cellfun(@length,currChan)/7001 < 0.1;
    end

    %%% Save Data 
    save(fullfile(dataPreprocessedFolder,analysisName,sprintf("%s_artif_rej.mat",currSubjName)),"short_data","data_artifact","good_trials","artifacts_byChan","ieds_byChan","iqr_byChan","range_byChan")
    

    %% --- VISUALIZE---- %%

    %%% Plot Artifact Rejected Data for trials 
    
    if ~exist(fullfile(saveFolder, analysisName, currSubjName))
        mkdir(fullfile(saveFolder, analysisName, currSubjName))
    end
    
    parfor chan = 1:length(raw_data.label)
        psFile = fullfile(saveFolder, analysisName, currSubjName, sprintf('%s_artifacts.ps', raw_data.label{chan}));

        %%
        for i = 1:numel(raw_data.trial)
            figure("Visible","off",'Units','normalized','Position',[0 0 0.7 1],'PaperOrientation','Landscape'); clf; hold on
            plot(raw_data.time{i}, raw_data.trial{i}(chan,:), 'k')

            %%% IED artefact Coloring
            art = sort(ieds_byChan{chan}{i}(:)');
            if ~isempty(art)
                yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
                for k = 1:numel(s)
                    patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'r','FaceAlpha',.3,'EdgeColor','none')
                end
            end

            %%% IQR artefact Coloring
            art = sort(iqr_byChan{chan}{i}(:)');
            if ~isempty(art)
                yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
                for k = 1:numel(s)
                    patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'g','FaceAlpha',.3,'EdgeColor','none')
                end
            end

            %%% Range artefact Coloring
            art = sort(range_byChan{chan}{i}(:)');
            if ~isempty(art)
                yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
                for k = 1:numel(s)
                    patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'b','FaceAlpha',.3,'EdgeColor','none')
                end
            end
            plot(raw_data.time{i}, raw_data.trial{i}(chan,:), 'k');

            title(sprintf('%s | Trial %d',raw_data.label{chan},i))
            if i==1; print(gcf,psFile,'-dpsc'); else; print(gcf,psFile,'-dpsc','-append'); end
            close(gcf);
        end
        %%
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psFile,".ps", ".pdf "),psFile));
        delete(psFile);
    end

    %%
    % (2) as a matrix of trials 
    figure('Visible','Off','Units','normalized','Position',[0 0 0.7 1],'PaperOrientation','Landscape');
    nChan = length(raw_data.label);
    nCols = ceil(sqrt(nChan));
    nRows = ceil(nChan / nCols);
    
    tl =tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    for chanNum = 1:length(raw_data.label)
        tmp      = cellfun(@(x) x(chanNum,:), raw_data.trial, 'UniformOutput', false);
        trialMat = vertcat(tmp{:});  % [nTrials x nSamples]
    
        nexttile
        imagesc(trialMat);
        colorbar();
        ylabel(raw_data.label{chanNum});
    end
    
    xlabel('Sample'); title(tl, sprintf('%s Trials x Time per Channel',currSubjName));
    exportgraphics(gcf, fullfile(saveFolder, analysisName, currSubjName, sprintf('%s_artifacts_matrix.png', currSubjName)), 'Resolution', 300);
    print(gcf,fullfile(saveFolder, analysisName,'artifacts_matrix.ps'),'-dpsc','-append'); 

    %%
    % (3) trials after artifact rejection (and interpolation for viz purposes


    % NaN interpolation
    data_interp = data_artifact;
    
    for chanNum = 1:numel(data_artifact.label)
        for tr = 1:numel(data_artifact.trial)
            x = data_artifact.trial{tr}(chanNum, :);
            nans = isnan(x);
            if any(nans) && ~all(nans)
                idx = 1:numel(x);
                x(nans) = interp1(idx(~nans), x(~nans), idx(nans), 'linear', 'extrap');
            end
            data_interp.trial{tr}(chanNum, :) = x;
        end
    end

    % Artifact Data Viz:
    figure('Visible','Off','Units','normalized','Position',[0 0 0.7 1],'PaperOrientation','Landscape');
    nChan = length(data_artifact.label); nCols = ceil(sqrt(nChan)); nRows = ceil(nChan / nCols);

    tl =tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    for chanNum = 1:length(data_artifact.label)
        tmp      = cellfun(@(x) x(chanNum,:), data_artifact.trial, 'UniformOutput', false);
        trialMat = vertcat(tmp{:});  % [nTrials x nSamples]
    
        nexttile
        imagesc(trialMat); colorbar();
        ylabel(data_artifact.label{chanNum});
    end
    
    xlabel('Sample'); title(tl, sprintf('%s Trials x Time per Channel',currSubjName));
    exportgraphics(gcf, fullfile(saveFolder, analysisName, currSubjName, sprintf('%s_artifacts_matrix_removed.png', currSubjName)), 'Resolution', 300);
    print(gcf,fullfile(saveFolder, analysisName,'artifacts_matrix_removed.ps'),'-dpsc','-append'); 


     % Intepolated Data Viz:
    figure('Visible','Off','Units','normalized','Position',[0 0 0.7 1],'PaperOrientation','Landscape');
    nChan = length(data_interp.label); nCols = ceil(sqrt(nChan)); nRows = ceil(nChan / nCols);

    tl =tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    for chanNum = 1:length(data_interp.label)
        tmp      = cellfun(@(x) x(chanNum,:), data_interp.trial, 'UniformOutput', false);
        trialMat = vertcat(tmp{:});  % [nTrials x nSamples]
    
        nexttile
        imagesc(trialMat); colorbar();
        ylabel(data_interp.label{chanNum});
    end
    
    xlabel('Sample'); title(tl, sprintf('%s Trials x Time per Channel',currSubjName));
    exportgraphics(gcf, fullfile(saveFolder, analysisName, currSubjName, sprintf('%s_artifacts_matrix_interp.png', currSubjName)), 'Resolution', 300);
    print(gcf,fullfile(saveFolder, analysisName,'artifacts_matrix_interp.ps'),'-dpsc','-append'); 


    disp('Subject_done!')

end
%%
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(fullfile(saveFolder, analysisName,'artifacts_matrix_interp.ps'),".ps", ".pdf "),fullfile(saveFolder, analysisName,'artifacts_matrix_interp.ps')));
delete(fullfile(saveFolder, analysisName,'artifacts_matrix_interp.ps'));

system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(fullfile(saveFolder, analysisName,'artifacts_matrix.ps'),".ps", ".pdf "),fullfile(saveFolder, analysisName,'artifacts_matrix.ps')));
delete(fullfile(saveFolder, analysisName,'artifacts_matrix.ps'));

system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(fullfile(saveFolder, analysisName,'artifacts_matrix_removed.ps'),".ps", ".pdf "),fullfile(saveFolder, analysisName,'artifacts_matrix_removed.ps')));
delete(fullfile(saveFolder, analysisName,'artifacts_matrix_removed.ps'));

%% Alt rejection - Faulty, not good enough1
%[rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool_mod(currSubjName,1,raw_data,1);

%%
