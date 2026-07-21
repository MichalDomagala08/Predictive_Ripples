
% params
params.hipOnly = true;

%%% Spike Detection Parameters %%%
%%% Spike Detection Parameters %%%
params.spikeCtsThresh  = 1;      % Zawsze 1 na czas debugowania!
params.spikePeakWin    = 0.3;   % Zwiększ okno dopasowania (hp) do 150ms
params.spikeZThresh    = 3;    % Obniż próg (z-score na dużych danych rzadko dobija do 4 dla rozlazłych fal)
params.spikeAmpScale   = 0.8;    % Obniż sumaryczny wymóg Peak-to-Trough
params.spikeMNegPeakW  = 200;    % KLUCZ: Pozwól negatywnej fazie trwać do 300ms
params.spikeTrackPeaks = false;   % Szukaj od pozytywnego (tak jak na Twoim obrazku)
params.spikeWindow     = 100;    % Zwiększ margines wycinania artefaktu wokół IED

%%% Ripple Detection Parameters %%%
params.lowpassfreq     = 80;  % cutoffs for filtering ripple band
params.highpassfreq    = 120;
params.cluster_join    = 7;   % timepoint distance on which clusters will be joined %15ms                                            
params.maxRippleLength = 100; % Maximum length (in samples of a ripple - Current 100ms)
params.maximum_peak    = 4;   % SDs of maximum peak from mean       
params.num_peaks       = 3;   % minimum peaks that are detected in ripple                                                             
params.detection_peaks = "peaks"; % thrsh - lowerbound threshold; 
                                  % peaks - minimum peak attained; 


% Paths
codeFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions';
dataPreprocessedFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\preprocessed';
dataFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\data\reref'; %Folder where preprocessed Electrophysiological data are bieng housed
saveFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\results'; %Folder where preprocessed Electrophysiological data are bieng housed
imageFolder =    'D:\Documents_Dell\Predictive_Ripples_2025\data\electrodes_viz';
addpath("functions")

warning('off', 'MATLAB:printerOutput:PostScriptNotSupported');

%%% Preprocessing and Preparing 


% Electrode info
T = readtable('D:\Documents_Dell\Predictive_Ripples_2025\data\Electrodes.xlsx','ReadRowNames',true, 'Sheet', 'electrodes');
if params.hipOnly
    electrodeFields = {'antHP',	'postHP'};
else
    electrodeFields = {'ITC',	'Vx', 'antHP',	'postHP'};
end



% Data 
files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractAfter(x, 'sub'), '.mat'), files, 'UniformOutput', false);

final_goodChannels = {};
%%%% --- READ single subjeft data ---- %%%%
iter = 1;
raw_data = load(fullfile(dataFolder,files{iter}),'data').data.data_eeg;
currSubjName = subjNames{iter};


% Select Rippleband channels only:   
goodChannels = [];
for i = 1:length(electrodeFields)
    channel_string = T.(electrodeFields{i}){currSubjName};
    if ~isempty(channel_string)
        goodChannels = [goodChannels; split(channel_string,',')];
    end
end

  % For Each Subject Get the Candidate Hippocampal Channel: Verify Later
  % (THis is just a placeholder! When I'll get channel pairs as good
  % channels THEN i will patch it!)
if ~isempty(goodChannels)
    final_goodChannels{iter} = raw_data.label(cellfun(@(x) any(ismember(strsplit(x, '-'), strtrim(goodChannels))), raw_data.label));
end


% Select only Hip Channels:
chanHipIdx = ismember(raw_data.label,final_goodChannels{iter});
raw_data.label = raw_data.label(chanHipIdx);
for tr = 1:length(raw_data.trial)
    raw_data.trial{tr} = raw_data.trial{tr}(chanHipIdx,:);
end

%%% Filtering to Hilber Freqs! 
[data_ripples,data_viz] = func_rippleband_filtering(raw_data,params);

disp("aa")





dd =  cell2mat(raw_data.trial);
[Ypk,Xpk,Wpk,Ppk] = findpeaks(zscore(dd(5,:)),'MinPeakHeight',3);
[Ypkn,Xpkn,Wpkn,Ppkn] = findpeaks(-zscore(dd(5,:)),'MinPeakHeight',3);

aa = find(abs(zscore(dd(5,:))) > 3);

figure(); plot(dd(5,:))

% Extract trial data: Creating Trial x Time for each channel 
tmp = cellfun(@(x) x(6,:), raw_data.trial, 'UniformOutput', false);
trialTime = vertcat(tmp{:});
set(0, 'DefaultFigureVisible', 'on');
figure();
imagesc(trialTime);
colorbar();

%%
[artifacts_byChan,data_artf,currentSpikes,ied_timestamps,badChans] =  func_artifact_rejection_main(raw_data,params);
short_data = raw_data;


%%% Cut first and last 1s from trials (to avoid edge artefacts) & Cut out Bad channels 
data_artf.label = data_artf.label(setdiff(1:size(data_artf.trial{1},1),badChans));
short_data.label = short_data.label(setdiff(1:size(short_data.trial{1},1),badChans));

for i = 1:length(data_artf.trial)
    data_artf.trial{i} = data_artf.trial{i}(setdiff(1:size(data_artf.trial{i},1),badChans),1001:length(data_artf.trial{i})-1000);
    short_data.trial{i} = short_data.trial{i}(setdiff(1:size(short_data.trial{i},1),badChans),1001:length(short_data.trial{i})-1000);

    data_artf.time{i} = data_artf.time{i}(1001:length(data_artf.time{i})-1000);
    short_data.time{i} = short_data.time{i}(1001:length(short_data.time{i})-1000);
end


%%% Find trials that had more than 25% of time artefactual in ANY OF THE
%%% CHANNELS (for later rejection)
good_trials = cell(1,length(artifacts_byChan));
for cells = 1:length(artifacts_byChan)
    currChan = artifacts_byChan{cells};
    good_trials{cells} =  cellfun(@length,currChan)/7001 < 0.25;
end


%%
%%% Plot Artifact Rejected Data for trials 

if ~exist(fullfile(saveFolder, 'ArtifactRejection', currSubjName))
    mkdir(fullfile(saveFolder, 'ArtifactRejection', currSubjName))
end

% (1) As a timecourse for every chan
parfor chan = 1:length(short_data.label)
    psFile = fullfile(saveFolder, 'ArtifactRejection', currSubjName, sprintf('%s_artifacts.ps', short_data.label{chan}));
    for i = 1:numel(short_data.trial)
        figure("Visible","off"); clf; hold on
        art = sort(artifacts_byChan{chan}{i}(:)');
        art = art(art>1000 & art <6001)-1000;
        if ~isempty(art)
            yl = ylim; gaps = find(diff(art)>1);
            s = art([1,gaps+1]); e = art([gaps,end]);
            for k = 1:numel(s)
                patch(short_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'r','FaceAlpha',.3,'EdgeColor','none')
            end
        end
        plot(short_data.time{i}, short_data.trial{i}(chan,:), 'k')

        title(sprintf('%s | Trial %d',short_data.label{chan},i))
        if i==1; print(gcf,psFile,'-dpsc'); else; print(gcf,psFile,'-dpsc','-append'); end
        close(gcf)
    end
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psFile,".ps", ".pdf "),psFile))
    delete(psFile)
end

%%
% (2) as a matrix of trials 
figure('Visible','On');
nChan = length(short_data.label);
nCols = ceil(sqrt(nChan));
nRows = ceil(nChan / nCols);

tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
for chanNum = 1:length(short_data.label)
    tmp      = cellfun(@(x) x(chanNum,:), short_data.trial, 'UniformOutput', false);
    trialMat = vertcat(tmp{:});  % [nTrials x nSamples]

    nexttile
    imagesc(trialMat);
    colorbar();
    ylabel(short_data.label{chanNum});

    save(fullfile(saveFolder, 'ArtifactRejection', currSubjName, sprintf('%s_artifacts_matrix.ps', currSubjName)))
end

xlabel('Sample'); title(t, 'Trials x Time per Channel');

%%

function plot_review_trial(data_artf)
    i = 1;
    while i >= 1 && i <= 80
        figure(1); clf;
        plot(data_artf.time{i}, data_artf.trial{i})
        legend(data_artf.label);
        title(sprintf('Trial %d / 80', i));
    
        [~, ~, key] = ginput(1);  % czeka na kliknięcie/klawisz
    
        % Odczytaj ostatni klawisz
        key = get(gcf, 'CurrentKey');
        if strcmp(key, 'rightarrow')
            i = min(i + 1, 80);
        elseif strcmp(key, 'leftarrow')
            i = max(i - 1, 1);
        elseif strcmp(key, 'escape')
            break;
        end
    end
end
%% Alt rejection - Faulty, not good enough1
%[rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool_mod(currSubjName,1,raw_data,1);

%%
