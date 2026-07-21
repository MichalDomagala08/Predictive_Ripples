
%%% Parameters:

% Ripple Detection Parameters:


% Other Generation parameters: 

% - ERP timestamps 

params.presentation = 2001; % 500 *4s +1
params.erpEnd = 0.75;    % in seconds relative to presentation time 
params.erpBeg = -0.25; % in seconds relative to presentation time 
params.hipOnly = true;

% Paths
codeFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions';
dataPreprocessedFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\preprocessed';
dataFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\data\reref'; %Folder where preprocessed Electrophysiological data are bieng housed
saveFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\results'; %Folder where preprocessed Electrophysiological data are bieng housed
imageFolder =    'D:\Documents_Dell\Predictive_Ripples_2025\data\electrodes_viz';
addpath("functions")

warning('off', 'MATLAB:printerOutput:PostScriptNotSupported');

%%% Preprocessing and Preparing 



% Data 
files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractAfter(x, 'sub'), '.mat'), files, 'UniformOutput', false);



% Electrode info
T = readtable('D:\Documents_Dell\Predictive_Ripples_2025\data\Electrodes.xlsx','ReadRowNames',true, 'Sheet', 'electrodes');

if params.hipOnly
    electrodeFields = {'antHP',	'postHP'};
else
    electrodeFields = {'ITC',	'Vx', 'antHP',	'postHP'};
end

result_channels = {};
labels = {};

warning('off', 'all');
%%% Create ERP Plots for Every channel! 
parfor iter = 1:length(subjNames)
    currFile = load(fullfile(dataFolder,files{iter}),'data').data;
    currSubjName = subjNames{iter};
    labels{iter} = currFile.data_eeg.label; % Load labels to a different structure

    % W pętli subjectów, przed rozpoczęciem przetwarzania kanałów:
    fprintf('Subject: %s ...', currSubjName);
 

     % Load Image info:
    imageStructs = dir(fullfile(imageFolder, currSubjName, '*.jpg'));
    images = {imageStructs.name}';  % tylko nazwy plików    
    %%
    fprintf('\n ...Current Channel:              ');
    for lab = 1:length(currFile.data_eeg.label)

        % W pętli kanałów:
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f/%3.0f', lab, length(currFile.data_eeg.label));
        drawnow limitrate;

        plot_rippleERP_parfor(lab,currFile.data_eeg.label{lab},currSubjName,currFile.data_eeg.trial,images,params,saveFolder,imageFolder)
        fprintf(' Done\n');
    end 
    %%
end

%%
%%% Transform PS to PDFS for each patient; get hippocampal channel candidate names 
for iter = 13:length(subjNames)
    currSubjName = subjNames{iter};

      
    goodChannels = [];
    for i = 1:length(electrodeFields)
        channel_string = T.(electrodeFields{i}){currSubjName};
        if ~isempty(channel_string)
            goodChannels = [goodChannels; split(channel_string,',')];
        end
    end
    
    % For Each Subject Get the Candidate Hippocampal Channel: Verify Later 
    if ~isempty(goodChannels)
        result_channels{iter} = labels{iter}(cellfun(@(x) any(ismember(strsplit(x, '-'), strtrim(goodChannels))), labels{iter}));
    end
    
    psPath = fullfile(saveFolder,'ERPs',sprintf('ERP_%s.ps',currSubjName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psPath,".ps", ".pdf "),psPath))
    %delete(psPath)

end

%%


function plot_rippleERP_parfor(chanNum, currChannelName, currSubjName, trialData, images, params, saveFolder, imageFolder)
    % A Quick function that plots Ripples and is catered for Parfor

    
    origChannelsCurr = split(currChannelName,'-');
    
    % Extract trial data: Creating Trial x Time for each channel 
    tmp = cellfun(@(x) x(chanNum,:), trialData, 'UniformOutput', false);
    trialTime = vertcat(tmp{:});
    
    % Find images - by findingi Slices. Getting only First. Getting 2 images for both channels
    pattern1 = string(regexprep(origChannelsCurr{1}, '([A-Za-z])(\d)', '$1_$2')) + "(Slices|\.)";
    matchingIdx1 = find(~cellfun('isempty', regexp(images, pattern1, 'once')), 1, 'first');
    
    pattern2 = string(regexprep(origChannelsCurr{2}, '([A-Za-z])(\d)', '$1_$2')) + "(Slices|\.)";
    matchingIdx2 = find(~cellfun('isempty', regexp(images, pattern2, 'once')), 1, 'first');
    
    if isempty(matchingIdx1) || isempty(matchingIdx2); return;  end % Pomijamy kanały bez obrazów
    
    imgName1 = images{matchingIdx1};
    imgName2 = images{matchingIdx2};
    
    % CHECK: if it is string and Not  Cell arrays! 
    if iscell(imgName1); imgName1 = imgName1{1}; end
    if iscell(imgName2); imgName2 = imgName2{1}; end
    
    % Create figure
    fig = figure('Visible','off','Units','normalized','Position',[0 0 0.7 1]);
    fig.PaperOrientation = 'portrait';
    fig.PaperUnits = 'normalized';
    fig.PaperPosition = [0 0 1 1];
    
    % Plot ERP
    subplot(3,1,1);
    try
        stdshade(trialTime(:,params.presentation+params.erpBeg*500:params.presentation+params.erpEnd*500),...
            0.5,[0.8500 0.3250 0.0980]);
    catch
        % Ignoruj błędy w stdshade
    end
    xticks([1,126,251,376,501,626,751]);
    xticklabels([-0.5,-0.25,0,0.25,0.5,0.75,1]);
    xline(251);
    xlabel("time (s)");
    ylabel("Voltage (uV)");
    title(sprintf("ERP and Locs: %s",currChannelName));

    % Plot images - TERAZ BEZ PODWÓJNEGO CIĄGU
    subplot(3,1,2);
    imgPath1 = fullfile(imageFolder, currSubjName, char(imgName1));  % char() konwertuje do stringa
    if exist(imgPath1, 'file')
        try
            imshow(imread(imgPath1));
            title(origChannelsCurr{1});
        catch
            title(sprintf('%s (błąd)', origChannelsCurr{1}));
        end
    else
        title(sprintf('%s (brak)', origChannelsCurr{1}));
    end
    
    subplot(3,1,3);
    imgPath2 = fullfile(imageFolder, currSubjName, char(imgName2));  % char() konwertuje do stringa
    if exist(imgPath2, 'file')
        try
            imshow(imread(imgPath2));
            title(origChannelsCurr{2});
        catch
            title(sprintf('%s (błąd)', origChannelsCurr{2}));
        end
    else
        title(sprintf('%s (brak)', origChannelsCurr{2}));
    end
    
    % Save
    erpDir = fullfile(saveFolder,'ERPs');
    if ~exist(erpDir, 'dir'), mkdir(erpDir); end
    print(fig, fullfile(erpDir, sprintf('ERP_%s.ps',currSubjName)), '-dpsc', '-append', '-r300');
    close(fig);
end


function plot_rippleERP_old(chanNum,currChannelName,currSubjName,trialData,images,params,saveFolder,imageFolder)

    origChannelsCurr = split(currChannelName,'-');

    % Cut out Matrix with Trials x Time for a given Channel 
    tmp = cellfun(@(x) x(chanNum,:),  trialData, 'UniformOutput', false);
    trialTime = vertcat(tmp{:});
    
    % Find a relevant image
    pattern = string(regexprep(origChannelsCurr{1}, '([A-Za-z])(\d)', '$1_$2')) + "(Slices|\.)";
    idx1 = ~cellfun('isempty', regexp(images, pattern, 'once'));
    
    pattern = string(regexprep(origChannelsCurr{2}, '([A-Za-z])(\d)', '$1_$2')) + "(Slices|\.)";
    idx2 = ~cellfun('isempty', regexp(images, pattern, 'once'));
    
    % Plot Image
    fig = figure('Visible','off','Units','normalized','Position',[0 0 0.7 1]);   % tall portrait window
    
    fig.PaperOrientation = 'portrait'; % Make it portrait for printing
    fig.PaperUnits = 'normalized';
    fig.PaperPosition = [0 0 1 1];
    
    
    subplot(3,1,1);
    stdshade(trialTime(:,params.presentation+params.erpBeg*500:params.presentation+params.erpEnd*500),0.5,[0.8500 0.3250 0.0980]);
    xticks([1,126,251,376,501,626,751]);
    xticklabels([-0.5,-0.25,0,0.25,0.5,0.75,1]);
    xline(251);
    xlabel("time (s)");
    ylabel("Voltage (uV)");
    title(sprintf("ERP and Locs: %s",currChannelName));

    subplot(3,1,2);
    imshow(imread( fullfile(imageFolder, currSubjName, images{idx1})));
    title(origChannelsCurr{1});
    
    subplot(3,1,3);
    imshow(imread(fullfile(imageFolder, currSubjName, images{idx2})));
    title(origChannelsCurr{2});
    
    
    
    % Save as PostScript page
    if ~exist(fullfile(saveFolder,'ERPs'), 'dir'); mkdir(fullfile(saveFolder,'ERPs')); end

    print(fig, fullfile(saveFolder,'ERPs',sprintf('ERP_%s.ps',currSubjName)), '-dpsc',  '-append', '-r300');
    close(fig);
end

