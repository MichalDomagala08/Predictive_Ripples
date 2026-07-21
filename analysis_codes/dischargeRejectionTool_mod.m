function [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool_mod(patientID,sessionNum,data_norm,plotFlag)
%DISCHARGREJECTIONTOOL Tool for rejecting trials with epileptiform discharges.
%
%   [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool
%       (patientID,sessionNum,nevFile) finds outliers based on their
%       maximum and mimimum values.
%
%   optional input argument: plotFlag == 0 will not plot rejected trials.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)

% default argument for plotting.
if nargin==3
    plotFlag = 1;
end




%% timing (seconds)
pre = 4;
post = 10;


%% parsing channel labels
macroLabels = data_norm.label;
numBFs = length(macroLabels);
% BFlabels =



%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% initializing LFPmat
LFPlabels = macroLabels;
LFPmat = zeros(size(data_norm.trial{1},2),length(data_norm.trial),size(data_norm.trial{1},1));
% for loop to save multiple epochs

    
for ch = 1:length(LFPlabels)
    offset = 3 * std(LFPmat(:,:,ch), [], 'all');  % automatycznie dopasowany do amplitudy
    if plotFlag
        % figure window.
        ah = figure('Color',[0 0 0]);
        hold on
    end
    
    for tt = 1:length(data_norm.trial)
        
        
        %% [20160622] constructing LFP tensor.
        LFPmat(:,tt,ch) = data_norm.trial{tt}(ch,:);
 
        
        %% [20160622] time vectors
        tsec = data_norm.time{tt};
        
        if plotFlag
            %% [20160622] plotting LFP from each trial brushing to reject
            plot(tsec,LFPmat(:,tt,ch)+((tt-1)*offset),'color',rgb('white'))
            axis tight off
        end
        
    end
    
%         keyboard 
    
    [rejectTrialsperChannel(ch).rejectTheseTrials,iprange,fence] = outliers(squeeze(range(LFPmat(:,:,ch))));
    rejectTrialsperChannel(ch).channelLabel = deblank(LFPlabels{ch});
    
    % saving all of the bad trials
    if isequal(ch,1)
        trialsAcrossChannels = outliers(squeeze(range(LFPmat(:,:,ch))));
    else
        trialsAcrossChannels = cat(2,trialsAcrossChannels,outliers(squeeze(range(LFPmat(:,:,ch)))));
    end
    
    %%
    if plotFlag
        for rjct = rejectTrialsperChannel(ch).rejectTheseTrials
            %% [20160622] plotting LFP from each trial brushing to reject
           if rjct ~= 0 
                plot(tsec,LFPmat(:,rjct,ch)+((rjct-1)*offset),'color',rgb('springgreen'))
           end
        end
        set(gcf, 'WindowState', 'maximized');
        title(sprintf('Channel: %s. || Rejected Trials in Green',rejectTrialsperChannel(ch).channelLabel),'color',rgb('springgreen'))
        
        
        %% saving figures.
        display('saving to Elliot"s dropbox. Thanks!')
        dbPath = fullfile('D:\Documents_Dell\Predictive_Ripples_2025\results\ArtifactRejection',patientID);
        mkdir(dbPath)
        fName = sprintf('%s/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
        saveas(gcf,fName, 'pdf')
        close(gcf)
   
        
        display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
        hold off
        
    end
end % looping over lfp channels

trialsAcrossChannels = unique(trialsAcrossChannels);
