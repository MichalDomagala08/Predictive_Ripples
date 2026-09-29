function outfile = preprocess_bipolar(iSubj, datadir, outdir, notch_hw)
% PREPROCESS_BIPOLAR Apply bipolar re-referencing to iEEG data
%
%   preprocess_bipolar(iSubj, datadir, outdir) loads raw iEEG data for
%   subject index iSubj, applies bipolar re-referencing within each
%   electrode, and saves the result.
%   preprocess_bipolar(iSubj, datadir, outdir, notch_hw) also removes line
%   noise (cleanTrialEEG, +-notch_hw Hz at 60/120/180/240 Hz) before saving,
%   so the file is written once, already filtered.
%   Every trial is demeaned per channel before saving (cfg.trial_demean).
%
%   Bipolar referencing: For each electrode, computes contact_n - contact_{n+1}
%   e.g., RDa1-RDa2, RDa2-RDa3, etc. This reduces channels by one per electrode.
%
%   Input:
%       iSubj    - Subject index (1 to number of .mat files in data directory)
%       datadir  - Directory containing raw .mat files
%       outdir   - Directory for saving re-referenced output files
%       notch_hw - optional, half-width of the line-noise stop bands in Hz
%                  ([] or omitted: no filtering)
%
%   Output:
%       outfile  - the file written. It is saved under a temporary name and
%                  renamed once complete, so an interrupted run never leaves a
%                  half-written file under the final name.

    if nargin < 4, notch_hw = []; end

    % Create output directory if it doesn't exist
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Get list of subject files
    files = dir(fullfile(datadir, '*.mat'));
    if iSubj < 1 || iSubj > length(files)
        error('iSubj must be between 1 and %d', length(files));
    end

    fname = files(iSubj).name;
    fprintf('Processing subject %d/%d: %s\n', iSubj, length(files), fname);

    % Load raw data
    raw = load(fullfile(datadir, fname));
    data_eeg = raw.data.data_eeg;

    % Parse electrode labels to group contacts
    labels = data_eeg.label;
    nChannels = length(labels);

    % Extract electrode name and contact number from each label
    electrodes = struct();
    for iCh = 1:nChannels
        label = labels{iCh};
        % Match pattern: letters followed by numbers
        tokens = regexp(label, '^([A-Za-z]+)(\d+)$', 'tokens');
        if ~isempty(tokens)
            elecName = tokens{1}{1};
            contactNum = str2double(tokens{1}{2});

            if ~isfield(electrodes, elecName)
                electrodes.(elecName) = [];
            end
            electrodes.(elecName) = [electrodes.(elecName); contactNum, iCh];
        else
            warning('Could not parse label: %s', label);
        end
    end

    % Sort contacts within each electrode and prepare bipolar pairs
    elecNames = fieldnames(electrodes);
    nElectrodes = length(elecNames);

    bipolarPairs = [];  % [idx1, idx2] for each bipolar channel
    bipolarLabels = {};
    electrodeInfo = struct();

    for iElec = 1:nElectrodes
        elecName = elecNames{iElec};
        contacts = electrodes.(elecName);

        % Sort by contact number
        [~, sortIdx] = sort(contacts(:, 1));
        contacts = contacts(sortIdx, :);

        nContacts = size(contacts, 1);
        nBipolar = nContacts - 1;

        % Store electrode info
        electrodeInfo.(elecName) = nBipolar;

        % Create bipolar pairs (contact_n - contact_{n+1})
        for iPair = 1:nBipolar
            idx1 = contacts(iPair, 2);
            idx2 = contacts(iPair + 1, 2);
            bipolarPairs = [bipolarPairs; idx1, idx2];

            label1 = labels{idx1};
            label2 = labels{idx2};
            bipolarLabels{end + 1} = sprintf('%s-%s', label1, label2);
        end
    end

    nBipolarChannels = size(bipolarPairs, 1);
    fprintf('  Original channels: %d\n', nChannels);
    fprintf('  Bipolar channels: %d\n', nBipolarChannels);
    fprintf('  Electrodes: %d\n', nElectrodes);

    % Apply bipolar referencing to each trial
    nTrials = length(data_eeg.trial);
    bipolarTrials = cell(1, nTrials);

    for iTrial = 1:nTrials
        trialData = data_eeg.trial{iTrial};  % [nChannels x nSamples]
        bipolarData = trialData(bipolarPairs(:, 1), :) - trialData(bipolarPairs(:, 2), :);
        bipolarTrials{iTrial} = bipolarData;
        
        trialTime = data_eeg.time{iTrial}(1,:) - data_eeg.time{iTrial}(1,2001); % time relative to onset of the image
        timeTrials{iTrial} = trialTime;
    end

    % Create output structure (preserving data_eeg format)
    data_eeg_reref = struct();
    data_eeg_reref.trial = bipolarTrials;
%     data_eeg_reref.time = data_eeg.time;
    data_eeg_reref.time = timeTrials; % time relative to onset of the image
    data_eeg_reref.label = bipolarLabels(:);  % Column cell array
    data_eeg_reref.fsample = data_eeg.fsample;
    data_eeg_reref.electrodes = electrodeInfo;

    % Update cfg with referencing info
    if isfield(data_eeg, 'cfg')
        data_eeg_reref.cfg = data_eeg.cfg;
    else
        data_eeg_reref.cfg = struct();
    end
    data_eeg_reref.cfg.reref = 'bipolar';
    data_eeg_reref.cfg.reref_method = 'adjacent_subtraction';
    data_eeg_reref.cfg.reref_description = 'contact_n minus contact_{n+1} within each electrode';
    data_eeg_reref.cfg.original_channels = nChannels;
    data_eeg_reref.cfg.original_labels = labels;
    
    % Update cfg with time-correction info
    data_eeg_reref.cfg.time_description = 'time aligned to image onset';

    % Save output (only data_eeg, not data_et)
    data = struct();
    data.data_eeg = data_eeg_reref;

    if ~isempty(notch_hw)
        data = cleanTrialEEG(data, notch_hw);
    end

    % Demean every recording trial, per channel (2026-09-24). Channel offsets
    % differ between trials (in subNS167 71% of channels shift by more than a
    % within-trial SD, with steps at block boundaries); left in, they inflate the
    % pooled baseline SD the analyses z-score with and add noise to every
    % single-trial value. The mean of the whole 14 s trial is removed - not a
    % linear trend, which could bend slow image-locked responses.
    for iTrial = 1:numel(data.data_eeg.trial)
        x = data.data_eeg.trial{iTrial};
        data.data_eeg.trial{iTrial} = x - mean(x, 2);
    end
    data.data_eeg.cfg.trial_demean = 'per recording trial and channel: mean of the whole trial removed (preprocess_bipolar, 2026-09-24)';   % not cfg.demean: FieldTrip's own field of that name comes with the raw cfg

    outfile = fullfile(outdir, ['reref_' fname]);
    partial = [outfile(1:end-4) '.partial.mat'];
    fprintf('  Saving to: %s\n', outfile);
    save(partial, 'data', '-v7.3');
    movefile(partial, outfile, 'f');

    fprintf('  Done.\n\n');
end
