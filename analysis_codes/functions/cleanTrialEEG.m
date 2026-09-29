function data = cleanTrialEEG(data, half_width)
% CLEANTRIALEEG  Line-noise removal on re-referenced trial-level EEG data.
%
%   data = cleanTrialEEG(data)
%   data = cleanTrialEEG(data, half_width)
%
% 4th-order Butterworth band-stop, zero-phase (filtfilt on second-order
% sections), at 60 Hz and its harmonics 120, 180, 240 Hz, +-half_width Hz
% (default 5), per trial. Records itself in data.data_eeg.cfg.bandstopFilter.
%
% Why +-5 Hz (2026-09-24, notes/trial_cleaning.md): the recordings carry
% transient line-noise bursts of ~200 ms, whose energy spreads +-5 Hz around
% 60 Hz and each harmonic; the earlier +-2 Hz notch left them, and a 58/62 Hz
% leak in the noisiest channels. The wider stop bands cost BHA (70-150 Hz) the
% bins around 120 Hz - see compute_bha_stats.
%
% Input/Output:
%   data       - struct with data.data_eeg (.trial, .fsample, .cfg), bipolar
%                trials straight from preprocess_bipolar
%   half_width - half-width of each stop band in Hz (default 5)
%
% A file that already carries a band-stop is left alone if it is the same one,
% and refused if it is a different one (e.g. the +-2 Hz of files made before
% 2026-09-24): regenerate it from raw (run_preprocess step 2) instead of
% filtering twice.

if nargin < 2 || isempty(half_width), half_width = 5; end
data_eeg  = data.data_eeg;
fs        = data_eeg.fsample;
line_freq = 60;
harm      = 4;
order     = 4;

if isfield(data_eeg.cfg, 'bandstopFilter')
    old_hw = getopt_field(data_eeg.cfg.bandstopFilter, 'half_width', 2);   % files before 2026-09-24: +-2 Hz
    if old_hw == half_width
        fprintf('Band-stop (+-%g Hz) already applied, skipping\n', half_width);
        return;
    end
    error(['This file already carries a +-%g Hz band-stop; regenerate it from raw ' ...
           '(run_preprocess step 2) before applying +-%g Hz.'], old_hw, half_width);
end

% FieldTrip's external/signal filtfilt shadows MATLAB's and cannot filter
% second-order sections: step it off the path for this call, restore on exit.
old_path = path;
restore  = onCleanup(@() path(old_path));
while ~contains(which('filtfilt'), fullfile('toolbox', 'signal'))
    shadow = fileparts(which('filtfilt'));
    if isempty(shadow), error('filtfilt (Signal Processing Toolbox) not found.'); end
    rmpath(shadow);
end

SOS = cell(1, harm); G = zeros(1, harm);
for h = 1:harm
    [z, p, k] = butter(order, [h*line_freq - half_width, h*line_freq + half_width] / (fs/2), 'stop');
    [SOS{h}, G(h)] = zp2sos(z, p, k);
end

nTrials = numel(data_eeg.trial);
for iTrial = 1:nTrials
    d = double(data_eeg.trial{iTrial})';          % samples x channels
    for h = 1:harm
        d = filtfilt(SOS{h}, G(h), d);
    end
    data_eeg.trial{iTrial} = d';
    if mod(iTrial, 10) == 0
        fprintf('  Filtered trial %d/%d\n', iTrial, nTrials);
    end
end

data_eeg.cfg.bandstopFilter = struct('line_freq', line_freq, 'harmonics', harm, 'order', order, ...
    'half_width', half_width, 'bands', (1:harm)'*line_freq + [-half_width half_width], ...
    'date', char(datetime('today', 'Format', 'yyyy-MM-dd')));
fprintf('Applied band-stop at %d Hz and %d harmonics, +-%g Hz, to %d trials\n', ...
        line_freq, harm - 1, half_width, nTrials);

data.data_eeg = data_eeg;
end
