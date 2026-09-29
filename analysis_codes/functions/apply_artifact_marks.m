function [data, mask] = apply_artifact_marks(data, marks, opts)
% APPLY_ARTIFACT_MARKS  NaN out marked samples of a plain FieldTrip raw structure.
%
%   data         = apply_artifact_marks(data, marks)
%   [data, mask] = apply_artifact_marks(data, marks, opts)
%
% Sample-level counterpart of cleanPerisaccadicEEG for the continuous trials
% the marks were computed on (mark_artifacts + annotate_marks). Instead of
% flagging whole channel-epochs, every marked sample is set to NaN in place;
% no fields are added to data, which stays a standard FieldTrip raw structure.
%
% Excluded by default, as .ok of cleanPerisaccadicEEG:
%   amp, diff, hf       artifact steps
%   ied_other           IEDs not locked to saccades
%   bad_channel         the whole channel
% opts.strict = true also excludes ied_locked (as .ok_strict).
%
% Inputs:
%   data  - FieldTrip raw (.trial, .time, .label, .fsample), the same trials
%           mark_artifacts ran on (same labels, trial count and length)
%   marks - from mark_artifacts (+ annotate_marks)
%   opts  - optional:
%           .strict       (false) also NaN saccade-locked IEDs
%           .steps        ({'amp','diff','hf'}) artifact steps to apply
%           .bad_channels ('nan') 'nan' whole channel, 'remove' drop it from
%                         data, 'keep' leave it untouched
%           .pad          (0 s) extra margin added around every marked stretch
%
% Output:
%   data  - with NaN in the marked samples
%   mask  - {1 x nTr} logical nCh x nS, true where data was set to NaN
%           (channels as in the output data)

if nargin < 3, opts = struct(); end
strict   = getopt_field(opts, 'strict', false);
steps    = getopt_field(opts, 'steps', {'amp', 'diff', 'hf'});
bad_mode = getopt_field(opts, 'bad_channels', 'nan');
pad      = getopt_field(opts, 'pad', 0);

% --- Consistency with the marks ---
fs  = data.fsample;
nTr = numel(data.trial);
nS  = marks.n_samples;
if abs(fs - marks.fs) > 0.5
    error('Sampling rate of the data (%g) and of the marks (%g) differ.', fs, marks.fs);
end
if ~isequal(data.label(:), marks.label(:))
    error('Channel labels of the data and of the marks differ.');
end
if nTr ~= marks.n_trials
    error('The data has %d trials, the marks %d.', nTr, marks.n_trials);
end
if any(cellfun(@(x) size(x, 2), data.trial) ~= nS)
    error('Trial length of the data differs from the marks (%d samples).', nS);
end
t0 = cellfun(@(tt) tt(1), data.time);
if any(abs(t0(:) - marks.t0(:)) > 0.5/fs)
    error('Trial start times of the data and of the marks differ.');
end

% --- One sparse mask of everything excluded ---
M = sparse(numel(marks.label), nTr*nS);
for s = 1:numel(steps)
    M = M | marks.(steps{s});
end
if isfield(marks, 'ied_other')
    M = M | marks.ied_other;
    if strict, M = M | marks.ied_locked; end
else
    warning('apply_artifact_marks:noAnnotation', ...
        'marks have no annotate_marks fields; all IEDs (marks.ied) are excluded.');
    M = M | marks.ied;
end

bad = false(1, numel(marks.label));
if isfield(marks, 'bad_channel'), bad = marks.bad_channel; end

% --- NaN the marked samples, trial by trial ---
nPad = round(pad * fs);
mask = cell(1, nTr);
for tr = 1:nTr
    m = full(M(:, (tr-1)*nS + (1:nS)));
    if nPad > 0
        m = conv2(double(m), ones(1, 2*nPad + 1), 'same') > 0;
    end
    if strcmp(bad_mode, 'nan'), m(bad, :) = true; end
    data.trial{tr}(m) = NaN;
    mask{tr} = m;
end

if strcmp(bad_mode, 'remove')
    keep = ~bad;
    data.label = data.label(keep);
    for tr = 1:nTr
        data.trial{tr} = data.trial{tr}(keep, :);
        mask{tr}       = mask{tr}(keep, :);
    end
    if isfield(data, 'elec'), warning('apply_artifact_marks:elec', ...
            'data.elec is not pruned to the kept channels.'); end
end

nCh = numel(marks.label);
fprintf('apply_artifact_marks: %.1f%% of samples set to NaN; %d/%d bad channels (%s)%s\n', ...
        100*mean(cellfun(@(x) mean(x, 'all'), mask)), sum(bad), nCh, bad_mode, ...
        repmat(', strict', 1, strict));
end
