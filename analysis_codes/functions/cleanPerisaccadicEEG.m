function data_perisac = cleanPerisaccadicEEG(data_perisac, marks)
% CLEANPERISACCADICEEG  Artifact flags for every channel-epoch of an epoch set.
%
%   data_perisac = cleanPerisaccadicEEG(data_perisac, marks)
%
% `marks` come from get_artifact_marks (mark_artifacts + annotate_marks) on the
% reref trials the epochs were cut from. A channel-epoch is flagged by a step if
% ANY of its samples - the whole epoch, padding included - carries that step's
% mark, so every analysis reading any part of an epoch sees the same trial set,
% and every epoch set of a subject is judged against the same marks. Bad
% channels come from the marks too (annotate_marks), so they are the same in
% every epoch set.
%
% Output, data_perisac.noiseInfo:
%   per channel-epoch (nTrials x nChannels logical, complete epochs only):
%     .amp .diff .hf          artifact steps
%     .ied_other              an IED not locked to saccades
%     .ied_locked             an IED in its channel's saccade-locked window
%     .ied                    either
%     .ok                     complete, no artifact, no other-IED, channel not bad
%     .ok_strict              .ok and no locked IED either
%   per channel (1 x nChannels):
%     .bad_channel, .bad_frac            (annotate_marks, all epoch sets)
%     .ied_locked_channel                the channel has a saccade-locked window
%     .saccadic_spike_channel, .saccadic_spike_z
%     .line_noise_channel
%     .ied_rate_per_min, .rejected_frac  (this set: complete epochs not .ok)
% data_perisac.params.artifactDetection records method and parameters
% (version 2) and is the done-flag; a file without .version (the per-epoch
% version 1) is cleaned again. See notes/trial_cleaning.md.

if isfield(data_perisac, 'artifactDetection') && ...
        getopt_field(data_perisac.artifactDetection, 'version', 1) >= 2
    fprintf('Artifact detection (v2) already applied, skipping\n');
    return;
end
if nargin < 2 || isempty(marks) || ~isfield(marks, 'bad_channel')
    error(['cleanPerisaccadicEEG needs the marks of get_artifact_marks ' ...
           '(mark_artifacts + annotate_marks) for the reref trials these epochs were cut from.']);
end

fs = data_perisac.fsample;
if abs(fs - marks.fs) > 0.5          % params.fs is rounded, fsample need not be
    error('Sampling rate of the epochs (%g) and of the marks (%g) differ.', fs, marks.fs);
end
if ~isequal(data_perisac.label(:), marks.label(:))
    error('Channel labels of the epochs and of the marks differ.');
end

% --- Mark complete epochs (no NaN padding) ---
nTrials = numel(data_perisac.trial);
if ~ismember('is_complete', data_perisac.eventInfo.Properties.VariableNames)
    is_complete = false(nTrials, 1);
    for iTrial = 1:nTrials
        is_complete(iTrial) = ~any(isnan(data_perisac.trial{iTrial}), 'all');
    end
    data_perisac.eventInfo.is_complete = is_complete;
    fprintf('Complete epochs: %d/%d\n', sum(is_complete), nTrials);
end
is_complete = data_perisac.eventInfo.is_complete;

% --- Each complete epoch's samples in the concatenated trials ---
nS   = marks.n_samples;
nT   = numel(data_perisac.time);
nCh  = numel(marks.label);
comp = find(is_complete);
ev   = data_perisac.eventInfo;

steps = {'amp', 'diff', 'hf', 'ied_other', 'ied_locked'};
for s = 1:numel(steps), F.(steps{s}) = false(nTrials, nCh); end
for e = comp'
    tr = ev.trial_number(e);
    i0 = round((ev.closest_sample_time(e) + data_perisac.time(1) - marks.t0(tr)) * fs) + 1;
    if i0 < 1 || i0 + nT - 1 > nS
        error('Epoch %d is complete but leaves its recording trial %d.', e, tr);
    end
    g = (tr - 1)*nS + i0 + (0:nT-1);
    for s = 1:numel(steps)
        F.(steps{s})(e, :) = any(marks.(steps{s})(:, g), 2)';
    end
end

artifact = F.amp | F.diff | F.hf;
ok = repmat(is_complete, 1, nCh) & ~artifact & ~F.ied_other & ~marks.bad_channel;

for s = 1:numel(steps), data_perisac.noiseInfo.(steps{s}) = F.(steps{s}); end
data_perisac.noiseInfo.ied       = F.ied_other | F.ied_locked;
data_perisac.noiseInfo.ok        = ok;
data_perisac.noiseInfo.ok_strict = ok & ~F.ied_locked;
data_perisac.noiseInfo.bad_channel            = marks.bad_channel;
data_perisac.noiseInfo.bad_frac               = marks.bad_frac;
data_perisac.noiseInfo.ied_locked_channel     = ~cellfun(@isempty, {marks.locking.window});
data_perisac.noiseInfo.saccadic_spike_channel = marks.saccadic_spike_channel;
data_perisac.noiseInfo.saccadic_spike_z       = marks.saccadic_spike_z;
data_perisac.noiseInfo.line_noise_channel     = marks.line_noise_channel;
data_perisac.noiseInfo.ied_rate_per_min       = marks.ied_rate_per_min;
data_perisac.noiseInfo.rejected_frac          = mean(~ok(comp, :), 1);

P = marks.params; A = marks.annot_params;
data_perisac.artifactDetection = struct( ...
    'version', 2, ...
    'date', P.date, ...
    'rule', ['channel-epoch flagged if any of its samples, padding included, is marked; marks computed once per ' ...
             'subject on the continuous reref trials (get_artifact_marks) and shared by all epoch sets'], ...
    'fence', sprintf(['per channel and step: mark samples whose statistic exceeds exp(Q3 + %g*IQR) of the log ' ...
                      'per-%g-s maxima (Tukey far-out fence), trial-demeaned signal'], P.fence_k, P.seg_len), ...
    'step_amp',   sprintf('|%.2f s moving mean|; whole window marked', P.amp_win), ...
    'step_diff',  sprintf('|first difference - median|; +-%.2f s marked', P.margin), ...
    'step_hf',    sprintf('|FieldTrip FIR high-pass %g Hz two-pass - median|; +-%.2f s marked', P.hf_cutoff, P.margin), ...
    'step_ied',   sprintf(['detect_ieds (after Michal''s findSpikeTimes_reviewed), both polarities, shape floor %g SD; ' ...
                           'threshold on depth+height set where a phase-randomised null gives FDR <= %g (this subject: %.2f SD); ' ...
                           '+-%.2f s marked'], ...
                          P.ied_z, P.ied_fdr, marks.ied_threshold, P.ied_margin), ...
    'ied_locking', sprintf(['per channel, IED-saccade onset cross-correlogram (%g s bins) vs its |lag| >= %g s baseline, ' ...
                            'Poisson + BH-FDR %g over channels x bins in [%g %g] s; IEDs in the locked window are ' ...
                            'ied_locked (excluded only by ok_strict), all others ied_other (excluded by ok)'], ...
                           A.bin, A.baseline_min, A.fdr, A.search), ...
    'saccadic_spike', sprintf(['per channel, max |saccade-onset-locked mean of the %g Hz high-passed signal| in [%g %g] s; ' ...
                               'flagged if above its SE-based Bonferroni z = %.2f (alpha %g) and its size (robust SD) ' ...
                               'beyond Q3 + %g IQR across channels; label only'], ...
                              A.sp_hp, A.sp_win, marks.saccadic_spike_crit, A.sp_alpha, A.sp_fence_k), ...
    'bad_channel', sprintf(['more than %g of the %.1f s segments tiling the trials hold an excluding mark ' ...
                            '(artifact or other IED); same for all epoch sets'], A.bad_chan_frac, A.bad_seg_len), ...
    'marks_params', P, 'annot_params', A);

nOk = sum(ok(:)); nTotal = numel(comp) * nCh;
pc = @(M) 100*mean(M(comp, :), 'all');
fprintf(['Artifact detection v2: amp %.1f%%  diff %.1f%%  hf %.1f%%  IED other %.1f%%  locked %.1f%%; ' ...
         '%d/%d bad channels; ok %d/%d (%.1f%%), ok_strict %.1f%%\n'], ...
        pc(F.amp), pc(F.diff), pc(F.hf), pc(F.ied_other), pc(F.ied_locked), ...
        sum(marks.bad_channel), nCh, nOk, nTotal, 100*nOk/nTotal, ...
        100*sum(data_perisac.noiseInfo.ok_strict(:))/nTotal);
end
