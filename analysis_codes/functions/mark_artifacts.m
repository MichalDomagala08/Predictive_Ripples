function marks = mark_artifacts(data_eeg, opts)
% MARK_ARTIFACTS  Sample-level artifact and IED marks on the continuous trials.
%
%   marks = mark_artifacts(data_eeg)
%   marks = mark_artifacts(data_eeg, opts)
%
% Runs once per subject on the reref trials, before any epoching and without
% any behavioural information, so every epoch set cut from these trials is
% judged against the same marks. annotate_marks then adds what needs saccade
% times (IED locking, saccadic-spike channels) and the bad-channel verdict, and
% cleanPerisaccadicEEG rejects a channel-epoch if any of its samples, padding
% included, is marked. get_artifact_marks does all of it once and caches it.
%
% ONE RULE for the three artifact steps. For each channel, the step's statistic
% s(t) >= 0 is computed on the trial-demeaned signal; its maximum is taken in
% every 1 s segment of the trials; a sample is marked where s(t) exceeds the
% Tukey "far out" fence of those maxima on a log scale:
%       fence = exp(Q3 + k*(Q3 - Q1)) of log(per-second maxima),  k = 3.
% So a mark means "more extreme than this channel's ordinary seconds get", with
% the same k for every step, whatever the channel's own level of transients.
%   amp   |0.25 s moving mean|; the whole window is marked
%   diff  |first difference - median|
%   hf    |FieldTrip FIR high-pass 225 Hz, two-pass - median| (Michal's filter)
%   diff and hf marks are widened by +-margin (0.2 s, as Michal).
% (A fourth step, 55-65 Hz burst envelope, was dropped on 2026-09-24: the +-5 Hz
% notch of cleanTrialEEG now removes the transient line-noise bursts it caught,
% and after it the step only re-caught broadband sharp events.)
%
% IEDs. detect_ieds (Michal's detector, both polarities) scores every candidate
% by trough depth + peak height in SD of the trial-demeaned channel. The same
% detector runs on a phase-randomised copy of every trial (same spectrum, each
% trial at its own robust SD, no transients); the threshold is the lowest score
% at which null/observed detections, pooled over channels, stay <= ied_fdr (1 %).
% Detections closer than ied_match_win are merged (one discharge, both passes).
% +-ied_margin (0.1 s) is marked around each.
%
% Line noise (a channel flag, not a mark): the 57-63 Hz power left inside the
% notch, over 40-45 and 75-80 Hz (log10); a channel is flagged when this is
% beyond the same Tukey fence across the subject's channels - its line noise is
% strong enough to survive the +-5 Hz band-stop of cleanTrialEEG.
%
% Inputs:
%   data_eeg - data.data_eeg of a reref file (.trial, .time, .fsample, .label)
%   opts     - optional: .fence_k (3), .seg_len (1 s), .amp_win (0.25 s),
%              .hf_cutoff (225), .margin (0.2 s),
%              .ied_fdr (0.01), .ied_z (4), .ied_max_width (50 samples),
%              .ied_match_win (0.1 s), .ied_margin (0.1 s), .ied_polarity
%              ('both'), .seed (1)
%
% Output marks:
%   .amp .diff .hf .ied           nCh x (nTr*nS) sparse logical, trials concatenated
%   .fence                        nCh x 3 (amp, diff, hf), in the statistic's units
%   .ied_times, .ied_scores       {1 x nCh}, global sample index / score of each IED
%   .ied_threshold, .ied_fdr_curve (grid, observed and null rate per channel-minute)
%   .line_excess, .line_noise_channel   1 x nCh
%   .n_samples, .n_trials, .fs, .t0, .label, .ied_rate_per_min, .params

if nargin < 2, opts = struct(); end
P.fence_k       = getopt_field(opts, 'fence_k', 3);
P.seg_len       = getopt_field(opts, 'seg_len', 1);
P.amp_win       = getopt_field(opts, 'amp_win', 0.25);
P.hf_cutoff     = getopt_field(opts, 'hf_cutoff', 225);
P.margin        = getopt_field(opts, 'margin', 0.2);
P.ied_fdr       = getopt_field(opts, 'ied_fdr', 0.01);   % 1 % since 2026-09-24 (5 % let oscillation bursts through, notes/trial_cleaning.md)
P.ied_z         = getopt_field(opts, 'ied_z', 4);
P.ied_max_width = getopt_field(opts, 'ied_max_width', 50);
P.ied_match_win = getopt_field(opts, 'ied_match_win', 0.1);
P.ied_margin    = getopt_field(opts, 'ied_margin', 0.1);
P.ied_polarity  = getopt_field(opts, 'ied_polarity', 'both');
P.seed          = getopt_field(opts, 'seed', 1);

fs  = data_eeg.fsample;
nTr = numel(data_eeg.trial);
[nCh, nS] = size(data_eeg.trial{1});
if any(cellfun(@(x) size(x, 2), data_eeg.trial) ~= nS)
    error('mark_artifacts: trials differ in length.');
end
N   = nTr * nS;
tix = @(t) (t-1)*nS + (1:nS);
seg = round(P.seg_len * fs);
nSeg = floor(nS / seg);

X = zeros(nCh, N);
for t = 1:nTr
    x = double(data_eeg.trial{t});
    X(:, tix(t)) = x - mean(x, 2);             % trial-demeaned
end

% statistic -> fence -> marks (dilated within each trial by `dil` samples)
    function [m, fence, M] = fence_marks(S, dil)
        M = zeros(nCh, nSeg*nTr);
        for tt = 1:nTr
            St = S(:, tix(tt));
            M(:, (tt-1)*nSeg + (1:nSeg)) = squeeze(max(reshape(St(:, 1:nSeg*seg), nCh, seg, nSeg), [], 2));
        end
        L = log(max(M, realmin));
        q = row_quartiles(L);
        fence = exp(q(:, 2) + P.fence_k * (q(:, 2) - q(:, 1)));
        m = false(nCh, N);
        for tt = 1:nTr
            o = S(:, tix(tt)) > fence;
            if dil > 0, o = conv2(double(o), ones(1, 2*dil + 1), 'same') > 0; end
            m(:, tix(tt)) = o;
        end
    end

fence = zeros(nCh, 3); seg_max = cell(1, 3);
w     = round(P.amp_win * fs);
dil   = round(P.margin * fs);

%% amp
S = zeros(nCh, N);
for t = 1:nTr, S(:, tix(t)) = abs(movmean(X(:, tix(t)), w, 2, 'Endpoints', 'shrink')); end
[amp, fence(:, 1), seg_max{1}] = fence_marks(S, floor(w/2));

%% diff
S = zeros(nCh, N);
for t = 1:nTr
    idx = tix(t);
    S(:, idx(1:end-1)) = diff(X(:, idx), 1, 2);
end
S = abs(S - median(S, 2));
[dif, fence(:, 2), seg_max{2}] = fence_marks(S, dil);

%% hf
for t = 1:nTr
    S(:, tix(t)) = ft_preproc_highpassfilter(X(:, tix(t)), fs, P.hf_cutoff, [], 'fir', 'twopass');
end
S = abs(S - median(S, 2));
[hf, fence(:, 3), seg_max{3}] = fence_marks(S, dil);

clear S

%% line-noise channel flag
nfw = round(2*fs);
[pxx, f] = pwelch(X(:, 1:min(N, 200*nfw))', nfw, nfw/2, nfw, fs);
band = @(lo, hi) mean(log10(pxx(f >= lo & f <= hi, :)), 1);
line_excess = band(57, 63) - mean([band(40, 45); band(75, 80)], 1);
qc = row_quartiles(line_excess);
line_noise_channel = line_excess > qc(2) + P.fence_k * (qc(2) - qc(1));

%% IEDs: candidates, phase-randomised null, FDR threshold
rng(P.seed);
iopt = struct('z', P.ied_z, 'amp_scale', 0, 'max_width', P.ied_max_width, ...
              'match_win', P.ied_match_win, 'polarity', P.ied_polarity);
hp = floor(P.ied_match_win * fs);
tc = cell(1, nCh); sc = cell(1, nCh); sn = cell(1, nCh);
for ch = 1:nCh
    x = X(ch, :);
    [t, ~, s] = detect_ieds((x - mean(x)) / std(x), fs, iopt);
    [t, s] = merge_events(t, s, hp);
    tc{ch} = t; sc{ch} = s;
    Xt = reshape(x, nS, nTr)';
    xs = real(ifft(abs(fft(Xt, [], 2)) .* exp(1i*angle(fft(randn(size(Xt)), [], 2))), [], 2));
    rs = @(A) median(abs(A - median(A, 2)), 2);
    xs = reshape((xs ./ rs(xs) .* rs(Xt))', 1, []);
    [t, ~, s] = detect_ieds((xs - mean(xs)) / std(xs), fs, iopt);
    [~, s] = merge_events(t, s, hp);
    sn{ch} = s;
end
grid = P.ied_z:0.25:40;
minutes = N / fs / 60;
so = [sc{:}]; snn = [sn{:}];
ro = arrayfun(@(g) sum(so >= g), grid) / (nCh * minutes);
rn = arrayfun(@(g) sum(snn >= g), grid) / (nCh * minutes);
fdr = rn ./ max(ro, eps);
iT = find(fliplr(cummax(fliplr(fdr))) <= P.ied_fdr, 1);
if isempty(iT), ied_thr = grid(end); else, ied_thr = grid(iT); end

ied = false(nCh, N);
ied_times = cell(1, nCh); ied_scores = cell(1, nCh);
mI = round(P.ied_margin * fs);
for ch = 1:nCh
    k = sc{ch} >= ied_thr;
    ied_times{ch} = tc{ch}(k); ied_scores{ch} = sc{ch}(k);
    ti = ied_times{ch};
    tr = ceil(ti / nS);
    lo_i = max(ti - mI, (tr-1)*nS + 1);
    hi_i = min(ti + mI, tr*nS);
    for j = 1:numel(ti), ied(ch, lo_i(j):hi_i(j)) = true; end
end
clear X

%% pack
marks.amp   = sparse(amp);
marks.diff  = sparse(dif);
marks.hf    = sparse(hf);
marks.ied   = sparse(ied);
marks.fence = fence;
marks.seg_max = cellfun(@single, seg_max, 'UniformOutput', false);   % per-second maxima the fences came from
marks.ied_times  = ied_times;
marks.ied_scores = ied_scores;
marks.ied_threshold = ied_thr;
marks.ied_fdr_curve = struct('grid', grid, 'observed', ro, 'null', rn, 'fdr', fdr);
marks.ied_rate_per_min = cellfun(@numel, ied_times) / minutes;
marks.line_excess = line_excess;
marks.line_noise_channel = line_noise_channel;
marks.n_samples = nS;
marks.n_trials  = nTr;
marks.fs        = fs;
marks.t0        = cellfun(@(tt) tt(1), data_eeg.time);
marks.label     = data_eeg.label(:);
P.date = char(datetime('today', 'Format', 'yyyy-MM-dd'));
marks.params = P;

pct = @(m) 100*nnz(m)/numel(m);
fprintf(['  marked samples: amp %.2f%%  diff %.2f%%  hf %.2f%%  ied %.2f%% ' ...
         '(IED threshold %.2f SD, median %.2f IEDs/min per channel); %d line-noise channels\n'], ...
        pct(amp), pct(dif), pct(hf), pct(ied), ied_thr, ...
        median(marks.ied_rate_per_min), sum(line_noise_channel));
end

function [t, s] = merge_events(t, s, gap)
% One event per cluster of detections closer than `gap` samples, at the
% highest-scoring one (a discharge found by both polarity passes counts once).
    if isempty(t), return; end
    [t, o] = sort(t); s = s(o);
    c = cumsum([true, diff(t) > gap]);
    keep = false(size(t));
    for k = 1:c(end)
        idx = find(c == k);
        [~, j] = max(s(idx));
        keep(idx(j)) = true;
    end
    t = t(keep); s = s(keep);
end

function q = row_quartiles(X)
% First and third quartile of each row, as quantile(X, [0.25 0.75], 2) computes
% them (linear interpolation at (k - 0.5)/n), without the Statistics Toolbox.
    n  = size(X, 2);
    Xs = sort(X, 2);
    q  = zeros(size(X, 1), 2);
    p  = [0.25 0.75];
    for k = 1:2
        h = n*p(k) + 0.5;
        lo = max(1, floor(h)); hi = min(n, lo + 1);
        q(:, k) = Xs(:, lo) + (h - lo) * (Xs(:, hi) - Xs(:, lo));
    end
end

function [t_ied, n_events, score] = detect_ieds(z, fs, opts)
% DETECT_IEDS  Interictal epileptiform discharges in one channel.
%
%   [t_ied, n_events] = detect_ieds(z, fs, opts)
%
% The detector of Michal's findSpikeTimes_reviewed (od_michala/cleaning_trials),
% with the parameters he runs it with (runScripts_review_check.m), rewritten
% without the buffering and run in both polarities: a bipolar channel's sign is
% arbitrary, and his pipeline computes the positive-first pass but never uses it.
%
% A candidate is a trough of depth >= opts.z and prominence >= opts.z whose width
% at half prominence is <= opts.max_width samples. It is kept if a peak of
% prominence >= opts.z (height >= 0) lies within +-opts.match_win s of it, with
% trough depth + peak height >= opts.amp_scale * opts.z. The positive pass is the
% same on -z.
%
% Inputs:
%   z    - 1 x N, the channel z-scored over the whole recording
%   fs   - sampling rate, Hz
%   opts - .z (4), .amp_scale (1; Michal's paper text says 12 SD, i.e. 3),
%          .max_width (50 samples), .match_win (0.1 s), .polarity ('both' |
%          'neg' | 'pos')
%
% Outputs:
%   t_ied    - sample indices of the detected troughs (both passes, sorted)
%   n_events - detections after merging those closer than match_win (one
%              discharge found by both passes counts once)
%   score    - per detection, trough depth + height of the highest matching
%              peak, in SD (what amp_scale*z is compared with)

if nargin < 3, opts = struct(); end
thr       = getopt_field(opts, 'z', 4);
amp_scale = getopt_field(opts, 'amp_scale', 1);
max_width = getopt_field(opts, 'max_width', 50);
match_win = getopt_field(opts, 'match_win', 0.1);
polarity  = getopt_field(opts, 'polarity', 'both');

z  = z(:)';
hp = floor(match_win * fs);
switch polarity
    case 'neg',  passes = {z};
    case 'pos',  passes = {-z};
    case 'both', passes = {z, -z};
    otherwise,   error('polarity must be ''neg'', ''pos'' or ''both''.');
end

w_state = warning('off', 'signal:findpeaks:largeMinPeakHeight');   % no trough that deep: fine
cleanup = onCleanup(@() warning(w_state));
t_ied = []; score = [];
for k = 1:numel(passes)
    y = passes{k};
    [pH, pI] = findpeaks(y,  'MinPeakProminence', thr, 'MinPeakHeight', 0);
    [nH, nI] = findpeaks(-y, 'MinPeakProminence', thr, 'MinPeakHeight', thr, ...
                         'MaxPeakWidth', max_width);
    if isempty(pI) || isempty(nI), continue; end
    P = -inf(size(y));  P(pI) = pH;
    best = movmax(P, [hp hp]);                 % highest matching peak near each sample
    keep = best(nI) + nH >= amp_scale * thr;   % -inf when no peak is in range
    t_ied = [t_ied, nI(keep)]; %#ok<AGROW>
    sc = best(nI) + nH;
    score = [score, sc(keep)]; %#ok<AGROW>
end
[t_ied, o] = sort(t_ied);
score = score(o);
n_events = sum(diff([-inf, t_ied]) > hp);
end


