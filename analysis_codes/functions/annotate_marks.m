function marks = annotate_marks(marks, data_eeg, content, opts)
% ANNOTATE_MARKS  What the artifact marks mean relative to saccades, and which
% channels are bad. The only step of the cleaning that sees saccade times; it
% labels, it does not change what mark_artifacts found.
%
%   marks = annotate_marks(marks, data_eeg, content)
%   marks = annotate_marks(marks, data_eeg, content, opts)
%
% 1. IEDs: saccade-locked vs other, per event.
%    Per channel, the cross-correlogram of IED times against saccade onsets
%    (same trial, lags -1..1 s, 20 ms bins) is compared with its baseline (mean
%    count per bin at |lag| >= 0.5 s): each bin in the search range
%    [-0.3, 0.5] s gets a Poisson p-value, and BH-FDR is applied over all bins of
%    all channels. A channel's locked window is the run of significant bins
%    around its most significant one. An IED with a saccade onset such that its
%    lag falls in that window is `locked`; every other IED is `other`.
%    Per channel: .window (s), .n_locked, .n_expected (in-window events the
%    baseline predicts, i.e. locked-tagged events that are probably ordinary IEDs).
% 2. Saccadic spike channels. Per channel, the saccade-onset-locked average of the
%    30 Hz high-passed signal, -20..+30 ms (the spike potential is fast; the
%    corneoretinal step every saccade carries is slow and drops out). Flagged
%    when its peak is significant (over its standard error, Bonferroni over time
%    points and channels at alpha 0.01) AND its size in robust SD of the
%    channel is beyond the Tukey far-out fence across the subject's channels
%    (Q3 + 3 IQR of the sizes). A label only.
% 3. Bad channels, one verdict for every epoch set: the trials are tiled into
%    bad_seg_len (2.6 s, the saccade union epoch) segments; a channel is bad when
%    more than bad_chan_frac (0.5) of them hold a mark that the default trial set
%    excludes (amp, diff, hf, other IED). Locked IEDs do not count.
%
% Inputs:
%   marks    - from mark_artifacts
%   data_eeg - the same reref trials (for the saccadic-spike average)
%   content  - saccade table (saccade_onset_time, trial_number)
%   opts     - .lag_max (1 s), .bin (0.02 s), .search ([-0.3 0.5] s),
%              .baseline_min (0.5 s), .fdr (0.05), .sp_win ([-0.02 0.03] s),
%              .sp_alpha (0.01), .bad_seg_len (2.6 s), .bad_chan_frac (0.5)
%
% Adds to marks: .ied_locked, .ied_other (sparse, like .ied), .ied_locked_times,
% .ied_other_times, .locking (per-channel struct array), .saccadic_spike_z,
% .saccadic_spike_channel, .bad_frac, .bad_channel, .annot_params.

if nargin < 4, opts = struct(); end
A.lag_max       = getopt_field(opts, 'lag_max', 1);
A.bin           = getopt_field(opts, 'bin', 0.02);
A.search        = getopt_field(opts, 'search', [-0.3 0.5]);
A.baseline_min  = getopt_field(opts, 'baseline_min', 0.5);
A.fdr           = getopt_field(opts, 'fdr', 0.05);
A.sp_win        = getopt_field(opts, 'sp_win', [-0.02 0.03]);
A.sp_alpha      = getopt_field(opts, 'sp_alpha', 0.01);
A.sp_hp         = getopt_field(opts, 'sp_hp', 30);
A.sp_fence_k    = getopt_field(opts, 'sp_fence_k', 3);
A.bad_seg_len   = getopt_field(opts, 'bad_seg_len', 2.6);
A.bad_chan_frac = getopt_field(opts, 'bad_chan_frac', 0.5);

fs = marks.fs; nS = marks.n_samples; nTr = marks.n_trials; N = nS*nTr;
nCh = numel(marks.label);

% saccade onsets as global samples
okc = ~isnan(content.saccade_onset_time) & content.trial_number >= 1 & content.trial_number <= nTr;
tr  = content.trial_number(okc);
g_on = (tr - 1)*nS + round((content.saccade_onset_time(okc) - marks.t0(tr)')*fs) + 1;
in_trial = g_on >= (tr-1)*nS + 1 & g_on <= tr*nS;
g_on = sort(g_on(in_trial)); tr_on = ceil(g_on / nS);

%% 1. IED locking
edges = -A.lag_max:A.bin:A.lag_max;
ctr = edges(1:end-1) + A.bin/2;
srch = ctr >= A.search(1) & ctr <= A.search(2);
base = abs(ctr) >= A.baseline_min;
C = zeros(nCh, numel(ctr)); lam = zeros(nCh, 1);
pv = nan(nCh, numel(ctr));
lagmax = round(A.lag_max * fs);
for ch = 1:nCh
    ti = marks.ied_times{ch};
    lags = [];
    for j = 1:numel(ti)
        k = find(tr_on == ceil(ti(j)/nS) & abs(ti(j) - g_on) <= lagmax);
        lags = [lags; (ti(j) - g_on(k)) / fs]; %#ok<AGROW>
    end
    C(ch, :) = histcounts(lags, edges);
    lam(ch) = mean(C(ch, base));
    if lam(ch) > 0
        k = C(ch, srch);
        p = ones(size(k));
        p(k >= 1) = gammainc(lam(ch), k(k >= 1));        % P(X >= k), Poisson(lam)
        pv(ch, srch) = p;
    end
end
qv = apply_bh_fdr(pv);

locked_mask = false(nCh, N); other_mask = false(nCh, N);
mI = round(marks.params.ied_margin * fs);
locking = struct('window', cell(1, nCh), 'n_locked', 0, 'n_other', 0, 'n_expected', 0, ...
                 'counts', [], 'baseline', 0);
lt = cell(1, nCh); ot = cell(1, nCh);
for ch = 1:nCh
    sig = qv(ch, :) < A.fdr;
    win = [];
    if any(sig)
        [~, b0] = min(qv(ch, :));
        a = b0; while a > 1 && sig(a-1), a = a - 1; end
        b = b0; while b < numel(sig) && sig(b+1), b = b + 1; end
        win = [edges(a), edges(b+1)];
    end
    ti = marks.ied_times{ch};
    is_locked = false(size(ti));
    if ~isempty(win)
        for j = 1:numel(ti)
            k = tr_on == ceil(ti(j)/nS);
            d = (ti(j) - g_on(k)) / fs;
            is_locked(j) = any(d >= win(1) & d < win(2));
        end
    end
    lt{ch} = ti(is_locked); ot{ch} = ti(~is_locked);
    nb = 0; if ~isempty(win), nb = round(diff(win)/A.bin); end
    locking(ch).window = win;
    locking(ch).n_locked = sum(is_locked);
    locking(ch).n_other = sum(~is_locked);
    locking(ch).n_expected = lam(ch) * nb;
    locking(ch).counts = C(ch, :);
    locking(ch).baseline = lam(ch);
    locked_mask(ch, :) = widen(lt{ch}, mI, nS, N);
    other_mask(ch, :)  = widen(ot{ch}, mI, nS, N);
end

%% 2. saccadic spike channels
% On the 30 Hz high-passed signal: the spike potential is a fast transient,
% while the corneoretinal step every saccade also carries is slow and drops out.
w0 = round(A.sp_win(1)*fs); w1 = round(A.sp_win(2)*fs);
ok_sp = g_on + w0 >= (tr_on-1)*nS + 1 & g_on + w1 <= tr_on*nS;
gs = g_on(ok_sp);
[bh, ah] = butter(4, A.sp_hp / (fs/2), 'high');
sp_z = zeros(1, nCh); sp_es = zeros(1, nCh);
for ch = 1:nCh
    H = zeros(1, N);
    for t = 1:nTr
        H((t-1)*nS + (1:nS)) = filtfilt(bh, ah, double(data_eeg.trial{t}(ch, :)));
    end
    m = mean(H(gs + (w0:w1)), 1);
    sp_z(ch)  = max(abs(m)) / (std(H) / sqrt(numel(gs)));
    sp_es(ch) = max(abs(m)) / (1.4826 * median(abs(H - median(H))));
end
z_crit = sqrt(2) * erfcinv(A.sp_alpha / (nCh * (w1 - w0 + 1)));   % two-sided Bonferroni
% linear scale: sizes pile up near zero, where a log scale would stretch the
% lower half and push the fence far out
qs = sort(sp_es); n = numel(qs);
qq = interp1((1:n) - 0.5, qs, n*[0.25 0.75], 'linear', 'extrap');  % quartiles, as quantile()
es_fence = qq(2) + A.sp_fence_k * (qq(2) - qq(1));
sp_chan = sp_z > z_crit & sp_es > es_fence;

%% 3. bad channels, one verdict for all epoch sets
excl = marks.amp | marks.diff | marks.hf | sparse(other_mask);
L = round(A.bad_seg_len * fs); nSegT = floor(nS / L);
hit = zeros(nCh, nSegT*nTr);
for t = 1:nTr
    E = full(excl(:, (t-1)*nS + (1:nSegT*L)));
    hit(:, (t-1)*nSegT + (1:nSegT)) = squeeze(any(reshape(E, nCh, L, nSegT), 2));
end
bad_frac = mean(hit, 2)';
bad_chan = bad_frac > A.bad_chan_frac;

%% pack
marks.ied_locked = sparse(locked_mask);
marks.ied_other  = sparse(other_mask);
marks.ied_locked_times = lt;
marks.ied_other_times  = ot;
marks.locking = locking;
marks.locking_bins = ctr;
marks.saccadic_spike_z = sp_z;
marks.saccadic_spike_crit = z_crit;
marks.saccadic_spike_es = sp_es;
marks.saccadic_spike_es_fence = es_fence;
marks.saccadic_spike_channel = sp_chan;
marks.bad_frac = bad_frac;
marks.bad_channel = bad_chan;
marks.annot_params = A;

fprintf(['  annotation: %d/%d channels with saccade-locked IEDs (%d locked, %d other IEDs); ' ...
         '%d saccadic-spike channels; %d bad channels\n'], ...
        sum(~cellfun(@isempty, {locking.window})), nCh, sum([locking.n_locked]), ...
        sum([locking.n_other]), sum(sp_chan), sum(bad_chan));
end

function m = widen(ti, mI, nS, N)
    m = false(1, N);
    tr = ceil(ti / nS);
    lo = max(ti - mI, (tr-1)*nS + 1); hi = min(ti + mI, tr*nS);
    for j = 1:numel(ti), m(lo(j):hi(j)) = true; end
end


function q = apply_bh_fdr(p_in)
% apply_bh_fdr  Benjamini-Hochberg FDR correction.
%
% q = apply_bh_fdr(p_in)
%
% Accepts a p-value array of any shape. NaN values are excluded from the
% BH family and remain NaN in the output. All non-NaN values are treated
% as one family (per-channel usage: pass the full p-value matrix for one
% electrode at a time).
%
% Returns q: same shape as p_in, with BH-adjusted p-values (q-values).
% A test is significant at level alpha if q < alpha.

q    = nan(size(p_in));
mask = ~isnan(p_in);
pv   = p_in(mask);
m    = numel(pv);
if m == 0, return; end

[sorted_p, ord] = sort(pv(:));
ranks    = (1:m)';
scaled   = sorted_p .* m ./ ranks;           % m/k * p_(k)
mono     = min(1, flip(cummin(flip(scaled)))); % step-up monotonisation
q_vec    = nan(m, 1);
q_vec(ord) = mono;
q(mask)  = q_vec;
end