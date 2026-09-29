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
