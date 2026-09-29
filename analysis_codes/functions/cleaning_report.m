function S = cleaning_report(data_name, opts)
% CLEANING_REPORT  Cleaning numbers and QC figures for one subject.
%
%   S = cleaning_report(data_name)
%   S = cleaning_report(data_name, opts)
%
% Summarises what the event filter and the artifact cleaning (get_artifact_marks
% -> get_events) remove for one subject, saves it, and draws two figures:
%   _cleaning_overview.png  where the losses come from: the event funnel, the
%                           per-channel exclusion by cause, a channel x trial map
%                           of the marks, the IED threshold and rates, key numbers
%   _cleaning_erps.png      six example channels chosen from the data (most
%                           excluded, highest IED rate, strongest saccadic spike,
%                           typical): saccade-locked ERP (onset by default) of every epoch, of the kept
%                           ones and of the excluded ones, and the single epochs
%                           grouped by fate
% Saccade events are those every analysis uses (onset and offset share one trial
% set): blinkSac == 0, not the fixation-0 row, complete. A channel-epoch is kept
% when noiseInfo.ok (other IEDs excluded, saccade-locked IEDs kept).
%
% Inputs:
%   data_name - column 1 of get_subject_map (e.g. 'subNS167')
%   opts      - .results_dir   (results\cleaning_report)
%               .figs_dir      (figs\cleaning_report)
%               .event_type    ('onset') event the example ERPs are locked to
%               .erp_win       ([-0.5 0.5] s), .baseline ([-0.5 -0.25] s)
%               .visible       (false) show the figures
%
% Output S (also saved as <results_dir>\<data_name>_cleaning.mat, and its
% scalar fields as one row of <results_dir>\cleaning_summary.csv, replaced if
% the subject is already there):
%   settings (marks version, IED threshold/FDR, blink window), event counts,
%   channel counts and flags, channel-epoch exclusion by cause, IED numbers,
%   image-epoch numbers, and S.per_channel (table, one row per channel).

if nargin < 2, opts = struct(); end
base        = 'D:\UJ\projects\intracPhotos\';
results_dir = getopt_field(opts, 'results_dir', fullfile(base, 'results', 'cleaning_report'));
figs_dir    = getopt_field(opts, 'figs_dir',    fullfile(base, 'figs', 'cleaning_report'));
event_type  = getopt_field(opts, 'event_type', 'onset');
erp_win     = getopt_field(opts, 'erp_win', [-0.5 0.5]);
bl_win      = getopt_field(opts, 'baseline', [-0.5 -0.25]);
vis         = getopt_field(opts, 'visible', false);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

%% data
[ep, c, E] = get_epochs(event_type, data_name);
marks = get_artifact_marks(data_name, [], []);
Ei    = get_events('image', data_name);
NI    = E.noiseInfo;
lab   = E.label(:); nCh = numel(lab);
fs    = marks.fs; nS = marks.n_samples; nTr = marks.n_trials;

is_sac   = ~isnan(c.saccade_onset_time);            % fixation-0 rows carry no saccade
fix0     = c.firstInTrial == 1;
blink    = logical(c.blinkSac) & is_sac & ~fix0;
comp     = E.eventInfo.is_complete;
before   = is_sac & ~fix0 & comp;                    % what the old filter kept (blinkSac was inert)
usable   = before & ~blink;                          % the event filter now
bad      = NI.bad_channel(:)';

%% numbers
S = struct();
S.data_name        = data_name;
S.date             = datestr(now, 'yyyy-mm-dd HH:MM');
S.marks_version    = getopt_field(marks, 'version', '');
S.ied_threshold_sd = marks.ied_threshold;
S.ied_fdr          = marks.params.ied_fdr;
S.blink_window_ms  = blink_proximity_ms();
S.fs               = fs;
S.n_channels       = nCh;
S.n_recording_trials = nTr;
S.recording_min    = nTr * nS / fs / 60;

S.n_saccades       = sum(is_sac & ~fix0);
S.n_blink_adjacent = sum(blink);
S.n_incomplete     = sum(is_sac & ~fix0 & ~blink & ~comp);
S.n_usable_events  = sum(usable);
S.pct_events_kept  = 100 * S.n_usable_events / S.n_saccades;

S.n_bad_channels        = sum(bad);
S.bad_channels          = strjoin(lab(bad), ' ');
S.n_line_noise_channels = sum(NI.line_noise_channel);
S.line_noise_channels   = strjoin(lab(logical(NI.line_noise_channel)), ' ');
S.n_saccadic_spike_channels = sum(NI.saccadic_spike_channel);
S.saccadic_spike_channels   = strjoin(lab(logical(NI.saccadic_spike_channel)), ' ');
S.n_ied_locked_channels = sum(NI.ied_locked_channel);

u = find(usable);
pc = @(M) 100 * mean(M(u, :), 1);                    % % of usable events, per channel
causes = {'amp', 'diff', 'hf', 'ied_other', 'ied_locked'};
for k = 1:numel(causes)
    v = pc(NI.(causes{k}));
    S.(['pct_' causes{k}]) = mean(v);                % pooled over channels (equal weight)
end
okp = pc(NI.ok); okps = pc(NI.ok_strict);
S.pct_ok                 = mean(okp);
S.pct_ok_strict          = mean(okps);
S.median_channel_ok_pct  = median(okp);
S.min_channel_ok_pct     = min(okp);
S.n_channels_ok_below_50 = sum(okp < 50);
S.channel_epochs_kept    = sum(sum(NI.ok(u, :)));
S.channel_epochs_usable  = numel(u) * nCh;

rate = marks.ied_rate_per_min(:)';
S.ied_total          = sum(cellfun(@numel, marks.ied_times));
S.ied_locked_total   = sum(cellfun(@numel, marks.ied_locked_times));
S.ied_other_total    = sum(cellfun(@numel, marks.ied_other_times));
S.ied_rate_median    = median(rate);
S.ied_rate_max       = max(rate);

Ni = Ei.noiseInfo; ci = Ei.eventInfo.is_complete;
S.n_image_trials     = numel(ci);
S.n_image_complete   = sum(ci);
S.pct_image_ok       = 100 * mean(mean(Ni.ok(ci, :), 1));

S.per_channel = table(lab, okp(:), okps(:), pc(NI.amp)', pc(NI.diff)', pc(NI.hf)', pc(NI.ied_other)', ...
    pc(NI.ied_locked)', rate(:), bad(:), logical(NI.line_noise_channel(:)), logical(NI.saccadic_spike_channel(:)), ...
    NI.saccadic_spike_z(:), logical(NI.ied_locked_channel(:)), ...
    'VariableNames', {'label', 'ok_pct', 'ok_strict_pct', 'amp_pct', 'diff_pct', 'hf_pct', 'ied_other_pct', ...
    'ied_locked_pct', 'ied_rate_per_min', 'bad', 'line_noise', 'saccadic_spike', 'saccadic_spike_z', 'ied_locked_channel'});

save(fullfile(results_dir, [data_name '_cleaning.mat']), 'S');
write_summary_csv(results_dir);

%% figure 1: overview
f1 = figure('Position', [20 20 1700 950], 'Visible', onoff(vis), 'Color', 'w');
tl = sprintf('%s - cleaning overview (marks %s, IED FDR %g%%, blink window %g ms)', data_name, S.marks_version, 100*S.ied_fdr, S.blink_window_ms);
sgtitle(tl, 'Interpreter', 'none', 'FontWeight', 'bold');

FS = 9;                                              % one font size for every panel

% a) event funnel
subplot(2, 3, 1);
vals = [S.n_saccades, S.n_saccades - S.n_blink_adjacent, S.n_usable_events];
names = {'saccades', {'not blink-', 'adjacent'}, {'complete', '= usable'}};
bar(1:3, vals, 0.6, 'FaceColor', [.55 .65 .85], 'EdgeColor', 'none'); hold on;
text(1:3, vals, compose('%d', vals), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FS);
lost = [NaN, S.n_blink_adjacent, S.n_incomplete];
text(2:3, vals(2:3) / 2, compose('-%d', lost(2:3)), 'HorizontalAlignment', 'center', 'FontSize', FS, 'Color', [.8 .1 .1]);
set(gca, 'XTick', 1:3, 'XTickLabel', {}); xlim([0.4 3.6]); ylim([0 max(vals) * 1.12]);
for k = 1:3, text(k, -0.03 * max(vals), names{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', FS); end
ylabel('saccade events');
title(sprintf('saccade events: %d of %d usable (%.0f%%)', S.n_usable_events, S.n_saccades, S.pct_events_kept));

% b) % of each channel's usable epochs excluded, channels sorted, causes stacked.
% An epoch with several causes is counted once, under the first in this order,
% so the stack height is the channel's total exclusion.
subplot(2, 3, [2 3]);
Xb = zeros(nCh, 5); rest = true(numel(u), nCh);
pri = {bad, 'amp', 'diff', 'hf', 'ied_other'};
for k = 1:5
    if k == 1, M = repmat(bad, numel(u), 1); else, M = NI.(pri{k})(u, :); end
    Xb(:, k) = 100 * mean(M & rest, 1)';
    rest = rest & ~M;
end
[~, o] = sort(sum(Xb, 2), 'descend');               % also orders the marked-time map
% Drawn with IED other at the bottom, so the rarer causes sit as caps on the
% outline instead of hiding as thin lines inside it (the counting order above
% is unchanged).
cause_names = {'bad channel', 'amp', 'diff', 'hf', 'IED other'};
cols = [0 0 0; .2 .45 .8; .85 .4 .1; .95 .7 .1; .55 .25 .75];
po = [5 4 3 2 1];
b = bar(Xb(o, po), 'stacked', 'BarWidth', 1, 'EdgeColor', 'none'); hold on;
for k = 1:5, b(k).FaceColor = cols(po(k), :); end
flagrow = {logical(NI.saccadic_spike_channel(:)), 'm', 'saccadic-spike channel'; ...
           logical(NI.ied_locked_channel(:)),     'g', 'channel with saccade-locked IEDs'};
hflag = gobjects(0); ynext = -4;
for k = 1:size(flagrow, 1)
    ii = find(flagrow{k, 1}(o));
    if isempty(ii), ii = NaN; end                    % keep the legend entry when no channel has the flag
    hflag(end+1) = plot(ii, ynext * ones(size(ii)), '^', 'Color', flagrow{k, 2}, 'MarkerFaceColor', flagrow{k, 2}, 'MarkerSize', 4); %#ok<AGROW>
    ynext = ynext - 4;
end
legend([fliplr(b) hflag], [cause_names, flagrow(:, 3)'], 'Location', 'northeast');   % legend top = top of the stack
xlim([0 nCh + 1]); ylim([ynext + 1 100]); set(gca, 'YTick', 0:10:100);   % fixed 0-100 %, comparable across subjects
xlabel('channel (sorted by % excluded)'); ylabel('% of the channel''s usable epochs excluded');
title(sprintf('epochs excluded per channel: mean %.1f%%, median %.1f%%, max %.1f%%; %d of %d channels > 50%%', ...
      100 - S.pct_ok, 100 - S.median_channel_ok_pct, 100 - S.min_channel_ok_pct, S.n_channels_ok_below_50, nCh));

% c) channel x recording trial map of marked time, channels in the order of b)
subplot(2, 3, 4);
Mk = logical(marks.amp) | logical(marks.diff) | logical(marks.hf) | logical(marks.ied_other);
Bk = kron(speye(nTr), ones(nS, 1));
frac = full(double(Mk) * Bk) / nS;
imagesc(100 * frac(o, :)); set(gca, 'CLim', [0 20]); colormap(gca, flipud(gray)); cb = colorbar; cb.Label.String = '% of trial marked';
xlabel('recording trial'); ylabel('channel (sorted as in the bar plot)');
title('marked time (amp, diff, hf, IED other)');

% d) IED threshold
subplot(2, 3, 5);
C = marks.ied_fdr_curve;
semilogy(C.grid, C.observed, 'b', C.grid, C.null, 'k--'); hold on;
xline(marks.ied_threshold, 'r', sprintf('threshold %.2f', marks.ied_threshold));
xlim([C.grid(1) min(C.grid(end), 20)]);
xlabel('IED detector score'); ylabel('detections per channel-minute');
legend({'observed', 'phase-randomised null'});
title(sprintf('IEDs: %d detected (%d saccade-locked, %d other); rate median %.2f, max %.1f /min', ...
      S.ied_total, S.ied_locked_total, S.ied_other_total, S.ied_rate_median, S.ied_rate_max));

% e) key numbers
subplot(2, 3, 6); axis off;
txt = {sprintf('channels: %d (bad %d: %s)', nCh, S.n_bad_channels, emptydash(S.bad_channels)), ...
       sprintf('saccadic-spike channels: %d', S.n_saccadic_spike_channels), ...
       sprintf('channels with saccade-locked IEDs: %d', S.n_ied_locked_channels), '', ...
       sprintf('saccades %d -> usable %d (%.0f%%)', S.n_saccades, S.n_usable_events, S.pct_events_kept), ...
       sprintf('   blink-adjacent %d, incomplete %d', S.n_blink_adjacent, S.n_incomplete), '', ...
       'epochs flagged (% of usable, any cause):', ...
       sprintf('   amp %.2f | diff %.2f | hf %.2f', S.pct_amp, S.pct_diff, S.pct_hf), ...
       sprintf('   IED other %.2f | IED locked %.2f (kept)', S.pct_ied_other, S.pct_ied_locked), ...
       sprintf('   kept: ok %.1f%%, ok\\_strict %.1f%%', S.pct_ok, S.pct_ok_strict), '', ...
       sprintf('image epochs: %d/%d complete, %.1f%% kept per channel', S.n_image_complete, S.n_image_trials, S.pct_image_ok), ...
       sprintf('recording: %d trials, %.1f min', nTr, S.recording_min)};
text(0, 1, txt, 'VerticalAlignment', 'top', 'FontSize', FS, 'FontName', 'Consolas');

set(findall(f1, 'Type', 'axes'), 'FontSize', FS, 'TitleFontSizeMultiplier', 1.1, 'LabelFontSizeMultiplier', 1, 'TitleFontWeight', 'bold');
set(findall(f1, 'Type', 'legend'), 'FontSize', FS - 1);
exportgraphics(f1, fullfile(figs_dir, [data_name '_cleaning_overview.png']), 'Resolution', 110);
if ~vis, close(f1); end

%% figure 2: example channels, before vs after
okp_all = 100 * mean(NI.ok(before, :), 1);           % includes the blink loss
cand = find(~bad);
pick = []; why = {};
[~, o] = sort(okp_all(cand), 'ascend');
for k = 1:2, pick(end+1) = cand(o(k)); why{end+1} = 'most excluded'; end %#ok<AGROW>
r = rate; r(pick) = -Inf; r(bad) = -Inf; [~, k] = max(r);
pick(end+1) = k; why{end+1} = 'highest IED rate';
z = NI.saccadic_spike_z(:)'; z(pick) = -Inf; z(bad) = -Inf; [~, k] = max(z);
pick(end+1) = k; why{end+1} = 'strongest saccadic spike';
rem = setdiff(cand, pick); [~, o] = sort(abs(okp_all(rem) - median(okp_all(cand))));
for k = 1:2, pick(end+1) = rem(o(k)); why{end+1} = 'typical'; end %#ok<AGROW>
if any(bad), pick(end) = find(bad, 1); why{end} = 'bad channel'; end

t   = ep.time; ti = t >= erp_win(1) - 1e-9 & t <= erp_win(2) + 1e-9; tt = t(ti) * 1000;
bli = t(ti) >= bl_win(1) - 1e-9 & t(ti) <= bl_win(2) + 1e-9;
eb  = find(before);
cA = [.8 .1 .1]; cB = [.95 .6 .1]; cK = [.1 .3 .8];  % artifact/IED, blink-adjacent, kept
xlab = sprintf('time from saccade %s (ms)', event_type);
f2 = figure('Position', [20 20 1800 900], 'Visible', onoff(vis), 'Color', 'w');
sgtitle(sprintf('%s - %s-locked ERPs before and after cleaning (baseline %g..%g ms)', data_name, event_type, 1000*bl_win), ...
        'Interpreter', 'none', 'FontWeight', 'bold');
for k = 1:numel(pick)
    ch = pick(k);
    X = cell2mat(cellfun(@(x) double(x(ch, ti)), ep.trial(eb), 'UniformOutput', false)') * 1e6;   % V -> uV
    X = X - mean(X(:, bli), 2);
    kept  = usable(eb) & NI.ok(eb, ch);
    exA   = usable(eb) & ~NI.ok(eb, ch);             % excluded by the channel's marks
    exB   = blink(eb);                               % excluded as blink-adjacent
    subplot(2, numel(pick), k); hold on;
    h = gobjects(0); lg = {};
    [h(end+1), lg{end+1}] = erpline(tt, X, true(size(kept)), [.6 .6 .6], 2.5, 'all (before)');
    if sum(exA) >= 3, [h(end+1), lg{end+1}] = erpline(tt, X, exA, cA, 1, 'excluded: artifact/IED'); end
    if sum(exB) >= 3, [h(end+1), lg{end+1}] = erpline(tt, X, exB, cB, 1, 'excluded: blink-adjacent'); end
    [h(end+1), lg{end+1}] = erpline(tt, X, kept, cK, 1.5, 'kept (after)');
    xline(0, ':'); xlim(tt([1 end])); xlabel(xlab);
    title({sprintf('%s (%s)', lab{ch}, why{k}), sprintf('kept %d of %d', sum(kept), numel(kept))}, 'Interpreter', 'none', 'FontSize', 8);
    if k == 1, ylabel('\muV (mean \pm SEM)'); end
    legend(h, lg, 'FontSize', 6, 'Location', 'best', 'Box', 'off');

    subplot(2, numel(pick), numel(pick) + k);
    grp = 1*exA + 2*(exB & ~exA) + 3*kept;           % same order as the lines above
    [gs, oo] = sort(grp);
    sd = 1.4826 * median(abs(X(:) - median(X(:))));
    w  = 0.05 * (tt(end) - tt(1));                   % width of the side bar
    imagesc(tt, 1:numel(oo), X(oo, :) / sd); set(gca, 'CLim', [-4 4]); colormap(gca, redblue());
    hold on;
    gcol = {cA, cB, cK};
    for g = 1:3
        r = find(gs == g);
        if isempty(r), continue; end
        patch(tt(1) - w + [0 w w 0], [r(1) r(1) r(end) r(end)] + [-.5 -.5 .5 .5], gcol{g}, 'EdgeColor', 'none');
    end
    for g = 1:2                                      % separators between the groups
        e = find(gs <= g, 1, 'last');
        if ~isempty(e) && e < numel(gs) && any(gs == g), yline(e + 0.5, 'k', 'LineWidth', 1); end
    end
    xline(0, ':k'); xlim([tt(1) - w, tt(end)]); set(gca, 'YDir', 'reverse');
    n = arrayfun(@(g) sum(gs == g), 1:3);
    title(sprintf('epochs: artifact/IED %d | blink %d | kept %d', n(1), n(2), n(3)), 'FontSize', 7);
    xlabel(xlab); if k == 1, ylabel('epoch (grouped)'); end
    if k == numel(pick)                              % colorbar beside the last panel, without shrinking it
        pos = get(gca, 'Position');
        cb = colorbar; set(gca, 'Position', pos); drawnow;
        pos = get(gca, 'Position');
        cb.Position = [pos(1) + pos(3) + 0.006, pos(2), 0.008, pos(4)];
        cb.Label.String = 'robust SD (scaled median absolute deviation, MAD)';   % 1.4826 x MAD, notes/trial_cleaning.md sec. 4
    end
end
exportgraphics(f2, fullfile(figs_dir, [data_name '_cleaning_erps.png']), 'Resolution', 110);
if ~vis, close(f2); end

fprintf('cleaning report %s: %d/%d events usable, %.1f%% channel-epochs kept; saved to %s\n', ...
        data_name, S.n_usable_events, S.n_saccades, S.pct_ok, results_dir);
end

%% helpers
function [h, lg] = erpline(tt, X, sel, col, lw, name)
    Y = X(sel, :); m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') / sqrt(size(Y, 1));
    fill([tt fliplr(tt)], [m + se, fliplr(m - se)], col, 'FaceAlpha', .15, 'EdgeColor', 'none');
    h = plot(tt, m, 'Color', col, 'LineWidth', lw);
    lg = sprintf('%s (%d)', name, size(Y, 1));
end

function write_summary_csv(results_dir)
% The CSV is rebuilt from every <data_name>_cleaning.mat in the folder, so it is
% never read back (readtable would turn the date column into datetimes) and a
% rerun of one subject replaces just its row. Fields missing from an older file
% are left empty.
    d = dir(fullfile(results_dir, '*_cleaning.mat'));
    rows = cell(numel(d), 1);
    for i = 1:numel(d)
        L = load(fullfile(d(i).folder, d(i).name), 'S'); S = L.S;
        fn = fieldnames(S);
        keep = fn(cellfun(@(x) (isnumeric(S.(x)) && isscalar(S.(x))) || ischar(S.(x)), fn));
        rows{i} = rmfield(S, setdiff(fn, keep));
    end
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows), return; end
    fl = cellfun(@fieldnames, rows, 'UniformOutput', false);
    names = unique(vertcat(fl{:}), 'stable');
    for i = 1:numel(rows)
        for k = 1:numel(names)
            if ~isfield(rows{i}, names{k}), rows{i}.(names{k}) = ''; end
        end
        rows{i} = orderfields(rows{i}, names);
    end
    T = struct2table(vertcat(rows{:}), 'AsArray', true);
    writetable(sortrows(T, 'data_name'), fullfile(results_dir, 'cleaning_summary.csv'));
end

function s = onoff(v)
    if v, s = 'on'; else, s = 'off'; end
end

function s = emptydash(s)
    if isempty(s), s = '-'; end
end

function m = redblue()
    n = 128; a = linspace(0, 1, n)';
    m = [[a a ones(n,1)]; [ones(n,1) flipud(a) flipud(a)]];
end
