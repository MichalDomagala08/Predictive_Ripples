function marks = get_artifact_marks(subjectId, data, content, overwrite)
% GET_ARTIFACT_MARKS  A subject's artifact marks, computed once and cached.
%
%   marks = get_artifact_marks(subjectId, data, content)
%   marks = get_artifact_marks(subjectId, data, content, overwrite)
%
% mark_artifacts + annotate_marks on the reref trials, saved to
% data_path('artifact_marks')artifact_marks_<subjectId>.mat and read back on later
% calls, so every epoch set of a subject (onset, offset, image) is cleaned
% against the same marks, not merely identically computed ones.
%
% Inputs:
%   subjectId - e.g. 'subns127_02'
%   data      - loaded reref file (data.data_eeg); may be [] when the cache exists
%   content   - saccade table (saccade_onset_time, trial_number); may be [] then too
%   overwrite - recompute even if cached (default false)

if nargin < 4, overwrite = false; end
% Bump when mark_artifacts / annotate_marks change what they mark: a cache
% written by other code is then recomputed instead of silently reused.
MARKS_VERSION = '2026-09-24d';   % IED threshold at FDR 1 % (was 5 %); + saccade-table fingerprint
markdir = 'D:\Documents_Dell\Predictive_Ripples_2025\preprocessed\RippleDensity_Franz_improved_saccadic\marks';
fpath   = [markdir 'artifact_marks_' subjectId '.mat'];

if ~overwrite && exist(fpath, 'file') == 2
    S = load(fpath, 'marks');
    if strcmp(getopt_field(S.marks, 'version', ''), MARKS_VERSION)
        marks = S.marks;
        fprintf('  artifact marks read from %s\n', fpath);
        return;
    end
    fprintf('  cached marks are from other code (%s), recomputing\n', getopt_field(S.marks, 'version', 'no version'));
end
if isempty(data) || isempty(content)
    error('No cached marks for %s: pass the reref data and the saccade content.', subjectId);
end
if ~exist(markdir, 'dir'), mkdir(markdir); end

t_m = tic;
marks = mark_artifacts(data);
marks = annotate_marks(marks, data, content);
marks.subject = subjectId;
marks.version = MARKS_VERSION;
% The annotation used this saccade table: a re-extracted table must re-mark
% (get_events compares this fingerprint with the current table).
marks.content_fingerprint = content_fingerprint(content);
save(fpath, 'marks', '-v7.3');
fprintf('  artifact marks computed in %.0f s, saved to %s\n', toc(t_m), fpath);
end

function fp = content_fingerprint(content)
% Row count and sum of saccade onset times: changes whenever the table is
% re-extracted, costs nothing to compute. get_events recomputes it the same way.
    fp = [height(content), sum(content.saccade_onset_time, 'omitnan')];
end
