# Analysis Pipeline — `run_*` scripts

All scripts live in `mfiles/`. Run them from the project root or from within MATLAB after `cd` to the project root (each script does this automatically).

`mfiles/` must be on the MATLAB path: the scripts `cd` to the project root and then call helpers by name, without `addpath`.

### Naming convention

| Prefix | Role | Writes to |
|---|---|---|
| `run_*` | Analysis — compute a measure and save it | `data/` or `results/` |
| `run_plot_*` | Visualisation — read pre-computed results and draw | `figs/` |

The split matters because the two are run at different times and for different reasons: analysis is slow and rerun rarely, plotting is fast and rerun constantly while tuning a figure. Where a plotting driver has a matching plotter function, the names mirror each other — `run_plot_erp_by_region.m` drives `plot_erp_by_region()`. The exceptions are `run_plot_split_vs_reg.m` (drives `plot_erp_by_region_multi` with `show_regression = true`) and the brain-map and percentage-significant drivers, which have no single `plot_*` counterpart.

### Shared helpers

| Helper | Purpose |
|---|---|
| `get_subject_map()` | Canonical 21×4 subject table: `data_name \| folder_name \| pic_name \| elec_locs_name`. Every `run_*` / `run_plot_*` script gets its subject list from here — add subjects in this one file, and append rather than reorder, since row indices are used as `subjects_to_run`. |
| `compute_itc(...)` | Wavelet ITC for one channel; also returns the per-trial phasors used by the cached-ITC pipeline. |
| `test_per_timepoint(mat1[, mat2][, alpha])` | Per-timepoint significance. One matrix → one-sample t-test vs 0; two → Wilcoxon rank-sum. Always returns `[p_vals, sig_mask]`. |
| `apply_bh_fdr(p_in)` | Benjamini-Hochberg correction (see Stage 1b). |
| `getopt_field(s, field, default)` | `s.(field)` if present and non-empty, else `default`. Used to give `opts`-struct arguments their defaults. |
| `data_path(type)` | The folder each data type is read from and written to: `raw`, `reref`, `artifact_marks`, `power_bha`, `power_and_phases_itc`. The one place that says which disc holds what: to move a type, move the folder and edit its line. Copies on the other disc are backups the code does not read. The per-epoch stores (`pad_perisac`, `bha_power`, `itc_phases`, `perisac`, `wide_perisac*`) were retired on 2026-09-24; asking for them is an error. |
| `get_events(event_type, data_name[, use_rand])` | One subject's events of one type (`onset` / `offset` / `image`, or the randomised control), each located in the continuous reref trials (`eventInfo.sample`), with its artifact flags (`noiseInfo`) — without loading any EEG. See *Events and epochs on demand*. |
| `get_epochs(event_type, data_name[, use_rand])` | The `data_perisac` struct the stored epoch files used to hold, cut on demand from the reref trials; also returns the event table and `get_events`' output. Every epoch reader calls this. |
| `get_power_bha(E, ch, time_rel[, rows])` / `get_bha_channel(E, ch, rows)` | BHA of chosen events from the continuous `power_bha` cache, with the per-epoch edge NaNs and detrend of the old cache (the second returns the old `.valid_idx` / `.bha_mat` struct). |
| `open_power_and_phases_itc(data_name, elec_name, event_type, dataset, time_win)` | The per-event view of the continuous wavelet cache, in the old `itc_phases` shape (`.channel(ch)` → `.valid_idx`, `.itc`, `.phases`, `.power`). |
| `get_artifact_marks(subj, data, content)` | A subject's artifact marks (`mark_artifacts` + `annotate_marks`), computed once and cached in `data/artifact_marks/` (see *Artifact detection* below). |
| `mark_artifacts(data_eeg[, opts])` | Sample-level artifact and IED marks on the continuous reref trials, blind to behaviour. |
| `annotate_marks(marks, data_eeg, content[, opts])` | Saccade-locked vs other IEDs, saccadic-spike channels, bad channels. |
| `detect_ieds(z, fs[, opts])` | IED detector of Michal's `findSpikeTimes_reviewed`, both polarities; returns each detection's score. |
| `selection_label(opts)` | A readable name for an electrode selection (`region ITC+antHP`, `list test_set`) and a filename-safe tag for it. Reads the same `selection_mode` fields as `rsa_resolve_targets`; used by the group summaries, which draw one figure per electrode **set** and so must name that set. |
| `return_saccade_defaults()` | The one place the return-saccade threshold lives: `thresh_px` (78 px = 2 dva at 39.2 px/deg), `px_per_deg`, `thresh_dva`. `add_return_saccades_to_content` takes its `opts.thresh_px` default from here, and the driver, the preview figure and the CSV exporter all read it, so a table built at one threshold can never be paired with a figure drawn at another. Changing it means the stored columns are stale — rerun `add_return_saccades.m` with `overwrite = true`, then `export_return_saccades_csv.m`. |
| `get_pics_lookups([pics_base])` | The two electrode-image maps (`elec_name -> PICS dir`, `elec_name -> pic_name`) that every `run_plot_*` driver showing anatomy passes to its plotter. The freesurfer root lives here only — do not re-hardcode it in a driver. |
| `get_cond_field_map([name])` | Short predictor name (`aws`, `lumAbsDiff`, `amp`, …) → column name in the content table; errors on an unknown name. No argument returns the whole map. Every script that resolves a `cond_var` (or `opts.predictors`) goes through it — add a variable here, not in a driver. |
| `make_condition_split(content, cond_var, split_type)` | Assigns each content row to a condition: returns `cond_vec` (condition index per row, NaN = in none) and `condition_names`. Holds the `tertile` / `median` / `quantile<N>` / `flags` rules for `run_erp_per_subject.m`; called with an empty `content` it only validates the names, so the driver errors before loading a subject. (`run_itc_by_condition.m` and `run_itc_band_by_condition.m` still carry their own copies of the first three rules and do not know `flags`.) |
| `cmap_rdbu([n])` | ColorBrewer RdBu diverging colormap (blue − white − red) for data on a symmetric scale around 0. Used by the ITC difference panels, the RdBu brain maps, the RSA similarity matrices and the regression beta strips. |
| `brain_map_tag(opts)` | Canonical file and folder names for brain-map outputs. Both `make_brain_map` and `run_plot_brain_map_animation` name files through it, so rendering and lookup cannot diverge. |

#### Shared by the by-region plotter family

All four by-region plotters — `plot_erp_by_region`, `plot_erp_by_region_multi`, `plot_erp_overall` and `plot_regression_by_region` — differ only in how they load stats and lay out one row. Everything around that is shared, so a change to selection, output naming or paging belongs in one of these:

| Helper | Purpose |
|---|---|
| `region_plot_begin(elec_locs_path, region_col, output_dir, options)` | Stage 1: normalises `region_col`, resolves `region_name` / `region_name_safe`, selects the subject-electrode pairs, and returns the output paths. Returns `reg` with `.pairs .region_name .region_name_safe .pages_dir .combined_dir`; empty `.pairs` means nothing to plot. Directories are **not** created here, so a run that finds no data leaves no empty folders. |
| `region_plot_finish(page_pdfs, combined_dir, combined_name, region_name, n_pairs)` | Stage 3 epilogue: merges the per-page PDFs (or copies a single page) and prints the closing summary. |
| `cond_suffix_str(cond_var, split_type)` | Filename suffix `_<cond_var>_<split_type>`, or `''` when either is empty. |
| `new_page_figure(page_title, n_rows, n_tile_cols)` / `save_page_pdf(fig, path, pg, n_pages)` | Bracket one PDF page: off-screen landscape figure with a tiled layout, then print, close and report. |
| `load_result_cached(cache, mat_path, mat_name, label, warn_missing)` | Stage 2: load the `result` variable from a per-subject `.mat` once, memoising misses as `[]` so each file is read (and reported) only once regardless of how many electrodes reference it. |
| `unpack_regression_stats(ch, time, alpha_sig)` | Pull `.r2 .beta .pvals .sig_mask .n` out of a channel's regression block; `[]` when the channel has none. |

#### ITC panel renderers

`render_itc_panel` (mean resultant length), `render_itc_diff_panel` (high − low, diverging RdBu) and `render_rayleigh_mask_panel` (binary significance) are thin wrappers over `render_itc_map`, which handles the empty-data branch, colormap, colour limits, colorbar, labels and title. `render_itc_map` draws through `itc_imagesc`, which owns the frequency axis: `'linear'` resamples the log-spaced rows onto a uniform Hz grid via `itc_interp_linear` (use `'nearest'` for binary masks), `'log'` plots raw rows against pixel index and relabels the ticks. Add a new ITC panel type as another wrapper, not another copy.

---

## Stage 0 — Preprocessing

```
run_preprocess.m
```

Covers the steps below (all in one file, sections commented out as needed). Since 2026-09-24 nothing is stored per epoch or per event type: Stage 0 ends with the continuous reref trials and the saccade tables, and everything event-locked is cut from them on demand.

| Step | Action | Input | Output |
|------|--------|-------|--------|
| 1 | Parse eye-tracking, compute saccade content: `make_contentEKM` → `ekmToContentPhotos` (Engbert–Kliegl detector `m2b3_ekm`, PSO handling `dealWithPSO`, both copied into `mfiles/` on 2026-09-24 from the McGill toolboxes). Gaze is smoothed by a causal 5-sample moving average over valid samples only (fixed 2026-09-24: the old zero-state `filter` faked a first saccade at ~3 ms and dropped fixation 0 in 75 % of trials — `notes/first_fixation_filter_fix.md`). `blinkSac` = an invalid eye sample within `blink_proximity_ms()` of the saccade, distance in `dataloss_dist_ms` (since 2026-09-24; the EyeLink-era test for a blink *inside* the saccade could never fire on this Tobii data — window derived and existing tables patched by `standalone/add_blink_proximity.m`). The derived columns are re-added by `standalone/rebuild_content.m` | `data_path('raw')` (`data.data_et`) | `data/sacextr_ekm/content_EKM_*.mat` |
| 2 | Add image statistics (AWS / luminance / edge) to content | content .mat | content .mat (appended) |
| 2b | Bipolar re-reference, line-noise removal and trial demeaning in one pass (script step "2)"): `preprocess_bipolar(iSubj, rawdir, rerefdir, 5)` re-references, calls `cleanTrialEEG` — 4th-order Butterworth band-stop, zero-phase, ±5 Hz at 60/120/180/240 Hz — then removes each recording trial's mean per channel (`cfg.trial_demean`), and saves each subject once (temporary name, renamed when complete). A subject whose file is already notched and demeaned is skipped, so an interrupted run can be restarted. (Before 2026-09-24: a separate ±2 Hz pass, no demeaning.) | `data_path('raw')` | `data_path('reref')`, overwritten |
| — | Steps 3–4 (cut and clean stored epoch sets: `wide_perisac`, then `run_pad_perisac` / `run_image_epochs`) are **retired**; see below | | |

**Why demean each trial.** The analyses z-score against one pooled baseline mean and SD per channel. Channel offsets differ between recording trials (in `subNS167` 71 % of channels shift by more than a within-trial SD, with steps at block boundaries); left in, they inflate that SD and add noise to every single-trial value, and they bias `abs_erp` and any predictor that correlates with trial order. The whole-trial mean is removed, not a linear trend, which could bend slow image-locked responses.

### Events and epochs on demand — `get_events` / `get_epochs`

```
get_events(event_type, data_name[, use_rand])   events + where each sits in the continuous trials + artifact flags (no EEG loaded)
get_epochs(event_type, data_name[, use_rand])   the data_perisac struct the stored files used to hold, cut from data_path('reref')
```

The stored epoch sets (`pad_perisac`, 53 GB) were retired on 2026-09-24. Every reader calls `get_epochs` instead of loading a file, and gets the same struct: `.trial` (single), `.time`, `.label`, `.electrodes`, `.eventInfo` (with `is_complete` and, new, `sample`, the event's index in the concatenated trials), `.noiseInfo`, `.params`. The rules are the ones `run_pad_perisac` / `run_image_epochs` applied, reproduced sample-for-sample (checked with `standalone/check_redesign.m`):

- **Saccade events** (`'onset'` / `'offset'`): each event is snapped to the nearest sample of its recording trial; the view is ±1.25 s around it (the analyses read ±0.5 s; the rest is spectral padding — BHA 0.05 s, ITC 0.48 s at 3 cycles / 3 Hz). The artifact flags are computed once over the **union** `[onset − 1.25 s, offset + 1.25 s]`, so onset and offset results rest on the same trial set; a saccade longer than 100 ms, or without a duration, is incomplete (0.14 % of events).
- **Image events** (`'image'`): one per recording trial, at image onset, view [−1.25 2.25] s (the analyses read [−0.5 1.5], `time_win_image`). The event table is a per-trial table (`trial_number`, `name`, `event_time`), so the saccade-level filters do not apply: the event filter is `is_complete` ∩ `noiseInfo.ok`. `run_erp_per_subject` writes `<subj>_<measure>_overall.mat` for it (no condition split).
- **Wider is not free**: a channel-epoch is flagged if *any* of its samples is marked, so rejection grows with the window. Measured with the version-1 detector on one subject: 90.4 % of channel-epochs usable at ±0.5 s, 88.5 % at ±1.05, 85.7 % at ±2.0.

### Artifact detection (version 2, 2026-09-24)

```
get_artifact_marks(subj, data, content)  ──►  data_path('artifact_marks')artifact_marks_<subj>.mat
   = mark_artifacts (blind to behaviour)  +  annotate_marks (saccade-related labels, bad channels)
```

Bad samples are marked once per subject on the continuous reref trials and cached; `get_events` builds the cache on first use, and rebuilds it when the saccade table it was annotated with has changed (a fingerprint of the table is stored with the marks) or the marking code has (`MARKS_VERSION`). Every event set of that subject (onset, offset, image) is flagged against the same marks, and bad channels are the same in all of them.

The artifact steps (amplitude, first difference, FIR high-pass at 225 Hz) share one rule: per channel, a sample is marked where the step's statistic exceeds the Tukey far-out fence (Q3 + 3·IQR) of the log per-second maxima, computed on trial-demeaned data. (A 55–65 Hz burst step was dropped on 2026-09-24: the ±5 Hz notch of step 2b removes the transient line-noise bursts it was catching.) IEDs come from `detect_ieds` (Michal's detector, both polarities), with a per-subject threshold at FDR ≤ 1 % against a phase-randomised null (5 % until 2026-09-24).

`annotate_marks` splits IEDs into `ied_locked` (inside the channel's saccade-locked window, estimated from the correlogram) and `ied_other`. The default `noiseInfo.ok` excludes artifacts and other IEDs; `ok_strict` also excludes locked IEDs. `annotate_marks` also sets the channel flags `saccadic_spike_channel`, `line_noise_channel` and `bad_channel`.

Full description, reasons and the effect on three recordings: `notes/trial_cleaning.md`.

**Cleaning report — `run_cleaning_report.m` → `cleaning_report(data_name, opts)`.** Per subject: what the event filter and the cleaning remove, saved for reporting, and two QC figures for eyeballing. Output: `results/cleaning_report/<data_name>_cleaning.mat` (struct `S` with the settings — marks version, IED threshold and FDR, blink window — event counts, channel flags, channel-epoch exclusion by cause, IED numbers, image-epoch numbers, and `S.per_channel`, one row per channel); `results/cleaning_report/cleaning_summary.csv` (the scalar fields, one row per subject, replaced on rerun); `figs/cleaning_report/<data_name>_cleaning_overview.png` (event funnel, per-channel exclusion by cause with flagged channels, channel × recording-trial map of marked time, IED threshold curve, key numbers) and `_cleaning_erps.png` (six channels picked from the data — two most excluded, highest IED rate, strongest saccadic spike, two typical, or the bad channel — offset-locked ERP of all / kept / excluded epochs, and the single epochs grouped by fate). ~40 s per subject once its marks exist.

### Randomised-control dataset

```
random_sac.m
```

Creates the shuffled event times the control analyses are compared against: for each subject and trial, replaces `saccade_onset_time` / `saccade_offset_time` with uniform random values in [0, 6] s, preserving the NaN for the first (index-0) entry per trial. Writes `data/sacextr_rand/`. `get_events(..., use_rand = 1)` reads it, so `use_rand = 1` in the drivers and the ITC `ctrl` dataset need nothing else — no separate epoch set or cache.

⚠️ **The current `data/sacextr_rand/` tables were made from the saccade tables before the 2026-09-24 filter fix** (e.g. `subNS128_02`: 609 rows vs 515 now). Rerun `random_sac.m` before any control analysis.

---

## Stage 1 — Compute neural measures (per subject)

Pipelines A and B produce the same output format and can feed the same visualisation scripts. Pipelines C and D are regression-based, each with its own output files and plotting scripts: C regresses on one continuous predictor, D on several at once.

### Continuous spectral caches — `power_bha`, `power_and_phases_itc`

```
run_power_bha.m             ──►  data_path('power_bha')power_bha_<data_name>.mat                       (~0.5 GB per subject)
run_power_and_phases_itc.m  ──►  data_path('power_and_phases_itc')power_and_phases_itc_<data_name>.mat  (~35 MB per channel, ~180 GB in all)
```

One cache per subject, computed on the whole 14 s trials and shared by every event type and by the randomised control — replacing the per-epoch, per-event-type `bha_power` (19 GB) and `itc_phases` (~270 GB **per event type**, and again for the control) caches, retired on 2026-09-24.

- **`power_bha`**: `bha{ch}`, trials × samples — `ft_specest_mtmconvol`, Hanning, 0.1 s window, 70:3:150 Hz, power averaged over frequency, as `compute_bha_stats`. Read with `get_power_bha(E, ch, time_rel, rows)`, which cuts the events and applies what the per-epoch estimate did on top (the edge samples of the padded window set to NaN, a per-epoch linear detrend), so readers get the old cache's values; `get_bha_channel` returns them in the old struct (`.valid_idx`, `.bha_mat`) over the old cache's window, for `build_brain_feature_mat` / `mreg_fit`.
- **`power_and_phases_itc`**: `coef{ch}`, trials × 15 frequencies × samples, complex single — `ft_specest_wavelet`, 3 cycles, 3–30 Hz log-spaced, as `compute_itc` — over [−0.75 6.85] s of each trial (saccades fall within 0–6 s; windows up to ±0.75 s fit). Read with `open_power_and_phases_itc(data_name, elec_name, event_type, dataset, time_win)`, which returns the old cache's shape (`.time`, `.freqs`, `.n_channels`, `.channel(ch)` → `.valid_idx`, `.itc`, `.phases`, `.power`). Because the transform ran on the whole trial, no epoch edge reaches the window any more.
- BHA readers: `run_erp_per_subject`, `run_regression_per_subject`, `mreg_fit`, `rsa_brain`. ITC readers: `run_itc_vs_rand`, `run_itc_by_condition`, `run_itc_band_by_condition` (`run_erp_per_subject` computes ITC itself from the epochs).

### Pipeline A — Simple (ERP / BHA / ITC in one script)

```
run_erp_per_subject.m
```

Set `measure_type` at the top:

```
'erp'      → results/<eventType>/erp/
'abs_erp'  → results/<eventType>/abs_erp/
'bha'      → results/<eventType>/bha/
'itc'      → results/<eventType>/itc/
```

Set `use_rand = 1` to process the randomised-control dataset instead of real data (output goes to `…/<measure_type>_rand/`).

Epochs come from `get_epochs` (cut on demand from the reref trials, ±1.25 s around a saccade event, [−1.25 2.25] s around image onset) — there is no epoch-set choice any more. Onset and offset share one trial set.

`cond_var = ''` turns the split off entirely: every trial surviving the event filter goes into one `'all'` condition and the file is `<subj>_<measure>_overall.mat`, as for `eventType = 'image'`. The saccade event filter (`blinkSac`, `firstInTrial`, `is_complete`) and the per-channel `noiseInfo.ok` still apply, so it is the unsplit counterpart of the splits below and not a different trial set. `split_type` is then ignored.

`split_type` is `'tertile'` (bottom vs top third, middle dropped), `'median'`, or `'quantile<N>'` — N equal-count bins of `cond_var`, bin 1 lowest, no trials dropped. The bins are there to *plot* a graded response (see `run_plot_erp_by_region.m`), not to test one: the between-condition Wilcoxon runs only for two conditions, so with more than two `cond_pvals` stays empty and no significance bar is drawn. `overall` (the t-test vs zero) is unaffected. The bin count travels in the filename (`…_erp_amp_quantile5.mat`), so a plotting driver reaches the file by setting the same `split_type` string; the drivers that assume a low / high **pair** — the brain maps, the percentage-significant scripts, `run_plot_erp_by_region_multi.m`, `run_plot_split_vs_reg.m` — should stay on `tertile` / `median`. Ties at a quantile edge merge bins, with a warning and correspondingly fewer condition names.

`split_type = 'flags'` compares **event classes** rather than bins of a continuous variable: `cond_var` names two or more 0/1 columns joined by `+`, and each column becomes one condition, in the order given — `cond_var = 'returned_to+return_1back'` gives the fixation that is later returned to (condition 1) and the fixation that returns to it (condition 2), one line each; any of the eight flags `add_return_saccades.m` writes can be listed, and `cond_var = 'returned_to+return_any'` pools the three return depths against their anchors. The conditions are made exclusive: a row with more than one listed flag set, or with another listed flag NaN (undefined), is dropped, so the between-condition Wilcoxon compares independent samples. For the returned-to / 1-back pair that is ~11% of the flagged rows — a fixation is often both a return and the anchor of a later one — and more for the pooled `returned_to+return_any` pair, since deeper returns overlap the anchors more. The epoch is the saccade *into* the flagged fixation, at either `eventType`. With two flags the result is a pair like `tertile`, so every downstream driver reads it (condition 2 minus condition 1 in difference panels); the file is `…_bha_returned_to+return_1back_flags.mat`, and the plotting drivers take the same two strings.

ITC inside this script runs wavelet convolution in-memory for each subject — suitable for a single subject or when RAM is not a bottleneck.

For `measure_type='itc'` (real data only): also writes an overall result to `results/<eventType>/itc_overall/<subject>_itc_overall.mat` — **only if that file does not already exist**. This gives `overall.rayleigh_pvals` but leaves `overall.vs_rand_pvals` empty. If Pipeline B has already produced the file (with vs_rand results), it is not overwritten.

---

### Pipeline B — Cached ITC (3-step pipeline)

Same per-subject logic as Pipeline A, but separates wavelet computation for the following permutation tests.
The wavelet phases are saved to disk after Step 1, so Steps 2 and 3 can be re-run independently
without redoing the expensive convolution — useful when iterating on permutation settings or
trying different condition splits.

```
Step 1 ─ run_power_and_phases_itc.m   (once per subject, for every event type and the control)
            └─ data_path('power_and_phases_itc')power_and_phases_itc_<data_name>.mat

Step 2 ─ run_itc_vs_rand.m    (reads the cache for the real events and the randomised control, permutation test)
            └─ results/<eventType>/itc_overall/<subject>_itc_overall.mat

Step 3 ─ run_itc_by_condition.m   (reads the cache + overall result)
            └─ results/<eventType>/itc/<subject>_itc_<cond>_<split>.mat
```

#### What Step 1 caches (since 2026-09-24)

The Morlet coefficients of every channel and recording trial — `coef{ch}`, trials × 15 frequencies × samples, complex single, over [−0.75 6.85] s of each trial — so phase and power come from the same numbers, for any event type, without a per-event or per-dataset cache. Readers open it with `open_power_and_phases_itc(data_name, elec_name, event_type, dataset, time_win)`, which returns what the old cache held: `.time`, `.freqs`, `.n_channels`, and `.channel(ch)` → `.label`, `.valid_idx` (event filter `blinkSac == 0`, `firstInTrial ~= 1`, complete, and the channel's `noiseInfo.ok`), `.itc`, `.phases` (single, radians), `.power`. `dataset = 'ctrl'` reads the randomised-control events from the same cache.

The old per-event cache (`itc_phases`, retired) measured 270 GiB for the onset set alone, with the offset set and the control each needing as much again; the continuous cache is ~180 GB for everything.

Steps 2 and 3 are fast (no wavelet recomputation). Re-run Step 3 cheaply for any
`cond_var` / `split_type` combination without redoing the wavelet transform.

Step 2 adds `vs_rand_pvals` and `overall_rand`; if `itc_overall/` already exists it merges those fields and keeps the existing Rayleigh results. If the file is absent it creates it with both Rayleigh and vs_rand.
Step 3 requires `itc_overall/` to exist — it will not recompute overall ITC.

#### `itc_overall/` — canonical overall result

`results/<eventType>/itc_overall/<subject>_itc_overall.mat` is the single source of truth for overall ITC. Two pipelines write it:

| Source | `rayleigh_pvals` | `vs_rand_pvals` | `overall_rand` | Behaviour if file exists |
|--------|-----------------|-----------------|----------------|--------------------------|
| Pipeline A (`run_erp_per_subject.m`) | ✓ | empty | empty | skips (does not overwrite) |
| Pipeline B Step 2 (`run_itc_vs_rand.m`) | ✓ (if absent) | ✓ | ✓ | merges vs_rand only; keeps existing Rayleigh |

Step 2 writes the Rayleigh results itself when the file is absent, from the same trials as `vs_rand_pvals`, so running Step 2 is the way to populate `itc_overall/`.

`itc/` files contain **condition-split data only** — they do not embed the overall result.

---

### Pipeline C — Continuous regression (no condition split)

```
run_regression_per_subject.m
```

Instead of splitting trials into low/high groups, regresses the neural measure on `cond_var` as a **continuous predictor**, per timepoint, using all valid trials. Complements Pipeline A rather than replacing it.

- `measure_type`: `'erp'` / `'abs_erp'` / `'bha'` — **ITC is not supported**.
- `'abs_erp'` regresses on the **rectified single trials**, `|z|`, matching Pipeline D. ⚠️ Changed on 2026-07-28: it previously fitted the *signed* signal (the abs was applied only to the display mean), so its betas were identical to `'erp'`'s. New files carry `result.rectified = true`; **abs_erp regression files without that field hold signed-ERP betas and should be refitted.** Filenames did not change. See `notes/abs_erp_noise_confound.md` for why the rectified measure is confounded by single-trial noise in a way the signed one is not.
- Output: `results/<eventType>/<measure_type>/<elec_name>_<measure_type>_<cond_var>_regression.mat` — the same folder as the split results, distinguished by the `_regression` suffix (`_regression_rand` when `use_rand = 1`).
- Per channel, `.regression` holds `.r2`, `.beta`, `.pvals` (F-test), `.pvals_fdr`, and `.n`, each 1 × n\_time.
- `standardize_predictor` controls whether `.beta` is in predictor units or SDs.

---

### Pipeline D — Multiple regression (several predictors at once)

```
run_mreg.m  → mreg_fit(opts)
     └─ results/<eventType>/mreg/<subject>_<chan>_<measure>_<win_tag>__<pred_tag>.mat
     └─ figs/mreg/<eventType>/mreg__<measure>_<win_tag>__<pred_tag>__<subject>_<chan>.png
```

Pipeline C with all the predictors in **one** model, fitted point by point:

```
measure(t) ~ b0(t) + b1(t)·pred1 + b2(t)·pred2 + … + bk(t)·predk
```

**Main effects only** — no interaction terms. Each beta is that predictor's unique contribution with the others held constant, which is the difference from running Pipeline C once per predictor.

- `opts.predictors` — any list of the short names in `get_cond_field_map` (e.g. `{'amp','lum','lumAbsDiff','aws','awsAbsDiff','dg2cb','dg2cbAbsDiff','meaning','meaningAbsDiff'}`). Adding or removing one changes nothing else.
- `opts.transform` — optional struct, one field per predictor to transform **before** standardising: `'log'`, `'log1p'`, `'sqrt'`, `'none'`, or a numeric Box-Cox power λ (`(x^λ − 1)/λ`, log at λ = 0). Anything not named is left alone. This is where `inspect_predictors`' recommendations go, e.g. `struct('amp','log','lumAbsDiff',0.25,'awsAbsDiff',0.25,'meaningAbsDiff',0.25)`. A column holding exact zeros (the `*AbsDiff` family) is shifted by its smallest positive value first — printed, and stored as `.transform_shift`; a column that goes negative is an error rather than a silent shift. The transform is applied once per subject, to the raw column, after the missingness is settled and before any per-channel subsetting, so the shift is one constant for every channel of that subject. **It also lands in the filename** (`amp~log+lum+lumAbsDiff~p0.25`), so a transformed model cannot overwrite the untransformed one — and `run_plot_mreg` needs the same struct to find the file again. A transformed predictor's beta reads per SD of the *transformed* variable: for a log, per proportional change rather than per unit.
- `opts.standardize` — `'zscore'` (default) puts every predictor on a common scale, so the betas of different predictors are directly comparable; `'center'` and `'none'` are there for when they should not be.
- `opts.measure` — one of `'erp'` / `'abs_erp'` / `'bha'`; `opts.time_win` is the window fitted, `opts.baseline_win` the z-scoring baseline of the neural data. `'abs_erp'` fits `|z|` **per trial** — a response-magnitude measure, blind to polarity, whose betas are not sign-comparable to `'erp'`'s. Pipeline C rectifies the same way (since 2026-07-28), so the two agree.
- `opts.selection_mode` — `single` / `region` / `all`, resolved by the same `rsa_resolve_targets` the RSA pipeline uses.

**Trial set.** `rsa_event_filter`'s filter ∩ the channel's `noiseInfo.ok` ∩ the trials where **every** predictor is present. One trial set for the whole model, so the betas are comparable and the model is not refitted on shifting data.

**Saved per electrode** (`result`): `.beta .se .tvals .pvals .pvals_fdr` (n_pred × n_time, intercept excluded), `.intercept .intercept_p .intercept_p_fdr`, the measure being fitted (`.measure_mean .measure_sem`, the trial mean and its SEM per timepoint), whole-model `.r2 .r2_adj .model_F .model_pvals .model_pvals_fdr`, `.time .n .df .valid_idx`, and the provenance needed to redraw or audit the fit — `.predictors .pred_fields .standardize .pred_mean .pred_sd .pred_corr .vif`, plus `.transform .transform_shift .pred_labels` (what each predictor was transformed by; `.pred_labels` is what the figure writes on the rows, e.g. `amp [log]`).

`.measure_mean` is free: with centred predictors it is exactly the intercept. It is stored anyway so the figure stays correct under `standardize = 'none'`, where the intercept is not the mean. **Except for `abs_erp`**, where `.measure_mean` is the *rectified group mean* `|mean(signed ERP)|` with the SEM of the signed mean (`.measure_mean_mode = 'abs_of_mean'`, as in `run_erp_per_subject`) — the plotted row is then an evoked response rather than `mean(|ERP|)`, which is dominated by single-trial noise and sits well above zero everywhere. The betas underneath are still fitted on the rectified trials, so for `abs_erp` the mean row is *not* the intercept.

`.trial_noise .noise_corr .noise_corr_p` are the noise-confound check. `.trial_noise` is each trial's SD inside `baseline_win` — taken from the untransformed signal, before any evoked response, so a genuinely large response does not read as a noisy trial — and `.noise_corr` is each predictor's correlation with it (NaN if the proxy is unavailable or constant). A NOTE prints for `|r| >= opts.noise_corr_warn` (default `0.1`). It means different things per measure: on `erp`/`bha` a predictor that tracks trial noise only makes the residuals heteroscedastic (betas unbiased, classical SEs approximate), but on `abs_erp` it **biases the betas themselves**, since `E|y|` grows with the trial's noise SD. See `notes/abs_erp_noise_confound.md`.

`.vif` (variance inflation factor) and `.pred_corr` are the collinearity check: with correlated predictors the betas stay unbiased but get imprecise. `mreg_fit` prints a note for any VIF > 5, and `mreg_ols` errors outright on an exactly rank-deficient design instead of returning an arbitrary solution.

**FDR.** `opts.fdr_scope = 'predictor'` (default) runs BH over the timepoints of each predictor separately — one family of n_time tests per predictor, no correction across predictors; `'all'` treats the whole predictor × time matrix as a single family. Uncorrected p-values are always kept alongside, and the figure draws both levels.

---

## Stage 1b — FDR correction

Benjamini-Hochberg (BH) False Discovery Rate correction can be applied either **in-line** (automatically during Stage 1) or **post-hoc** via a standalone script on existing files. Both paths store corrected p-values as `_fdr` twin fields alongside the originals, which are never overwritten.

### In-line (automatic)

Each Stage 1 calculation script calls `apply_bh_fdr` immediately after computing p-values. No extra steps are required for newly processed data.

| Script | Fields with in-line FDR |
|---|---|
| `run_erp_per_subject.m` (erp/bha) | `ch.overall.pvals_fdr`, `ch.cond_pvals_fdr` |
| `run_erp_per_subject.m` (itc) | `ch.overall.rayleigh_pvals_fdr` |
| `run_itc_vs_rand.m` | `ch.overall.rayleigh_pvals_fdr`, `ch.overall.vs_rand_pvals_fdr` |
| `run_itc_by_condition.m` | `ch.cond_pvals_fdr` |
| `run_itc_band_by_condition.m` | `ch.band_cond_pvals_fdr` |

### Post-hoc (standalone script)

```
run_fdr_correction.m
```

Reads existing `.mat` result files and adds `_fdr` fields in-place. Use this to retroactively correct files produced before in-line FDR was added, or to re-correct with a different time window or frequency band.

**Key parameters:**

| Parameter | Meaning |
|---|---|
| `measure_type` | `'erp'` / `'abs_erp'` / `'bha'` / `'itc'` / `'itc_overall'` / `'itc_band'` |
| `fdr_win` | `[]` = all timepoints; `[t1 t2]` = restrict FDR family to this window (s) |
| `fdr_freq_band` | `[]` = all frequencies; `[f1 f2]` = restrict to band (Hz), ITC only |

**FDR scope:** per-channel. Each electrode's p-value set (1 × n\_time for ERP/BHA/ITC-band; n\_freqs × n\_time for ITC) forms one independent BH family.

The script records `result.fdr_params` (window, band, date applied) in each updated file for provenance.

### Utility function

```
apply_bh_fdr(p_in)
```

Accepts a p-value array of any shape; NaN entries are excluded from the BH family and remain NaN in the output. Returns BH-adjusted p-values (q-values) of the same shape — a test is significant at level α if `q < α`.

---

## Stage 2 — Visualise results

All visualisation scripts read pre-computed `.mat` files from Stage 1.
No raw data loading is required.

### `run_plot_erp_overall.m`

Plots the **overall** (no condition split) ERP / BHA / ITC for electrodes in a brain region.

- Significance markers: t-test vs zero (ERP/BHA) or Rayleigh test (ITC).
- ITC: two panels per electrode — ITC heatmap | Rayleigh significance mask (p < α).
- Output: `figs/lineplots_and_heatmaps/<eventType>/<measure_type>_overall/`

### `run_plot_erp_by_region.m`

Plots **condition-split** results (low vs high saliency/luminance/edge) per electrode, grouped by brain region.

- ITC can be shown as a time-frequency heatmap or as a frequency-band line plot (`itc_plot_mode`).
- Requires pre-computed files from Stage 1 with matching `cond_var` / `split_type`.
- Handles any number of conditions, so it is the driver for a `quantile<N>` split. `line_cmap` and `show_error` are auto by default: two conditions keep the distinct hues and the SEM shading, more than two switch to an ordered colormap (`parula`) and bare lines, since the colour then carries the bin order and overlapping bands would hide the curves. Set either explicitly to override (`line_cmap` takes any colormap function name). The legend sits outside the axes on the right (`legend_loc`, default `'eastoutside'`), narrowing the panel rather than covering the baseline — with five bins an inset legend hides the traces under it.
- Output: `figs/lineplots_and_heatmaps/<eventType>/<measure_type>_<cond_var>/`  (e.g. `bha_aws`, `itc_aws`; analysis type first since 2026-09-24)

### `run_plot_erp_by_region_multi.m`

Plots **multiple measures side-by-side** (e.g. ERP + |ERP| + BHA) for each electrode in a region.
Good for a quick overview without running separate `run_plot_erp_by_region` passes.

- Output: `figs/lineplots_and_heatmaps/<eventType>/multi_<cond_var>/`  (e.g. `multi_aws`)

### `run_plot_regression_by_region.m`

Plots pre-computed regression R² per electrode, grouped by brain region.

- Requires `_regression.mat` files from `run_regression_per_subject.m` (Pipeline C).

### `run_plot_mreg.m`

Redraws the multiple-regression figure from a saved `mreg` result — one row per predictor, its beta colour-coded over time (diverging blue-white-red, symmetric around 0), with a thin bar underneath marking the significant timepoints.

The bar shows **both** levels at once: light grey where the uncorrected p < α, dark grey where the FDR q < α. Since BH q ≥ p, the dark stretches are always nested inside the light ones, so "suggestive" and "survives correction" are readable in one pass instead of two renders.

Two optional line rows sit on top, in the same row format (line + its own significance bar), for context on what the model is fitting: `show_mean` draws the trial-mean of the measure ± SEM (significance = one-sample t vs zero), `show_r2` the whole-model R² (significance = the omnibus F). A line row is `line_row_scale` predictor rows tall and the figure grows to fit, so the vertical rhythm stays even.

- `run_mreg.m` already draws and saves this figure when `opts.preview_plot = true`; this driver exists to retune the figure (α, colour limits, row height, time range, context rows) without refitting.
- The settings that name the file — `measure`, `time_win`, `predictors` / `model_tag` — must match the analysis run; everything else is free.
- Both go through `mreg_plot_betas(result, popts)`, so the preview and the redraw cannot drift apart.
- Output: `figs/mreg/<eventType>/`

### `run_plot_mreg_summary.m`

Collapses a whole electrode set into **one figure per predictor**, where `run_plot_mreg.m` draws one electrode at a time (predictors down the rows). Thin driver over `mreg_summary(opts)`.

```
[ beta ]  electrodes (y) x time (x), one shared colour scale
[ sig  ]  white = n.s. | light grey = p < alpha | dark grey = FDR q < alpha
```

- **One analysis per summary** — one measure, one window, one predictor set — named exactly as `run_plot_mreg` names it (`bha_t-200_400ms__amp~log+lum+aws+meaning+lumAbsDiff~log`). Nothing is recomputed; every electrode shown must already have been fitted by `run_mreg.m`.
- Electrodes are grouped into contiguous per-subject blocks (`opts.sort_subjects` = `'map'`, the `get_subject_map` order, or `'name'`), and a colour bar beside **each** panel marks the subject in the study-wide `get_subject_colors` palette — the right-hand map is too far from the labels to track a row back across the colorbar.
- Row labels are one per electrode up to `max_row_labels` (default 40) and one per subject block above it, so a large region stays readable.
- `clim_mode = 'predictor'` (default) gives each figure one colour scale over its electrodes — the "shared colorbar for all channels"; `'all'` shares one scale across the predictors too, so the figures compare directly. `clim_pct` (default 99) sets the limit from a percentile of |beta| rather than its maximum, so a single extreme electrode cannot wash the map out; whatever that clips is counted on the console rather than left for the eye to find.
- Selection is the pipeline's usual `selection_mode` (`region` / `list` are the point of it); an empty one takes every electrode on disk with results for this analysis. As in `rsa_summary`, a globbed filename is verified against the provenance inside the file, so an `erp` run does not pick up `abs_erp` and `NS140` does not pick up `NS140_02`.
- Prints a per-predictor table (mean |beta|, % of electrode x time significant uncorrected and after FDR, how many electrodes have any FDR-significant timepoint) and returns the whole set as `.beta` / `.sig` (n_elec x n_time x n_pred) for further analysis.
- Output: `figs/mreg/<eventType>/mreg_summary__<measure>_<win_tag>__<pred_tag>__<selection>_<predictor>.png`

### `run_plot_rsa.m`

Redraws RSA figures from saved results — the RSA counterpart of `run_plot_mreg.m`. The three stages each preview while computing; this driver re-renders afterwards, for restyling (colour limits, colormap) or for previews a batch run skipped because of `max_preview`.

| `plot_stage` | Reads | Draws |
|---|---|---|
| `'combined'` | `results/<eventType>/rsa/` | `rsa_plot_preview` — brain matrix, stimulus matrix, second-order scatter |
| `'brain'` | `results/<eventType>/rsa_brain/` | `rsa_plot_rdm` — one brain matrix |
| `'stim'` | `results/<eventType>/rsa_stim/` | `rsa_plot_rdm` — one stimulus matrix |

- Nothing is recomputed: the Stage 3 file embeds both matrices along with `rho` / `p`.
- Selection works as in `run_rsa.m` (`single` / `region` / `all`); an empty `selection_mode` redraws every result on disk matching the measure / window / stimulus tags. As in `rsa_combine`, a globbed filename is verified against the provenance inside the file before drawing, so an `erp` run does not pick up `abs_erp` or `erp+bha` results.
- Both sides take **lists** exactly as `rsa_combine` does (`opts.measures` a cell of sets, `opts.stim_tags` a cellstr), so one run redraws every pairing; `max_plot` still caps the figures opened.
- Figures are written under the same `rsa_tags('fig', …)` names the stages use, so a replot overwrites the figure it came from.
- Output: `figs/rsa/<eventType>/`

### `run_plot_rsa_summary.m`

Collapses a whole electrode set into **one** figure, where `run_plot_rsa.m` draws one electrode at a time: one brain measure and one time window, against **several stimulus features on the x axis**.

```
y = mean Fisher z of rho ± SEM over electrodes      x = stimulus set (edg, lum, aws, meaning, …)
```

- Thin driver over `rsa_summary(opts)`, which reads `results/<eventType>/rsa/` — nothing is recomputed, so every pairing on the axis must already exist from Stage 3.
- **One** measure set and **one** window per summary (`opts.measures`, `opts.time_win`); more than one of either is an error rather than a silently mixed average. The stimulus sets are the list.
- Correlations are averaged as **Fisher z** (`z = atanh(rho)`), since rho is not additive; the mean, SEM and the one-sample t-test vs 0 (the stars) are all in z. `plot_units = 'r'` plots `tanh(mean z)` instead, with an asymmetric SEM bar — the transform is not linear. The second-order metric is Spearman by default, for which Fisher z is the usual approximation rather than exact. Electrodes are weighted equally regardless of trial count (`S.n_trials` keeps the counts if a weighted average is ever wanted).
- A stimulus set with **no** results at all is dropped from the axis with a note. An electrode missing from one set is, by default (`require_all_stims = true`), dropped from every column, so all columns describe the same electrodes and the means are joined by a line; set it false to average each set over whatever it has (per-column n then differs, and the line is dropped with it).
- Selection is the pipeline's usual `selection_mode` (`region` / `list` are the point of it); an empty one summarises every electrode on disk matching the tags.
- Returns the summary struct (`.z_mean .z_sem .r_mean .n_elec .p`, the per-electrode `.z .rho` and the electrode list) as well as drawing it.
- Output: `figs/rsa/<eventType>/rsa_summary__<sets>__<measure_tag>_<win_tag>__<selection>.png`

### `run_plot_split_vs_reg.m`

Plots the split-based result and the regression R² **side-by-side** for each electrode, one row per electrode: `[anatomy image] [split panel] [R² panel]`, a pair per measure.

- Requires both `run_erp_per_subject.m` (split) **and** `run_regression_per_subject.m` (regression) output for the same `cond_var`.
- The quickest way to see whether a condition split and a continuous fit tell the same story.
- Calls the same plotter as `run_plot_erp_by_region_multi.m` with `opts.show_regression = true`; that option is the only difference between the two drivers' output. Output files are prefixed `split_vs_reg_` instead of `multi_`.

### `run_plot_percentage_significant_all.m`

Compound summary figure: what fraction of responsive electrodes show a condition difference at each time point, with the selected electrodes plotted on brain surfaces.

- **Selection**: electrodes must pass the overall significance threshold (`alpha_selection`) in at least `select_pct_threshold`% of time points within `select_win_ms`.
- **Condition difference**: checked at each time point using `cond_pvals` vs `alpha_conditions`.
- Top strip: % selected electrodes significant per time point (shaded line plot).
- Lower panels: 4 brain surface views (left/right × lateral/medial), electrodes colour-coded by subject.
- Requires pre-computed files from Stage 1 (ERP/BHA only — uses `overall.pvals` and `cond_pvals`).
- `use_fdr = true` (default) — uses `_fdr` fields if present; falls back silently to uncorrected if absent.
- Output filename includes `_fdr` / `_uncorr` token to avoid overwriting.
- Output: `figs/percent_significant/<eventType>/<measure_type>_<cond_var>/`  (e.g. `percent_significant/offset/bha_aws`; analysis type first since 2026-09-24)

### `run_plot_percent_responsive.m`

Which channels respond to a visual event — saccade offset, image onset, or both — one figure per measure (ERP and BHA in one run).

- **Responsive** = at least `min_sig_pct` % of the timepoints in `resp_win_ms` (default 50–500 ms after the event, `min_sig_pct = 10`), and at least one, with BH-FDR q < `alpha_fdr`; the FDR family is that window's timepoints, per channel and event. `min_sig_pct = 0` is the plain any-timepoint rule, and the two do not overwrite each other — a non-zero value adds a `_min<N>pct` token to the file names. The test is Stage 1's `overall.pvals` (t-test vs 0 on baseline z-scored trials). The stored `pvals_fdr` are not used: their family is the whole epoch, which is longer for image (`[-0.5 1.5]` s) than for saccades (`[-0.5 0.5]` s), so the definition would differ between events.
- The window cannot end after 500 ms while the saccade results stop there; a window outside a result's time axis is an error, not a clip.
- Reads `results/offset/<measure>/*_<offset_result_suffix>.mat` (any split — `.overall` does not depend on it) and `results/image/<measure>/*_overall.mat`. Only channels with stats for both events are counted.
- `focus_modes` draws each figure twice in the same layout: `all` (every responsive channel, coloured by category) and `both` (only the channels responsive to both events, the others as grey context dots). The `both` figure is written with a `_both` suffix; the Venn and the summary still describe the whole classification.
- Figure: area-proportional Venn + summary (counts, %, median first significant time); glass brain (sagittal / coronal / axial, as in `run_plot_wang15_electrodes`) by category, non-responsive as grey dots, categories layered rarest on top; glass brain of the responsive channels by subject (`get_subject_colors`). Bipolar channels sit at the mean of their two contacts' `.FSAVERAGE` coordinates.
- Output: `figs/percent_responsive/responsive_<measure>_<win>_fdr[_min<N>pct][_both].png`, plus the per-channel table `T` (`resp_saccade`, `resp_image`, `category`, first significant times, % of the window significant, trial counts, coordinates) in `results/percent_responsive/responsive_<measure>_<win>_fdr[_min<N>pct].mat`.

### `run_plot_responsive_visual_overlap.m`

How the Wang 2015 visual channels overlap with the responsive ones, one figure per measure (`erp` / `bha`) x Wang labelling (`max_dist_mm` 3 / 6) x responsiveness set (`saccade` = responsive to saccade offset, `both` = to saccade offset and image onset), in the layout of `run_plot_percent_responsive`: Venn (visual / responsive / overlap) + summary; glass brains coloured visual only / responsive only / overlap, the rest as grey dots; glass brains of the `subject_row` channels (`overlap`, default, or `union`) by subject.

- **Visual** = a bipolar channel whose first contact is in any Wang area (all areas pooled), matched as in `run_plot_wang15_overall` (`LOi8` -> `LOi8-LOi9`).
- **Responsive** comes from the `run_plot_percent_responsive` table named by `resp_file_tag` (default `50to500ms_fdr_min10pct`), coordinates included; nothing is recomputed. The universe is every channel in that table.
- The summary gives the rates that make the Venn readable: % of visual channels responsive vs % of non-visual, and % of responsive channels that are visual. Channels cluster within subjects, so these are descriptive, not a test.
- Output: `figs/percent_responsive/visual_overlap/overlap_<measure>_<resp>_wang15_d<N>mm_<tag>.png` (one per comparison), `overlap_overview_<tag>.png` (all eight Venns on one page) and `overlap_counts_<tag>.csv`.

### `run_responsive_counts_table.m`

Per-subject channel counts for the same selections, one table per Wang labelling (`wang_dists_mm`, 3 and 6): `erp_saccade`, `erp_both`, `bha_saccade`, `bha_both`, `visual` (all Wang areas pooled) and one column per area group (`area_cols`, default V1 V2 V3 hV4 VO PHC LO TO V3AB IPS FEF), plus a `total` row. Rows follow `get_subject_map`.

- Every column counts **bipolar channels** over the `run_plot_percent_responsive` table (stats for both events), a channel being visual when its first contact is in a Wang area, so the totals are the Venn counts of `run_plot_responsive_visual_overlap`. Contact counts are in `results/wang15/wang15_summary_<tag>.xlsx` and differ, since the deepest contact of a shaft starts no channel.
- A visual channel in a group without a column (SPL1 by default) still counts toward `visual` and is reported on the console.
- Output: `results/percent_responsive/channel_counts_<resp_file_tag>_wang15_d<N>mm.csv`, and both tables as sheets of `channel_counts_<resp_file_tag>.xlsx`.

#### Shared by the two responsiveness drivers

| Helper | Purpose |
|---|---|
| `plot_venn_glass_figure(F, G, out_png)` | The page both drivers draw: two-set Venn + summary, glass brains by category (context channels grey, categories layered in `F.draw_order`), glass brains by subject. Everything it shows comes in the `F` struct (see its header). |
| `glass_brain_views(fs_dir)` | fsaverage silhouettes and the three projections (sagittal / coronal / axial) the glass brains are drawn on — the ones `run_plot_wang15_electrodes` uses. Load once per run. |
| `draw_venn(ax, n_regions, colors, set_names[, font_size])` | Area-proportional two-set Venn; a crescent too thin for its count gets the count outside it. |

### Brain surface maps — `run_plot_brain_map.m` / `run_plot_brain_map_animation.m`

Colours each electrode on an fsaverage surface by its measure value at one time point.

```
run_plot_brain_map.m            one map at opts.target_time
run_plot_brain_map_animation.m  one map per time point, compiled into an MP4
        └─ both call make_brain_map(opts)
```

`make_brain_map` is a **function** — every setting comes from the `opts` struct the driver builds, so the two drivers cannot drift out of sync. Add a new setting to `make_brain_map`'s header comment and to `getopt_field`, not to the drivers.

- `measure_type`: `'erp'` / `'abs_erp'` / `'bha'` / `'itc_band'` (the latter also needs `band_name`).
- `cond_mode`: `0` = average over conditions, `-1` = cond 2 − cond 1, `1` / `2` = that condition alone.
- `filter_to_selection = true` restricts to a pre-saved selection from `run_itc_band_select.m` or `run_plot_percentage_significant_all.m`; the selection file is derived from the other settings unless `opts.selection_file` is set.
- `opts.do_plot = false` extracts values to the `.txt` only and skips rendering (no iELVis needed).
- Output: `figs/brain_maps/<eventType>/<measure_type>/` — `values_<tag>.txt` (electrode, value, RGB) and `brain_<tag>.png`, one pair per time point. The animation adds a per-analysis subfolder with `frames/` and the `.mp4`.
- `opts.skip_existing = true` resumes an interrupted render; frame lookup and frame naming both go through `brain_map_tag`, so a resumed run finds the frames it already made.

**`opts.target_time` may be a vector.** The per-subject `.mat` files and the electrode coordinates are then read once for the whole series instead of once per frame — the dominant cost of an animation. For the 76-frame default this is ~4 s of extraction rather than ~270 s.

**Colour limits.** `opts.color_limits = []` derives symmetric limits from the data, computed once over every time point, so a series shares one scale. Still set them explicitly when animating: with `[]` plus `skip_existing`, a resumed render sees only the remaining time points and would rescale them against the frames already on disk (`make_brain_map` warns in that case).

---

## Visual-area atlas — Wang et al. 2015

Assigns each electrode contact to a retinotopic visual area of the Wang, Mruczek, Arcaro & Kastner (2015) maximum-probability atlas (25 areas: V1v/d … FEF), counts contacts and patients per area, and plots them. Independent of the EEG pipelines — it only needs the FreeSurfer recons and iELVis electrode files.

```
data/atlases/wang15/{lh,rh}.wang15_mplbl.v1_0.mgz   (fsaverage atlas, from neuropythy)
        │
run_wang15_atlas.m ── wang15_cortex() / wang15_volume()
        │   ──► results/wang15/wang15_contacts_<tag>.csv/.mat, wang15_summary_<tag>.xlsx,
        │       elec_locs_wang15_<tag>.xlsx, volumes/<sub>.wang15_mplbl.mgz
        ▼
run_plot_wang15_electrodes.m ──► figs/wang15/surface_<tag>.png, glass_<tag>.png,
                                  native/<sub>/<sub>_<elec>_<tag>.png
```

- **Labelling is done in native space.** The atlas is carried to each subject through `sphere.reg` (nearest vertex, as `mri_surf2surf` does), each vertex is sampled at several depths between the white and pial surfaces, and a contact (iELVis `.PIAL` coordinates) takes the label of the nearest sample. A contact is assigned only if that nearest grey matter is within `max_dist_mm` (default 3) **and** carries a Wang label. `.FSAVERAGE` coordinates are used only for drawing — for depth electrodes they come from an affine transform.
- `tag = d<max_dist_mm>mm`, so different thresholds do not overwrite each other. The contacts table keeps `dist_cortex_mm`, `nearest_label`, `nearest_wang_area`, `dist_nearest_wang_mm` and `aparc_aseg` for every contact, so another criterion can be applied without rerunning.
- **Counts**: "datasets" are rows of `get_subject_map()`; "patients" merge datasets of one person by their `NSxxx` number (NS140 + NS140_02, NS174_02 + NS174_03). Contacts are counted per dataset.
- `elec_locs_wang15_<tag>.xlsx` has the `elec_locs_clean.xlsx` layout (`sub` + one column per area group, plus `early_visual` and `any_wang`), so it can be passed as `elec_locs_path` to the by-region plotters. It lists contacts, not bipolar pairs.
- `wang15_areas()` is the one place that decodes label numbers into names, area groups (V1 = V1v+V1d, IPS = IPS0–5, …) and colours.
- `run_plot_wang15_overall.m` plots the ERP / BHA of those contacts, one small panel per contact and one page per area, from any per-subject Stage 1 file (`result_suffix`, default `amp_quantile5`, computed from the pad epochs for every subject with visual contacts). `plot_split = false` draws the overall response (mean ± SEM, all trials, bar where p < `alpha_zero` vs 0); `true` draws one line per condition of the file's split (parula, no SEM). Contacts map to the bipolar channel they start (`LOi8` → `LOi8-LOi9`), so the last contact of a shaft has none. `event_type = 'image'` plots the image-locked results (`result_suffix = 'overall'`, `plot_split = false`, since those files have no split). Output: `figs/wang15/overall/<event_type>/wang15_<overall|result_suffix>_<measure>_<group_by>_<tag>.pdf`.
- The volumes in `results/wang15/volumes/` overlay the subject's `mri/T1.mgz` in freeview (colour by lookup table or heat).
- `knnsearch` is called through `KDTreeSearcher`: FieldTrip's `external/stats/knnsearch.m` shadows the Statistics Toolbox function and is brute force.

---

## ITC frequency-band pipeline

Band-specific ITC analysis: selects electrodes responsive in a given frequency band (e.g. alpha 8–13 Hz), runs a condition permutation test on band-averaged ITC, and produces the same compound summary figure as above.

```
Step 1 ─ run_itc_band_select.m
            Reads  results/<eventType>/itc_overall/<elec>_itc_overall.mat
                   (overall.rayleigh_pvals — condition-independent)
            Writes results/<eventType>/itc_band/<band>_<eventType>_selection.mat

Step 2 ─ run_itc_band_by_condition.m
            Reads the continuous wavelet cache (open_power_and_phases_itc), slices to band freqs, permutes
            Writes results/<eventType>/itc_band/<elec>_itc_band_<band>_<cond>_<split>.mat

Step 3 ─ run_plot_percentage_significant_itc_band.m
            Loads selection + band condition results, produces compound figure
            Writes figs/percent_significant/<eventType>/itc_band_<cond_var>/
```

**Key parameters** (set at the top of each script):

| Parameter | Meaning |
|-----------|---------|
| `itc_freq_band` | `[f_low f_high]` Hz — e.g. `[8 13]` alpha, `[13 30]` beta, `[4 8]` theta |
| `select_win_ms` | Time window for electrode selection |
| `alpha_selection` | Rayleigh test threshold for selection |
| `select_pct_threshold` | % of time points in window that must be significant |
| `select_basis` | `'rayleigh_fdr'` (default) or `'rayleigh'`; `'vs_rand'` / `'vs_rand_fdr'` are stubs, not yet implemented |
| `alpha_conditions` | p-value threshold for condition difference in the figure |
| `use_fdr` | `true` (default) — use `band_cond_pvals_fdr`; `false` — use uncorrected `band_cond_pvals` |

**Notes:**
- Step 1 reads from `itc_overall/` (condition-independent); `cond_var` / `split_type` are not used for the source file.
- Step 2 reads the continuous wavelet cache (`run_power_and_phases_itc.m`, through `open_power_and_phases_itc`); a subject without it errors. Its frequencies and cycles are the cache's; `time_win` in Step 2 sets the window read.

---

## RSA pipeline (representational similarity analysis)

For one electrode, RSA builds two **trial-by-trial** similarity matrices — a **brain** matrix from neural time-series features and a **stimulus** matrix from image features around the gaze point — then correlates the two matrices (the **second-order** correlation). It runs for a single preview electrode, electrodes from a region file, or all electrodes, using the same selection as the by-region plotters.

A 3-stage decoupled pipeline, modelled on the cached-ITC pipeline. Each stage is an `opts`-driven engine function (`rsa_brain` / `rsa_stim` / `rsa_combine`) with a thin `run_*` driver; `run_rsa.m` sets one shared `opts` and runs all three in order.

```
Stage 1 ─ run_rsa_brain.m  → rsa_brain(opts)
            └─ results/<eventType>/rsa_brain/<elec>_<chan>_<measure_tag>_<win_tag>.mat

Stage 2 ─ run_rsa_stim.m   → rsa_stim(opts)
            └─ results/<eventType>/rsa_stim/<elec>_<chan>_<stim_tag>.mat

Stage 3 ─ run_rsa_combine.m → rsa_combine(opts)
            └─ results/<eventType>/rsa/<elec>_<chan>_<measure_tag>_<win_tag>__<stim_tag>.mat

run_rsa.m  → runs all three with one shared settings block.

Previews (any stage) → figs/rsa/<eventType>/*.png
```

### Brain side (Stage 1)

- `opts.measures` — any of `'erp'` / `'abs_erp'` / `'bha'`, **concatenated** in order into the per-trial feature vector. ERP from `get_epochs` (baseline z-scored, per-channel `noiseInfo.ok` filtered); BHA from the continuous `power_bha` cache through `get_bha_channel` (over the old cache's span, [-0.5 0.5] s ± 0.05).
- A cell of **sets** (`{{'erp'},{'erp','bha'}}`) builds one brain matrix per set in the same pass over the data; a plain cellstr stays one set. `rsa_specs` makes that reading, so every stage agrees on what a given `opts.measures` asks for.
- `opts.time_windows` — cell of `[t1 t2]` (s); Stage 1 **loops** them, one brain matrix per window.
- First-order matrix = `corr(features', opts.first_order_metric)` across trials (default `'Spearman'`).

### Stimulus side (Stage 2)

- `opts.stim_fun` — a pluggable function `feat = f(content, valid_idx, sopts)` returning `n_valid × feat`. Default `@rsa_stim_patch`: a round gaze patch (radius `opts.radius_px`) around `(Xpx_square, Ypx_square)`, vectorised. `opts.stimuli_dir` chooses the image set (`stimuli_square`, `stimuli_square_aws`, …).
- The stimulus feature is per-epoch, so there is **one** stimulus matrix per electrode, independent of the brain time windows.
- **One image set per run** — Stage 2 has to build it, so rerun the driver per set. `opts.stim_tags` is a Stage 3 setting and is rejected here rather than quietly building something else.

### Second order (Stage 3)

- Works off the Stage 1/2 files (no raw data): pairs each brain matrix with its electrode's stimulus matrix, **asserts the trial sets match**, and correlates the strict upper triangles with `opts.second_order_metric` (default `'Spearman'`).
- **Any brain matrix pairs with any stimulus matrix, and both sides take lists** — the run is their cross-product, which is what the split into stages buys: a matrix shared by several comparisons is built once and read once here.

  ```matlab
  opts.measures     = {{'erp'}, {'bha'}};                    % brain: two sets
  opts.time_windows = {[-0.2 0.4]};
  opts.stim_tags    = {'patch_r50_aws', 'patch_r50_edg'};    % stimulus: two sets
  %  -> 4 pairings per electrode, 4 output files
  ```

- Stimulus matrices are named by the **tag in their filename** (`opts.stim_tags`), which is all this stage needs — it opens no images. Leaving `stim_tags` out falls back on the Stage 2 recipe (`stim_fun` / `stimuli_dir` / `radius_px` / `grayscale`) as the single stimulus set, but re-deriving a name from a recipe is how a stale `stimuli_dir` in the Stage 3 driver silently pairs a brain matrix with the *wrong* stimulus set.
- A requested set that was never computed is an **error naming the tag and the folder**, raised before any pairing runs; an electrode merely missing from an existing set is a warning and a count in the summary.
- It honours the same `opts.selection_mode` as Stages 1–2, so a `single` run combines only that electrode. Without this it would recombine — and re-preview — every electrode left in `rsa_brain/` by earlier runs, which with a small `max_preview` means the preview you see may not even be the electrode you selected. **Omit** `selection_mode` to deliberately sweep every brain file matching the tags (useful after a batch run).

### Trial alignment

The per-electrode trial set is `valid_idx = find(event_filter & per_ch_ok(:,ch))`, where `event_filter` (from `rsa_event_filter`) is the `run_erp_per_subject` filter (`blinkSac==0 & firstInTrial~=1 & is_complete==1`) **plus** a finite-gaze requirement so every selected trial has a stimulus feature. Both Stage 1 and Stage 2 recompute this identically from the same events (`get_epochs` / `get_events` + `content`), so they align without depending on one another; Stage 3 asserts equality before correlating. Edge-clipped gaze patches are replicate-padded so no trial is dropped.

### Preview plotting

Set `opts.preview_plot = true` (meant for `selection_mode='single'`) to visualise the result. Each stage previews its own output — Stage 1 the brain matrix, Stage 2 the stimulus matrix, Stage 3 the combined figure (brain matrix + stimulus matrix + second-order scatter with rho / p). `run_rsa.m` draws only the Stage 3 combined figure (it suppresses the per-stage matrix figures to avoid duplicates). `opts.max_preview` caps how many figures a run may open, so previewing during a `region` / `all` run stays bounded.

Every drawn preview is also **saved** (`rsa_save_fig`) to `figs/rsa/<eventType>/`, as a 150-dpi PNG named by `rsa_tags('fig', …)`:

| Stage | Figure file |
|---|---|
| 1 (brain) | `rsa_brain__<measure_tag>_<win_tag>__<subject>_<chan>.png` |
| 2 (stim) | `rsa_stim__<stim_tag>__<subject>_<chan>.png` |
| 3 (combined) | `rsa__<stim_tag>__<measure_tag>_<win_tag>__<subject>_<chan>.png` |

e.g. `rsa__patch_r50_aws__erp_t-200_400ms__NS128_02_LDa1.png`. `run_plot_rsa.m` redraws any of these three from the saved `.mat` (see below). Semantic blocks are joined with `__` because subject names (`NS128_02`) and window tags (`t-200_400ms`) contain single underscores; putting the analysis settings first and the subject/electrode last makes runs with the same settings sort together. Saving is on by default — set `opts.save_preview = false` for on-screen only, `opts.save_preview_fig = true` to also write a MATLAB `.fig` next to each PNG, or `opts.figs_dir` to redirect the output folder. Since only drawn previews are saved, `max_preview` also caps how many PNGs a run writes.

### RSA helpers

| Helper | Purpose |
|---|---|
| `rsa_resolve_targets(opts)` | Turn `selection_mode` (`single` / `region` / `all`) into a `{data_name, elec_name, chan_prefix}` list; bridges elec_locs names to data names via `get_subject_map`. |
| `rsa_event_filter(data_perisac, content)` | Canonical RSA trial filter (`event_filter`, `per_ch_ok`) — the single source of truth shared by Stages 1 and 2. |
| `rsa_channel_idx(labels, chan_prefix)` | Resolve a channel prefix (or `'*all*'`) to indices into a label list. |
| `build_brain_feature_mat(...)` | Per-trial neural feature matrix for one channel: z-score + crop each measure, concatenate. |
| `rsa_first_order_rdm(F, metric)` | Trial-by-trial similarity matrix from a feature matrix. |
| `rsa_second_order(M1, M2, metric)` | Correlate the strict upper triangles of two matrices; returns `[rho, p]`. |
| `rsa_stim_patch(content, valid_idx, sopts)` | Default stimulus feature: vectorised round gaze patch. |
| `rsa_specs(kind, opts)` | Normalise `opts.measures` / `opts.stim_tags` into explicit **lists** of brain and stimulus specs, so no two stages can read the same opts differently. |
| `rsa_summary(opts)` | Group summary over the Stage 3 files: one brain matrix vs several stimulus sets, averaged across electrodes in Fisher z. Driven by `run_plot_rsa_summary.m` (see Stage 2). |
| `rsa_tags(kind, arg)` | Canonical filename tags (`measure` / `win` / `stim` / `stimfun` / `stimset` / `fig`) shared by all three stages. `stimset` marks the screen-space sets (`stimuli_aws` → `screen_aws`) so they cannot collide with the 800×800 square sets (`stimuli_square_aws` → `aws`). |
| `rsa_plot_rdm(rdm, ttl[, ax][, popts])` | Show one similarity matrix; `popts.clim` / `.cmap` override the default symmetric diverging scale. Returns the figure handle. |
| `rsa_plot_second_order(M1, M2, metric, rho, p[, ax])` | Scatter the two matrices' upper triangles with a fit line. |
| `rsa_plot_preview(brain, stim, rho, p, metric[, popts])` | Combined single-electrode preview figure (both matrices + scatter); `popts` passes colour settings to both matrix panels (`.clim`, `.clim_brain` / `.clim_stim`, `.clim_mode = 'shared'`, `.cmap`). |
| `rsa_save_fig(hfig, opts, fig_name)` | Write a preview to `figs/rsa/<eventType>/` (PNG at 150 dpi, optional `.fig`). |

---

## Multiple-regression helpers

| Helper | Purpose |
|---|---|
| `mreg_fit(opts)` | The engine: selection → trial set → design matrix → per-timepoint fit → save (+ optional preview). |
| `mreg_ols(Y, X)` | The statistics, independent of iEEG specifics: one QR of `X`, then betas, SEs, t, p, R², adjusted R² and the omnibus F for every column of `Y`. Verified against `fitlm`. |
| `mreg_predictor_matrix(content, predictors[, n])` | The raw predictor columns as an n_trials × n_pred matrix (content may be a table or a struct). |
| `mreg_transform(P, predictors, transform)` | Applies `opts.transform` to the raw columns: `log` / `log1p` / `sqrt` / `none` / a numeric Box-Cox power, shifting a zero-containing column off zero and refusing to shift a signed one. Returns the transformed matrix plus the specs, shifts, filename tags and figure labels actually used. |
| `mreg_plot_betas(result[, popts])` | The figure: one beta strip per predictor plus its significance bar (one electrode). |
| `mreg_summary(opts)` | Group summary over one analysis: per predictor, a beta heatmap and a significance heatmap of electrodes x time, with a subject colour bar. Driven by `run_plot_mreg_summary.m` (see Stage 2). |
| `mreg_tags(kind, arg[, transform])` | Canonical `tf` / `pred` / `win` / `file` / `fig` names, so the analysis and plotting drivers agree on filenames — including the transform, which is part of the `pred` tag. |

`mreg_fit` reuses the RSA pipeline's `rsa_resolve_targets`, `rsa_event_filter`, `rsa_channel_idx`, `build_brain_feature_mat` and `rsa_save_fig` (with `opts.figs_dir` pointed at `figs/mreg/`) — the electrode selection, trial filter and neural-feature extraction are deliberately the same code as RSA's brain stage.

---

## Standalone scripts — `mfiles/standalone/`

Run directly and by hand; nothing in the pipelines calls them. They are kept out of `mfiles/` so that folder holds only pipeline code.

**Path note:** these scripts call helpers that live in `mfiles/` (e.g. `inspect_channel` → `cleanPerisaccadicEEG`), and MATLAB does not add subfolders to the path automatically. To run one, either `cd` into `mfiles/standalone/` with `mfiles/` still on the path, or add both folders.

### Data generation / preparation

| Script | Purpose |
|---|---|
| `make_aws_maps.m` | Generates the AWS saliency maps in `data/stimuli_square_aws/`. |
| `make_stimuli_order.m` | For each subject and trial (1–80), identifies which stimulus was shown. |
| `generate_fsaverage_coords.m` | Produces the fsaverage electrode coordinates used by the brain-surface plots. |
| `rebuild_content.m` | Regenerates `data/sacextr_ekm` end to end: `make_contentEKM` (all subjects, overwrites), `concatenateAllSubjects`, `add_missing_imstats`, `add_saccade_sequence_to_content`, `add_return_saccades`, and `export_return_saccades_csv` once the ResNet columns exist. Does **not** rebuild `resnet50_*_dist` — rerun `train_model_r50px.ipynb` on the new content, then `add_patch_repr_dist.m`. |
| `add_missing_imstats.m` | Fills image-statistic columns a content file never got, so a subject stops being silently dropped from every model that uses one (`subNS128_02` had no `aws_*` / `lum_*`). Compares each `content_EKM_*.mat` against the target column list, computes only what is missing with the same `add_imstat_to_content` call `run_preprocess` uses, never touches an existing column, backs the file up first, and cross-checks the result against `content_all_subjects.mat`. |
| `add_patch_repr_dist.m` | Appends the `*_dist` columns of `od_michala/applying_models/patch_r50px/patch_repr_r50px_dist.csv` (from `train_model_r50px.ipynb`: cosine distance between a network's representation of the r = 50 px gaze patch and the previous fixation's) to every `content_EKM_*.mat` and `content_all_subjects.mat`. Rows are matched on (`name`, `trial_number`, `indx_sac`) and `Xpx_square` is cross-checked; existing columns are kept unless `overwrite = true`; files are backed up to `data/sacextr_ekm/backup_before_patch_dist/`. Rerun after adding a model to the notebook. |
| `add_return_saccades.m` | Marks return saccades and refixations, and appends eight 0/1 columns (`returned_to`, `return_1back`, `return_2back`, `return_3back`, `return_any`, `refixated`, `refixation`, `departure` — the n+1 of a 1-back return, the fixation that took the gaze away from the anchor) plus `amplitude_px` (euclidean step from the previous fixation in `Xpx_square` / `Ypx_square`) and `amplitude_px2deg` (the same over `px_per_deg`, the basis of the notebook's amplitude analyses) to every `content_EKM_*.mat` and `content_all_subjects.mat`. The gaze counts as having *left* a spot once it is >= `thresh_px` from it (default 78 px = 2 dva in `Xpx_square` space, see `notes/visual_angle_and_patch_sizes.md`). A fixation is a **1-back** return when the one fixation in between was away from the anchor two back (n+2 -> n), a **2-back** return when *both* fixations in between were away from the anchor three back (n+3 -> n), and a **3-back** return when *all three* were away from the anchor four back (n+4 -> n). `return_any` is the union of the three depths, the counterpart of `returned_to`; it is slightly smaller than the anchor count, since one fixation can return to two anchors at once (122 of 3,051 returns do). Every fixation between the anchor and the return has to be outside the circle, so the gaze demonstrably left rather than lingered; requiring it of n+1 as well means a spot that was refixated anchors the return to the **last** fixation of that run, not an earlier one. `returned_to` marks the anchor of a return at **any** of the three depths (3,175 fixations; 1,563 are returned to 1-back, 927 2-back, 685 3-back), `refixation` a fixation within the threshold of its predecessor, and `refixated` the start of a spot that the next fixation refixates. Trial boundaries are respected and missing gaze is NaN, not 0 -- both the fixation itself and any verdict that depends on it, unless the verdict is settled regardless. `opts.skip_invalid = true` steps over gazeless fixations instead of letting them block a return (it moves the rates by <1.5 pp). Existing columns are kept unless `overwrite = true`; files are backed up to `data/sacextr_ekm/backup_before_return_saccades/`. The work is done by the helper `add_return_saccades_to_content(content[, opts])`, which takes any fixation table. |
| `preview_return_saccades.m` | **Shows the 1-back and 2-back returns only** — it was not extended when `return_3back` / `return_any` were added, so a 3-back return draws no arrow and panel C has five rows, not seven. Draws one trial's scanpath and shows why each fixation got the flags `add_return_saccades.m` gave it: the returns with a dashed threshold circle round every `returned_to` anchor and a curved arrow from each return back to it, the refixation runs, and the five columns as a grid along `indx_sac` over the distance from the previous fixation. A return lands within one threshold of its anchor by definition, so their markers overlap; anchors are drawn larger and behind, so the pair reads as a small dot inside a bigger one, and labels are left on their own marker — this is a preview, and panel C is the record. Set `data_name` / `trial` at the top; writes `figs/return_saccades/<subject>_trial<N>.png`. |
| `export_return_saccades_csv.m` | Dumps the columns `notebooks/explore_return_saccades.ipynb` needs from `content_all_subjects.mat` to `data/return_saccades/return_saccades.csv` (regenerable, gitignored) — pandas cannot read a MATLAB table out of a .mat. Rerun after `add_return_saccades.m`. The notebook covers how often each event happens (per subject, mean ± SEM), how the five flags overlap, and the amplitude / latency / dwell distributions against a no-flag baseline. It also writes `params.json` beside the CSV with the threshold the columns were actually built with, so the notebook reads that instead of hardcoding a number that could drift. **Read its section 3 before trusting any amplitude contrast**: a saccade under the threshold *is* a refixation, so `refixation` sits entirely below 78 px while `refixated` and the no-flag baseline sit entirely above it — all three are truncated by the definition. Restricted to the >= 78 px stratum the 1-back returns and the anchors keep a genuinely shorter amplitude (4.2 and 4.5 deg vs 4.9, p < 0.01) while 2-back and `refixated` lose theirs. The notebook also checks eccentricity per class (every class sits nearer the image centre than baseline — the main confound), the image statistics at the fixated point, the QC flags per class, whether a short first look predicts a return, and a threshold sweep from 0.5 to 3 dva. The sweep reimplements the rule in Python and **asserts it reproduces the stored MATLAB columns exactly** before plotting anything. |

### Inspection / QA

| Script | Purpose |
|---|---|
| `inspect_channel.m` | One channel's cleaning (version 2): kept epochs next to the epochs each step rejected, each step's fence against its per-second maxima, IED and flag summary. |
| `inspect_epoch.m` | Single epochs of one channel with the samples each step marked (version 2). |
| `audit_cleaning.m` | Audit of the version-1 cleaning against Michal's (`audit_cleaning_subject`); figures from `plot_cleaning_audit.m`. |
| `diag_cleaning_thresholds.m`, `diag_threshold_methods.m` | Drift, tails of the per-epoch statistics, butterflies either side of a threshold, IED waveforms; null / excess / breakpoint threshold candidates. |
| `diag_ied_threshold.m` | IED threshold from a phase-randomised null (FDR), waveforms by score, detection rate around saccades. |
| `plot_cleaning_v2_effects.m` | Version 2 vs the stored version-1 `noiseInfo` on a few recordings, in memory: fences, newly rejected / kept epochs, bad channels, ERPs, IED locking, saccadic spike; table `results/cleaning_audit/v2_effects.csv`. |
| `preview_masks.m` | Checks that the circular masks used in `add_imstat_to_content` land where intended. Calls `preview_fixations`. |
| `preview_fixations.m` | Draws fixations over a stimulus; used by `preview_masks` and runnable on its own. |
| `preview_image.m` | Quick look at a stimulus and its feature maps. |
| `preview.m` | Scratch plotting sandbox (saccade-amplitude histograms and similar). |
| `demo_gaze_angle_to_center.m` | Validates `gaze_angle_to_center` against `notes/Tobii_coordinates.xlsx`. |
| `plot_corr_amp_awsdiff.m` | Scatter + linear fit: saccade amplitude vs AWS difference. |
| `inspect_predictors.m` | Pre-flight for Pipeline D: takes an `opts.predictors` set and the same electrode selection as `run_mreg`, and reports what the design looks like before a model is fitted — Pearson / Spearman / **partial** correlation matrices, VIF and the eigenvalues of the correlation matrix, per-variable skewness with the best Box-Cox power, the trimmed and log-log \|r\| of every pair, the trials each predictor costs, and how stable each correlation is across subjects. It then applies the recommended transforms and measures the design again (skew, VIF, r, partial r before vs after). Loads content files only (no epochs). Four figures to `figs/mreg/predictors/`, and the whole report to `results/mreg_predictors/predreport__<predictors>.txt`. |
| `explore.m` | Scratch pad for poking at a single subject's raw `.mat`. |

### Environment

| Script | Purpose |
|---|---|
| `recompile_fieldtrip_mex.m` | Rebuilds FieldTrip MEX files for the current MATLAB version on Windows. Run after a MATLAB upgrade. |

---

## Key parameters (shared across scripts)

| Parameter | Values | Meaning |
|-----------|--------|---------|
| `eventType` | `'onset'` / `'offset'` / `'image'` | Epoch aligned to saccade start or end, or to image onset (one epoch per trial, no condition split; see `run_image_epochs.m`) |
| `measure_type` | `'erp'` / `'abs_erp'` / `'bha'` / `'itc'` | Neural measure to compute or plot |
| `cond_var` | `'aws'` / `'lum'` / `'edg'` / `'dg2cb'` / `'meaning'` / `'amp'` / `'vel'` / `'lat'`; next-saccade variants `'ampNext'` / `'velNext'` / `'latNext'`; difference variants `'awsAbsDiff'` / `'lumAbsDiff'` / `'dg2cbAbsDiff'` / `'meaningAbsDiff'`; network patch distances `'resnetL4Dist'` / `'resnetL3Dist'` (cosine distance between the ResNet50 layer4 / layer3 representations of the r = 50 px patch and the previous fixation's); with `split_type = 'flags'`, 0/1 return-saccade columns joined by `+` (`'returned_to+return_1back'`) | Image property or saccade metric used for the condition split, or the continuous predictor in the regression pipelines. `get_cond_field_map` defines the full list for every script; a name is usable wherever the content file actually has that column (`edg` is missing from most). |
| `predictors` | cellstr of `cond_var` names | The predictors of one multiple-regression model (Pipeline D). All enter the same fit, as main effects only. |
| `split_type` | `'tertile'` / `'median'` / `'quantile<N>'` / `'flags'` | How trials are divided; `flags` makes one condition per listed 0/1 column (Pipeline A only). `tertile` and `median` give the low / high pair every downstream script expects; `quantile5` gives 5 equal-count bins, for plotting a graded set of curves (Stage 1 + `run_plot_erp_by_region.m` only). |
| `use_rand` | `0` / `1` | Real data vs randomised-event-times control |
| `region_col` | e.g. `'antHP'`, `{'ITC','antHP'}` | Brain region(s) to plot |

---

## Full pipeline at a glance

```
run_preprocess.m  ──► data_path('reref')  (bipolar, ±5 Hz notch, demeaned)     random_sac.m ──► data/sacextr_rand/
make_contentEKM / rebuild_content ──► data/sacextr_ekm/                        (read by get_events, use_rand = 1)
      │
      ├── get_events / get_epochs   ── events + epochs cut on demand, flags from get_artifact_marks
      │                                (data_path('artifact_marks'), built on first use)
      ├── run_power_bha.m            ──► data_path('power_bha')              [continuous, all event types]
      ├── run_power_and_phases_itc.m ──► data_path('power_and_phases_itc')   [continuous, all event types]
      │
      ├──── Pipeline A ── run_erp_per_subject.m      (condition split)
      │        measure_type: erp / abs_erp / bha ──► results/.../erp|abs_erp|bha/
      │        measure_type: itc ─────────────────► results/.../itc/          (conditions only)
      │                                             results/.../itc_overall/  (rayleigh; only if absent)
      │
      ├──── Pipeline B ── cached ITC
      │        (run_power_and_phases_itc.m, above)
      │        run_itc_vs_rand.m                   ──► results/.../itc_overall/  (adds vs_rand; keeps rayleigh)
      │        run_itc_by_condition.m              ──► results/.../itc/          (conditions only)
      │
      ├──── Pipeline C ── run_regression_per_subject.m   (one continuous predictor; no split)
      │        measure_type: erp / abs_erp / bha ──► results/.../<measure>/*_<cond>_regression.mat
      │
      └──── Pipeline D ── run_mreg.m ── mreg_fit()       (several predictors in one model)
               measure: erp / abs_erp / bha ─────────► results/.../mreg/*.mat

                                    ▼  Stage 2 — visualise

 results/.../<measure>/  ──► run_plot_erp_overall.m           ──► figs/lineplots_and_heatmaps/…/<measure>_overall/
                        ──► run_plot_erp_by_region.m          ──► figs/lineplots_and_heatmaps/…/<measure>_<cond>/
                        ──► run_plot_erp_by_region_multi.m    ──► figs/lineplots_and_heatmaps/…/multi_<cond>/
                        ──► run_plot_percentage_significant_all.m
                                                         ──► figs/percent_significant/…/<measure>_<cond>/
                        ──► run_plot_brain_map.m / run_plot_brain_map_animation.m
                                 └─ make_brain_map()     ──► figs/brain_maps/…/<measure>/  (.txt + .png [+ .mp4])

 *_regression.mat       ──► run_plot_regression_by_region.m       (R² only)
                        ──► run_plot_split_vs_reg.m   (split + R² side-by-side)

 results/.../mreg/      ──► run_plot_mreg.m ── mreg_plot_betas() ──► figs/mreg/…/  (one electrode)
                        ──► run_plot_mreg_summary.m ── mreg_summary()
                                                        ──► figs/mreg/…/  (electrodes x time, per predictor)

 itc_overall/           ──► run_itc_band_select.m  ──► itc_band/<band>_selection.mat
                                 └─► run_itc_band_by_condition.m
                                         └─► run_plot_percentage_significant_itc_band.m
                                                 └─► figs/percent_significant/…/itc_band_<cond>/
```
