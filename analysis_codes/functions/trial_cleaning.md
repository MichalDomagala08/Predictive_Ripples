# Trial cleaning: Michal's procedure, version 1, version 2

*2026-09-23. Written from the code, not from its comments: where the two disagree, the code is what is described. Evidence: `mfiles/standalone/audit_cleaning.m`, `diag_cleaning_thresholds.m`, `diag_threshold_methods.m`, `check_cleaning_v2.m`; figures in `figs/cleaning_audit/`.*

Both our versions share the upstream steps: bipolar re-referencing within electrode shafts (`preprocess_bipolar`, from `data/raw`), 500 Hz, and line-noise removal in `cleanTrialEEG`. That is a 4th-order Butterworth band-stop at **60** Hz and 3 harmonics (120, 180, 240), zero-phase, per 14 s trial: **±2 Hz** in version 1, **±5 Hz** from 2026-09-24 (see *Tests* below), with the filter now inside `cleanTrialEEG` rather than Michal's `BstFilter`. 60 Hz is right for these recordings (US); the spectra show no 50 Hz peak. Rebuilding a reref file from raw with the ±2 Hz setting reproduces the existing one to 10⁻¹³ of its range (NS144_02).

## 1. Michal's procedure (`od_michala/cleaning_trials`)

Written for continuous hippocampal recordings and ripple detection; artifacts are marked as **samples** and set to NaN, and no epochs are rejected.

**Preprocessing** (`MW_0_Data_Loading.m`). Bipolar re-reference (FieldTrip, within groups) → resample to 500 Hz (`ft_resampledata`) → band-stop at **50**, 100, 150, 200 Hz, 4th-order Butterworth ±2 Hz, zero-phase (`BstFilter`). The paper's methods text gives the order resample → notch → re-reference.

**Artifacts** (`Artifacts_VAZ.m`), per channel over the whole recording, no demeaning:
- *IQR step* — non-overlapping segments of 251 samples starting every 250 (~0.5 s); a segment is marked whole when `|mean| >= |Q3| + 2.3·IQR` of the channel (after Vaz et al.). The test is on absolute values against |Q3|, so it is one-sided in effect.
- *High-frequency step* — FieldTrip FIR high-pass at 225 Hz, two-pass; z-scored with mean/SD; `|z| > 5`. At 500 Hz this passes roughly 170–250 Hz (gain 0.53 at 200 Hz, 0.08 at 150 Hz).
- *Diff step* — z-scored first difference, `|z| > 5`.
- The high-frequency and diff samples are clustered (consecutive samples) and each cluster padded **±100 samples (±200 ms)**. The last cluster of a channel is dropped by a loop bug.

**IEDs** (`findSpikeTimes_reviewed.m`, parameters in `runScripts_review_check.m`). Channel z-scored over the recording (mean/SD). A trough of depth ≥ 4 SD and prominence ≥ 4 SD, width ≤ 50 samples (100 ms, at half prominence), counts when a peak of prominence ≥ 4 SD lies within ±100 ms and trough + peak ≥ `amp_scale`·4 SD. With `amp_scale = 1` that last condition is automatic. **±50 samples (±100 ms)** are marked around each trough. Channels with fewer than **60** detections in the recording get none. The positive-first pass is computed but not used, so only negative-first discharges count. Ripple candidates were also checked against IEDs by eye.

**Paper vs code.** The paper says 12 SD for trough + peak (`amp_scale = 3`), a trough of ≤ 50 ms, a peak *after* the trough, and does not mention the 60-detection gate. The code comment says the IED margin is 1 s; it is ±100 ms. The paper cites Vaz et al. and a "standard approach" (its ref. 97) without full references in the excerpt we have.

## 2. Version 1 — our procedure until 2026-09-23 (`cleanPerisaccadicEEG`, per epoch set)

Run separately on each epoch set: the saccade union set (onset − 1.25 s to offset + 1.25 s, one pass shared by both alignments) and the image set ([−1.25, 2.25] s). Statistics per channel are pooled over the set's complete epochs. Perisaccadic epochs overlap heavily, so a sample is counted once per epoch it falls in.

- *Amplitude* — raw, **not demeaned** values; band `[Q1 − 2.3·IQR, Q3 + 2.3·IQR]`; flagged if any **0.25 s** moving mean inside the epoch leaves it. The saved description said "0.5 s" and "detrended samples"; both are wrong.
- *High-frequency* — 4th-order Butterworth high-pass at `min(250, fs/2 − 1)` = **249 Hz**, `filtfilt` per epoch; flagged if any `|x − mean| > 5 SD` (non-robust). At 500 Hz this passes only 249–250 Hz (gain 0.004 at 248 Hz) and flagged ~0.1 % of channel-epochs: effectively inert.
- *Diff* — raw first difference; flagged if any `|d − mean| > 5 SD` (non-robust).
- A channel-epoch is rejected if any sample of the epoch, padding included, trips any step. No IED step, no bad-channel rule.

**What the audit found** (8 recordings): 85.8 % of channel-epochs usable. The diff step does most of the rejecting (13.3 % of channel-epochs; 34 % in NS174_03). Its 5 SD threshold sits inside the bulk of the per-epoch maxima, not at the tail, and the epochs just above and just below it look the same. The events there are sharp biphasic transients of about 20 ms, not technical artifacts. 1–5 % of *kept* channel-epochs contain an IED; in NS167 and NS140_02, dropping them moves the onset ERP of 34 and 22 channels by more than 2 SEM. On average 15 channels per subject lose more than half their epochs, yet their remaining half is kept. Channel offsets drift between trials (in NS167, 71 % of channels by more than 1 within-trial SD, with steps at block boundaries), and the raw-value amplitude step is exposed to that. No step was locked to saccade onset.


## 3. Version 2 — from 2026-09-24 (`get_artifact_marks` + `cleanPerisaccadicEEG`)

**Structure.** Three steps, run once per subject on the continuous reref trials and cached in `data/artifact_marks/`:
1. `mark_artifacts` marks bad samples **without any behavioural information**.
2. `annotate_marks` is the only step that sees saccade times. It *labels* the marks: IEDs locked to saccades vs other IEDs, and saccadic-spike channels. It also sets the bad-channel verdict. It never changes what step 1 marked.
3. `cleanPerisaccadicEEG(data_perisac, marks)` flags a channel-epoch if **any of its samples, padding included, carries a mark**.

The onset, offset and image sets are all judged against the same cached marks, and bad channels are the same in every set. `get_artifact_marks(subjectId, data, content)` does steps 1–2 and caches the result. Since 2026-09-24 no epochs are stored: `get_events` / `get_epochs` compute the flags per event set on demand from the cached marks (building them on first use, and rebuilding them when the saccade table has changed).

**One rule for the artifact steps (Tukey far-out fence).** For each channel and step, the statistic's maximum is taken in every 1 s segment of the trials. A sample is marked where the statistic exceeds

  `fence = exp(Q3 + 3·IQR)` of the log per-second maxima.

So a mark means "more extreme than this channel's ordinary seconds get". The constant (k = 3, Tukey's "far out") is the same for every step, and the fence adapts to each channel's own level of transients. The signal is demeaned per trial first, which removes drift and block offsets.

| Step | Statistic | Marked | Rule |
|---|---|---|---|
| amplitude | abs(0.25 s moving mean) | the window | fence |
| diff | abs(first difference − median) | ± 0.2 s | fence |
| hf | abs(FieldTrip FIR high-pass 225 Hz, two-pass − median) (Michal's filter) | ± 0.2 s | fence |
| ~~burst~~ | *dropped 2026-09-24*: Hilbert envelope, 55–65 Hz (transient line noise the ±2 Hz notch could not remove); the ±5 Hz notch now removes those bursts | ± 0.2 s | fence |
| IED | `detect_ieds`: Michal's detector, both polarities; score = trough depth + peak height (SD) | ± 0.1 s | FDR ≤ 1 % against a phase-randomised null |
| saccadic spike (channel label) | saccade-locked mean of the 30 Hz high-passed signal, −20..+30 ms | – | significant (Bonferroni, α 0.01) **and** size beyond the fence across channels (Q3 + 3·IQR, linear) |
| line noise (channel label) | 57–63 Hz power left inside the ±5 Hz notch, over 40–45 / 75–80 Hz | – | fence across channels |
| bad channel | fraction of 2.6 s segments holding an excluding mark (artifact or other-IED) | whole channel | > 0.5 |

The ± 0.2 s margin is Michal's.


**IEDs.**
- Candidates must pass Michal's shape criteria: trough and matching peak each ≥ 4 SD prominent, trough ≤ 100 ms wide, peak within ± 100 ms. A candidate inside a 60 Hz burst is a burst, not an IED. Detections within 100 ms are merged into one event.
- The threshold on the score is the lowest at which a phase-randomised copy of the same trials, pooled over channels, gives an estimated FDR ≤ 1 % (5 % until 2026-09-24). It is set per subject.
- Why 1 %, not 5 %: the phase-randomised null keeps each channel's spectrum but flattens its amplitude bursts, so strong oscillation bursts also beat it. At 5 % (subNS167) the kept detections piled up just above the threshold, and on several Vx channels of the LHi shaft their average showed no spike, only a slightly larger oscillation peak; those channels lost 15–40 % of their epochs. At 1 % the threshold rises from 6.75 to 8.25 SD and the Vx loss falls from 9.7 % to 3.7 % of epochs. A null that keeps the bursts, or a morphology criterion, would be the fuller fix.
- The data-driven thresholds agree with Michal's code value (4 SD) far better than with his paper's 12 SD. Above about 12 SD the detector mostly catches 60 Hz bursts, which is why bursts are now their own step.

**Saccade-locked vs other IEDs (event level).**
- For each channel, IED times are cross-correlated with saccade onsets (20 ms bins, ±1 s). Each bin in −0.3..0.5 s gets a Poisson p-value against the channel's baseline (|lag| ≥ 0.5 s), and BH-FDR 5 % is applied over all channels × bins.
- A channel's **locked window** is the run of significant bins around its most significant one. An IED lying in that window after a saccade is `ied_locked`; every other IED is `ied_other`.
- `marks.locking(ch)` holds the window, the counts, and `n_expected`: how many in-window events the baseline predicts, i.e. how many of the "locked" ones are probably ordinary IEDs.
- In NS167, ROI3-ROI4, RHi1-RHi2 and LHi2-LHi3 have 37–44 locked events at +130–150 ms, where chance predicts 3–5. Their waveforms look like the other IEDs, which is why the two groups are kept apart rather than all rejected.

**What noiseInfo holds.**
- Per channel-epoch: `.amp .diff .hf .ied_other .ied_locked .ied`, and two verdicts:
  - `.ok` = complete, no artifact, no **other** IED, channel not bad;
  - `.ok_strict` = `.ok` and no locked IED either.
- Per channel: `.bad_channel .bad_frac .ied_locked_channel .saccadic_spike_channel .saccadic_spike_z .line_noise_channel .ied_rate_per_min .rejected_frac`.
- `params.artifactDetection` holds `version = 2`, the rule, a description of every step, the subject's IED threshold and all parameters. A file without `version` is version 1 and gets cleaned again.


## 4. Cleaning report (per subject, for reporting and QC)

`run_cleaning_report.m` → `cleaning_report(data_name)`. Per subject it saves the numbers (`results/cleaning_report/<data_name>_cleaning.mat`, with a per-channel table, and one row of `results/cleaning_report/cleaning_summary.csv`) and draws two figures in `figs/cleaning_report/`. It reads the cached marks, so once they exist a subject takes about 40 s.

**Which events count.** Saccade events are those every analysis uses (onset and offset share one trial set):
- not the fixation-0 row of a trial;
- not blink-adjacent: `blinkSac` is now set when an invalid eye sample lies within 10 ms of the saccade (`blink_proximity_ms`; the old test, for a blink *inside* the saccade, came from EyeLink-era code and could never fire on this Tobii data — `notes/first_fixation_filter_fix.md`, `standalone/add_blink_proximity.m`);
- complete (the whole epoch inside its recording trial, saccade ≤ 100 ms).

A channel-epoch is kept when `noiseInfo.ok`.

**`_cleaning_overview.png`.**
- Event funnel: saccades → not blink-adjacent → complete (= usable).
- Epochs kept per channel: a histogram over channels of the % of usable epochs each channel keeps; bars stacked by the channel's main cause of exclusion (bad channel, amp, diff, hf, IED other, or nothing excluded).
- Channel × recording-trial map of the % of time marked (amp, diff, hf, IED other).
- IED threshold: detections per channel-minute of the data and of the phase-randomised null against the detector score, with the FDR 1 % threshold.
- Key numbers.
- The line-noise channel label is not shown: it is relative (the channels whose 57–63 Hz power left after the ±5 Hz notch is an outlier among the subject's channels), so it always picks some channels, and it excludes nothing.

**`_cleaning_erps.png`.** Six channels picked from the data: the two most excluded, the highest IED rate, the strongest saccadic spike, two typical ones (or the bad channel, if there is one). For each, saccade-onset-locked:
- top: the ERP (µV, baseline −500 to −250 ms, mean ± SEM) of all epochs, of those excluded for an artifact or IED, of the blink-adjacent ones, and of the kept ones;
- bottom: every epoch, grouped in the same order (coloured bar at the side, black line between groups), colour = amplitude in robust SD units.

**Robust SD.** The spread of a signal measured with medians instead of means, so that a few huge values (artifacts, IEDs) hardly change it:
- MAD, the median absolute deviation: take the median of the values, the distance of every value from it, and the median of those distances — `MAD = median(|x − median(x)|)`.
- For normally distributed data MAD = 0.6745 × SD, so `robust SD = MAD / 0.6745 = 1.4826 × MAD` estimates the SD the data would have without their outliers. 1.4826 is this fixed constant, not a value fitted to our data.
- In the epoch images it is computed per channel, from all of that channel's epochs shown in the panel (kept and excluded together, −500 to 500 ms, after the baseline subtraction). Each channel's image is divided by it, so ±4 on the colour scale means four times that channel's typical spread, whatever its µV level.
- The same estimate appears elsewhere in the cleaning: each phase-randomised surrogate trial is rescaled to its original trial's robust SD before the IED detector runs on it (§ 3).
