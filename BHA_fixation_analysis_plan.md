# BHA × representational distance between fixations: analysis plan

Goal: check whether the distance between neural-network representations of consecutive fixations (DL_fixation_analysis, `pure_sim` / `incremental_sim`) is also reflected in hippocampal BHA, and whether ripples drive the effect.

Legend: ✅ keep as it is · 🔧 keep with a change · ⏸ later

---

## 0. Computing BHA (reference: Leszczyński et al. 2020, Sci. Adv.)

- [ ] Bipolar montage + 50/100/150 Hz notch (check whether the 150 Hz notch is there, since it sits on the edge of the band)
- [ ] Band-pass 70–150 Hz → Hilbert → **power = |analytic|²**
- [ ] Compute on the **continuous signal**, then cut around the fixations
- [ ] Normalize **each channel separately** (z-score over the session, or % change vs ITI)
- [ ] Control: multi-band (10 Hz sub-bands, z-score, average) to account for 1/f
- [ ] Watch the time base: `fixation_beg/end` is at ~500 Hz within the epoch, ripples use a different time base

---

## 1. BHA measure per fixation

| Measure | Rating | Role |
|---|---|---|
| **Before vs after** (end of fixation i−1 vs start of fixation i) | ✅ best | **Primary.** Matches `pure_sim(i)` directly (change i−1 → i) and the RSA result (ripples closer to the previous fixation) |
| Time segments (bins, e.g. 10 ms, locked to fixation onset) | 🔧 | **Secondary.** Instead of arbitrary segments, fit the model in every bin + cluster permutation or FDR → time course of the effect |
| Mean over the whole fixation | 🔧 | **Sanity check.** Fixed window (e.g. 50–250 ms) or the whole fixation + duration as a covariate |

- [ ] Define the windows: "before" = last ~200 ms of fixation i−1 before the saccade, "after" = 50–250 ms of fixation i (skip the first ~50 ms because of the saccadic spike potential)
- [ ] Decide between a difference (after − before) as the outcome and a long format with `window` (before/after) × distance. The long format is more flexible
- [ ] Short fixations: NaN in the missing bins, don't throw out the whole fixation

---

## 2. Channel mean vs per-channel model

| Variant | Rating | Role |
|---|---|---|
| Every channel separately (long format) + `(1 \| subject:channel)` | ✅ | **Primary.** Different numbers of channels per person are fine |
| Mean over hippocampal channels per fixation | ✅ | **Robustness + plots.** If it disagrees with the model, channels are probably heterogeneous |

- [ ] Channel = **random effect nested in subject**, not fixed (channel IDs don't mean the same thing across people)
- [ ] Take the same hippocampal channels as in ripple detection (PR_1–PR_3)

---

## 3. Ripples cut out vs not

🔧 **Keep both**, this is key. In the hippocampus, 70–150 Hz overlaps the ripple band, so without cutting them out BHA ≈ ripple.

- [ ] Variant A: raw BHA
- [ ] Variant B: ripple windows ±50 ms → NaN (per channel, from that channel's detection)
- [ ] Interpretation: effect in A and not in B → the effect is driven by ripples. Effect in both → BHA carries something independent

---

## 4. Ripples as a covariate

🔧 **Keep, but note that it answers a different question than point 3.** Point 3 removes the power of the ripple itself from the signal. Point 4 asks whether *the occurrence of a ripple* around the saccade explains the effect of distance (mediation).

- [ ] Use **ripple presence (0/1) or count in the window of that specific saccade / fixation**, not per trial
- [ ] Compare the `odleglosc` coefficient with and without `ripple` in the model (does it shrink?)
- [ ] Optionally a formal mediation: odleglosc → ripple → BHA
- [ ] Watch out for collinearity: the ripple was detected from the same signal (variant A + covariate is almost circular, so it makes more sense to combine 4 with variant B)

---

## 5. Other regions (ventral visual cortex)

⏸ **Later.** Fusiform / LOC / ITC fits the CNN distance more naturally.

- [ ] Channel selection criterion: BHA increase at image onset vs ITI (independent of the question → no double dipping)
- [ ] Only then add `* rejon` to the model

---

## Model

**Now (hippocampus only). Without `rejon`**, because with one region that term drops out:

```
BHA ~ odleglosc + indx_sac + czas_fiksacji + amplituda_sakady
      + (1 + odleglosc | subject) + (1 | subject:channel) + (1 | image)
```

**+ ripples (point 4):**

```
BHA ~ odleglosc + ripple + indx_sac + czas_fiksacji + amplituda_sakady + (...)
```

**Later (point 5):**

```
BHA ~ odleglosc * rejon + indx_sac + czas_fiksacji + amplituda_sakady
      + (1 + odleglosc | subject) + (1 | subject:channel) + (1 | image)
```

Notes on the model:

- [ ] **Saccade amplitude is required.** It correlates with the distance between crops and drives saccade-locked BHA (spike potential)
- [ ] `indx_sac`: check whether it's linear. If not, use a log or spline (`incremental_sim` has a strong trend with fixation order)
- [ ] If `(1 + odleglosc | subject)` doesn't converge → drop the slope's correlation `(1 + odleglosc || subject)` → a random intercept only
- [ ] Standardize the continuous predictors (easier comparison and convergence)
- [ ] Distance = 1 − cosine similarity

## Which distance metric (so the comparisons don't multiply)

- [ ] Choose a **primary** metric a priori (suggestion: `pure_sim`, because it matches "before vs after" 1:1)
- [ ] The rest (`incremental_sim`, dino vs resnet50, `_cr` layer vs last) as secondary, with FDR correction
- [ ] `_cr` vs last layer = low-level vs semantic dissociation. It's worth showing both in the results

---

## Order of work

1. [ ] BHA on the continuous signal + normalization (point 0)
2. [ ] Join with `fixation_data.csv` (time base!) → long table: fixation × channel
3. [ ] Before/after measure, primary model, variant A (raw)
4. [ ] Variant B (ripples cut out) + model with ripple as a covariate
5. [ ] Channel mean as a check
6. [ ] Time-resolved (bins + cluster permutation)
7. [ ] Mean over the fixation as a sanity check
8. [ ] ⏸ Other regions
