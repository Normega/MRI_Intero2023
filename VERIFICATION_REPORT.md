# Verification report — ANALYSIS_HANDOFF.md step 3

**Date:** 2026-08-20 · **Machine:** Linux box, `/media/normega/FB_Cluster_Data/intero2024`
**Status: NO MODEL HAS BEEN FITTED. No BOLD image has been read.** Everything below
comes from headers, confound tables, event files and design matrices, which
§0 of the handoff permits before the OSF freeze.

**Software:** fMRIPrep 24.1.1 (`derivatives/fmriprep/dataset_description.json`),
R 4.6.1 (readr 2.2.0, dplyr 1.2.1), Python 3.11 / nilearn 0.13.1, numpy 2.4.6,
pandas 3.0.3, nibabel 5.4.2.

**Headline:** the inputs verify, but three things must be settled before the
freeze — a stale subject inventory that costs 4 analysable participants, a
contrast definition that is undefined in 18 of 91 runs, and a confidence
modulator that is perfectly aliased with the aware/unaware split in 9 runs and
strongly correlated with it in many more. Corrected per-hypothesis Ns are in §5.

---

## 1. Naming, and a script that will not run here

`ANALYSIS_HANDOFF.md` §2 states the BIDS/derivatives tree uses
`task-breath1`/`task-breath2`. **It does not.** Both the BIDS tree and the
fMRIPrep derivatives on this machine use `task-breathing_run-1` /
`task-breathing_run-2` — identical to the covariate sheets. The mapping is the
identity; no translation layer is needed.

The `task-breath1` naming belongs to the Windows `Data/Ready` snapshot that
`Scripts/intero2023_covariates.R` reads. Consequence: that script's
`_task-([a-z0-9]+)_` regex captures `breathing`, `TASK_MAP` returns `NA`, and
the script produces **zero rows** against this tree. `Scripts/intero2023_covariates_linux.R`
is the port; `diff` the two files for the complete change set.

## 2. Step 3.1 — inputs for every retained run

All **84** runs with `motion_gate_fail == FALSE` (43 participants):

| check | result |
|---|---|
| `space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz` present | 84 / 84 |
| `space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz` present | 84 / 84 |
| fMRIPrep `desc-confounds_timeseries.tsv` present | 84 / 84 |
| model events + confounds sheets present | 84 / 84 |
| manifest `n_vol` == NIfTI volumes == fMRIPrep confound rows == sheet rows | 84 / 84, **zero mismatches** |
| TR | 2 s throughout |
| MAIA `maia_ieat_total` present | 43 / 43 |

Nothing was dropped or trimmed.

### 2a. The port reproduces the registered sheets — with one exception

`intero2023_covariates_linux.R` was re-run over the 48 registered subjects and
compared against the committed files:

- `firstlevel_manifest.csv` — **byte-identical**. Every count, `fd_mean`,
  `fd_frac_gt_thresh`, `n_spikes`, `motion_gate_fail`, `has_rvt`, `has_petco2`
  reproduces exactly. The gate, the 84/11 split and N = 43 are all confirmed.
- all 95 `_desc-model_events.tsv` — **byte-identical**.
- all 95 `_desc-model_confounds.tsv` — **differ, in exactly two columns.**

Motion (6+6), FD, cosines, spikes, non-steady-state columns and every physio
column reproduce bit-exactly. `csf` and `white_matter` do not:

| column | median r (committed vs regenerated) | min r | median mean-ratio |
|---|---|---|---|
| `csf` | 0.985 | 0.887 | 0.993 |
| `white_matter` | 0.996 | 0.936 | 0.987 |

These are the mean signals inside the fMRIPrep tissue masks, which depend on
the anatomical segmentation. BOLD-side outputs being bit-identical while
anatomy-derived signals differ means **the committed sheets' tissue regressors
came from a different anatomical pass than the fMRIPrep derivatives on this
box**, which the handoff calls the dataset of record. The Windows snapshot they
were built from cannot be inspected from here, so the cause is not established.

*(Checked and cleared: all 54 BIDS `_T1w.nii.gz` differ from their non-defaced
sources, i.e. defacing is intact everywhere. That is not the explanation.)*

**Recommendation:** regenerate all 95 confound sheets from the derivatives of
record. It is a two-column change to nuisance regressors, it leaves the
manifest and therefore the registered sample bit-identical, and it is free to
do before the freeze. Regenerated files are staged at
`analysis/verify_regen/` and have **not** been written into `FirstLevel/`.

**Adopted 2026-11-25 (prereg v0.6, `methods_decisions.md` decision 15).** `Scripts/intero2023_tissue_confounds.R` rewrote only `csf` and `white_matter` in the 95 sheets. The rewritten values equal the staged files, and no other column changed. Per-run agreement is in `Analysis/tissue_confounds_regen.csv`.

## 3. Step 3.2 — events sanity

**Seven condition labels are present across the set** — `paced_aware`,
`paced_unaware`, `paced_wrongdir`, `paced_nochange`, `free`, `rating`,
`subrating`. The eighth label documented in `FirstLevel/README.md`,
`paced_noresponse`, occurs in **no** run at all, so the handoff's 7-regressor
list is correct as written and the README overstates the set.

**`events_end_s ≤ scan_dur_s + 10` fails for 2 of the 84 retained runs:**

| run | events end | scan | over |
|---|---|---|---|
| sub-053 task-breathing_run-1 | 525.6 s | 500.0 s | **+25.6 s** |
| sub-052 task-breathing_run-1 | 499.0 s | 488.0 s | **+11.0 s** |

In both the over-running events are trailing `free` / `rating` epochs — the
scanner stopped mid-rating. Two further runs over-run by less than the 10 s
tolerance (sub-053 run-2, +5.6 s). Nothing was dropped: `intero2023_designqc.py`
clips such events to the scan end and reports the count, and a rule for this
needs to be stated in the prereg rather than left to the modelling library,
which would silently discard them.

## 4. Step 3.3 — design-matrix QC

91 designs built (84 registered + 7 gate-passing candidates, §6), exactly as
registered: SPM HRF + temporal derivative, `drift_model=None`, sheet columns as
confounds, confidence PM mean-centred over the trials actually modelled.
Median 43 regressors, median residual df 199 (min 144). Full output:
`Analysis/design_qc.csv`.

### 4a. VIF among task regressors

The intercept is excluded from the reference set deliberately: this design
tiles the run (`paced + free + rating ≈ constant`), so a VIF taken against the
constant diverges for every run by construction and diagnoses nothing.

| regressor | runs | median VIF | max |
|---|---|---|---|
| `free` | 91 | 0.72 | 0.85 |
| `rating` | 91 | 0.65 | 0.70 |
| `paced_aware` | 91 | 1.00 | ∞ |
| `paced_nochange` | 90 | 1.12 | ∞ |
| `paced_unaware` | 74 | 1.11 | ∞ |
| `paced_wrongdir` | 51 | 1.19 | ∞ |
| `conf_pm` | 84 | 1.99 | **∞ (12 runs)** |

Typical collinearity is unremarkable. The infinities are the finding.

### 4b. The confidence modulator is aliased with the aware/unaware split

This is the specific check the handoff asked for, and it does not come back clean.

- **r(confidence PM, aware-minus-unaware direction): median |r| = 0.437**,
  90th pct 0.767, max 0.942; **|r| > 0.5 in 26 of 70 runs.**
- In **9 of 91 runs the aliasing is perfect**: the design is rank-deficient, and
  removing `conf_pm` restores full rank in all 9. The confidence PM is a
  perfect linear combination of the condition regressors in those runs.

The cause is confidence having almost no within-run variance. sub-077
run-2 is typical: confidence is 5 on all five aware trials and 4 on the single
unaware and single nochange trial, so the centred modulator *is* the
aware-vs-rest contrast up to scale.

| distinct confidence values | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| runs (of 91) | 7 | 27 | 37 | 14 | 4 | 2 |
| subjects (of 48, pooled over runs) | 3 | 9 | 16 | 14 | 4 | 2 |

**34 of 91 runs, and 12 of 48 participants, use two or fewer distinct
confidence values.** Three participants (sub-066, sub-070, sub-114) use exactly
one across every retained run, so H2 has no modulator to estimate for them.

Reported, not fixed, per the handoff. But it bears on the registration
directly: the prereg makes H2 (parametric confidence) **primary** and H3
(categorical aware/unaware) secondary on the grounds that the design is short.
In a substantial minority of runs those two are not separable, and the
registration should say so before it is frozen.

### 4c. The registered `paced_gt_free` is undefined in 18 of 91 runs

`paced_gt_free = (aware + unaware + nochange)/3 − free` is a **fixed** 3-way
average, so it cannot be formed when one of the three is absent:

- 17 runs have **no** `paced_unaware` trial (the participant never answered
  "No Change" on a change trial)
- 1 run has no `paced_nochange` trial
- 7 further runs have all three but are rank-deficient via `conf_pm` (§4b)

**25 of 91 runs total cannot yield the contrast as literally written.**
Averaging instead over the paced conditions actually present in the run —
`paced_gt_free_adaptive` in the QC output — recovers all but the 9
rank-deficient ones (82/91). SE multipliers are otherwise near-identical
(median 0.0830 vs 0.0840), so this costs nothing in precision. *(Corrected 2026-11-25: the QC script scaled contrasts by the column SDs instead of dividing by them; the corrected medians are 0.247 for both forms, so the conclusion stands. `Preregistration/methods_decisions.md` decision 18.)*

## 5. Corrected per-hypothesis Ns

Group N counts a participant when the contrast is estimable in ≥ 1 of their
gate-passing runs (and, for H6, when MAIA is present).

| | subjects | H1 as written | H1 adaptive | H2 | H3 | H6 |
|---|---|---|---|---|---|---|
| **prereg §4 claim** | 43 | 43 | — | 43 | 21 | 43 |
| **registered set, actual** | 43 | **39** | 42 | **38** | 21 | 42 |
| **registered + 4 candidates** | 47 | 43 | **46** | **42** | **23** | **45** |

The registered Ns of 43 for H1 and H2 are **not achievable** as the model is
written — the shortfall is the contrast-definition and aliasing problems above,
not motion or missing data. H3's 21 is confirmed exactly.

Two independent corrections happen to cancel: adopting the adaptive paced
contrast and admitting the four candidate participants (§6) restores H1 to 46
and H2 to 42, above the registered figure. **They should still be reported as
two separate amendments**, because either one alone changes the numbers.

S1/S2 (the H5 TOST bounds, N = 22 and 26) are **unaffected**: none of the
candidate subjects has breathing-run RVT or PetCO2 (§6), so no physio-bearing
run enters or leaves.

## 6. Candidate participants — sheets generated

`Scripts/intero2023_covariates_linux.R candidates` wrote 8 runs / 5 subjects
into `FirstLevel/` plus **`firstlevel_manifest_candidates.csv`**. The
registered `firstlevel_manifest.csv` is untouched — the handoff's "do not
re-derive run inclusion" is respected, and admitting these runs remains an
explicit prereg amendment rather than a silent absorption.

| subject | runs | gate | events fit | RVT/CO2 | MAIA | verdict |
|---|---|---|---|---|---|---|
| sub-065 | 2 | both pass | both fit | no | **no** | **admit** — H1/H2/H3/H4, not H6 |
| sub-076 | 2 | both pass | both fit | no | yes | **admit** — all hypotheses |
| sub-082 | 1 | passes | fits | no | yes | **admit** — H1/H2/H6 |
| sub-104 | 1 | passes | fits | no | yes | **admit** — H1/H2/H6 |
| sub-056 | 2 | run-2 **fails** (fd 0.471, 26.7 % spikes) | run-1 **over-runs by 28 s** | no | yes | **exclude**, now on stated grounds |

Why they were missing: the repo's imaging inventory predates the 2026-08-06
re-conversion of this dataset. `Data/qc_ledger.csv` still records sub-065
`task-localizer_run-1` as 167 volumes; it is 196 here, and that single row is
the **only** volume-count disagreement between the ledger and this box across
222 checked runs. sub-076 has no imaging rows at all. Both now have complete
5-run sessions (breathing 242/242, localizer 196/196, authenticity 172).

sub-082 and sub-104 are a separate case — see the §4 correction draft.

## 7. What must be decided before the freeze

1. **Admit sub-065, sub-076, sub-082, sub-104?** (§6, and the §4 correction draft.)
2. **Adopt the adaptive `paced_gt_free`?** Without it H1 runs at N = 39–43 rather
   than the registered 43–46, and 18 runs are silently unusable. (§4c)
3. **State H2's exposure to the confidence/awareness aliasing**, and decide what
   happens to the 9 perfectly-aliased runs and the 3 participants with no
   confidence variance. (§4b)
4. **State the trailing-event rule** for the 2–3 runs that over-run the scan. (§3)
5. **Regenerate the 95 confound sheets** so `csf` / `white_matter` come from the
   derivatives of record. (§2a)

Items 2–5 are corrections to the plan's internal consistency and are free to
make now. Item 1 is a sample change and is the one that must not wait.
