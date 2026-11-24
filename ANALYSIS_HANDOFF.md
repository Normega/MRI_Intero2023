# Analysis handoff — first-level fMRI (for the Claude session on the Linux box)

You are picking up the MRI_Intero2023 imaging analysis at the point where
everything upstream of the GLM is finished, frozen, and preregistered. Your
job is to implement the registered first-level and group analysis — nothing
else. This document is self-contained, but the two authoritative references
travel with it in this repository:

- `Preregistration/intero2023_imaging_prereg.md` — the full plan. **Binding.**
- `Preregistration/osf_template_answers.md` — the OSF form text.

## 0. Gate: do not fit anything yet

The registration must be **frozen on OSF before the first GLM is fitted**.
Ask the user to confirm the OSF registration is frozen before you run any
model. Design-matrix construction and input verification (steps 1–2) are fine
before the freeze; parameter estimation is not.

After the freeze, any departure from the plan goes in
`Preregistration/deviations.md`: date, what changed, why, and whether the
affected result had been seen. Post-hoc analyses are labelled exploratory,
no exceptions.

## 1. Context in five sentences

Participants (N = 47 analysable, 90 breathing runs; prereg v0.3) breathed
with a visual pacer whose rate increased, decreased, or held steady (step vs
ramp salience), interleaved with free breathing, and reported the change
direction with confidence and arousal ratings. Behaviour and physiology are
fully analysed (that foreknowledge is disclosed in the registration):
identification is well above chance, confidence is metacognitively valid, and
"No Change" responses function as "I don't know". The imaging hypotheses are
H1–H6 in the prereg; H1 is a direct replication of Farb, Zuo & Price (2023,
eNeuro) — insula/opercular DEACTIVATION under exogenous pacing (paced < free,
direction committed) — and H6 is a **two-sided conceptual extension** of that
paper's MAIA–ACC result (sphere at MNI −4, +12, +36). The functional images
have never been touched beyond fMRIPrep; keep it that way until step 3.

## 2. Inputs

### On this machine (Linux): the imaging

- BIDS dataset + **fMRIPrep 24.1.1 derivatives** (dataset of record; recently
  re-run). Space: `MNI152NLin2009cAsym`, 2 mm, TR = 2 s.
- **Task naming matches the covariate sheets.** The BIDS/derivatives tree of
  record uses `task-breathing_run-1` / `task-breathing_run-2`, identical to the
  sheets; the mapping is the identity. (An earlier version of this file said
  `task-breath1/breath2` — that was the obsolete Windows `Data/Ready` copy.)

### From the repository (transfer or clone from the Windows/Drive copy)

- `FirstLevel/<sub>_<ent>_desc-model_events.tsv` — onset, duration, condition
  (`paced_aware / paced_unaware / paced_wrongdir / paced_nochange / free /
  rating / subrating`), raw confidence and arousal on breathing epochs.
  Onset 0 = first scanner trigger = volume 1.
- `FirstLevel/<sub>_<ent>_desc-model_confounds.tsv` — one row per volume,
  row-aligned with the fMRIPrep confounds: motion 6+6, FD, CSF, WM, cosine*,
  spike_NN (FD > 0.5 mm), and physio (`rvt_rrf_z`, `rvt_z`, `resp_rate_c`,
  `petco2_hrf_lag00/10/20/30`), demeaned, NA→0. Conventions: `FirstLevel/README.md`.
- **The run list is `FirstLevel/firstlevel_manifest.csv` +
  `FirstLevel/firstlevel_manifest_candidates.csv`**, `motion_gate_fail ==
  FALSE`, excluding sub-056 → **90 runs, 47 participants** (prereg §4). The
  candidates file carries the four reinstated participants (sub-065, -076,
  -082, -104). Do not re-derive inclusion from anything else.
- Per-hypothesis participant sets: `analysis/design_qc_twomodel.csv`
  (H3's 23 are also listed in prereg §3).
- `Data/maia_scores.csv` — `maia_ieat_total` (0–5 metric) for H6.
- `Data/prereg_sample.csv`, `Data/qc_ledger.csv` — reference only.

## 3. Verification before modelling (do these, report results)

1. Every manifest-retained run has a `desc-preproc_bold` + `desc-brain_mask`
   in MNI space and a confounds TSV whose **row count equals `n_vol`** in the
   manifest. Flag any mismatch; do not silently drop or trim.
   Watch for sub-065 and sub-076 specifically (late additions to the dataset)
   and confirm sub-082 / sub-104 ARE in the run list with one breathing run
   each (reinstated in v0.2; the "shifted run labels" exclusion applied only
   to the obsolete Windows copy).
2. Events sanity: per run, `events_end_s ≤ scan_dur_s + 10`, **except** the
   three registered clipped runs (prereg §5: sub-052 r1, sub-053 r1, sub-053
   r2) — clip trailing events to `n_vol × TR`, drop any with duration ≤ 0,
   and log the count. All seven condition labels present across the set.
3. Design-matrix QC (allowed pre-freeze: this touches no BOLD data): build
   all design matrices for Models A/B/C; confirm 0/90 rank-deficient in every
   model (`Scripts/intero2023_twomodel_qc.py` reproduces this). High
   collinearity is a finding to report, not to fix ad hoc.

## 4. The registered model, exactly

**First level, per run** (SPM-style canonical HRF + temporal derivative;
6 mm FWHM smoothing; drift handled by the cosine columns — do not add another
high-pass):

Three models per run (prereg §7), identical nuisance set, boxcars over
epoch durations:

- **Model A (awareness → H1, H3, H6):** `paced_aware`, `paced_unaware`,
  `paced_wrongdir`, `paced_nochange`, `free`, `rating`, `subrating`. No
  parametric modulator.
- **Model B (confidence → H2 only):** all paced epochs (incl. No Change) as
  one `trial` regressor + mean-centred confidence modulator; `free`,
  `rating`, `subrating` as in A.
- **Model C (salience → H4):** paced change epochs split `high` / `low`
  salience; `paced_nochange`, `free`, `rating`, `subrating` as in A.
- Arousal is NOT a registered test. If run, it is a separate model of Model
  B's form, labelled exploratory — never added to A or C.
- All confound columns from the sheet. **Never global signal.**
- Sensitivity variants: base / +RVT (`rvt_rrf_z`) / +PetCO2 (4-lag basis,
  all four columns or none) / +both.

**No implicit baseline.** Events tile ~94% of every run. Interpret only
between-condition contrasts; never report or plot a condition-vs-rest
contrast or a single-condition beta map (prereg §5).

**Combining runs:** fixed effects within participant (average of run
contrasts weighted by design efficiency, or SPM-style FFX contrast across
runs). Run as a random factor is registered OUT. **For `confidence_mod`,
skip runs where confidence is constant** (modulator all zeros): sub-106
run-1 and sub-112 run-2 contribute their other run only; sub-066, sub-070,
sub-114 are out of H2 entirely.

**Contrasts (first level):**
| name | model | definition |
|---|---|---|
| `paced_gt_free` (C1) | A | unweighted mean of the paced conditions **present in the run** among aware / unaware / nochange, minus free *(H1: prediction NEGATIVE in insula/operculum, somatomotor)* |
| `aware_gt_unaware` (C2) | A | aware − unaware (H3) |
| `confidence_mod` (C3) | B | confidence parametric on `trial` (H2) |
| `salience_high_gt_low` (C5) | C | high − low salience change trials (H4) |

**Group level:** one-sample t on participant contrast images. Ns: H1 47,
H2 44, H3 23, H4 43 (≥ 3 trials per cell pooled across runs; sets in
`analysis/design_qc_twomodel.csv`). H6 (N = 46; sub-065 has no MAIA): add
`maia_ieat_total` (demeaned) as a covariate on `paced_gt_free`; the MAIA
effect is the test, **two-sided** (conceptual extension; prereg §3).

**ROIs:**
- Harvard–Oxford (50% threshold): L anterior insula, R anterior insula,
  dACC/MCC, posterior insula, SMA, central operculum.
- H2/H3/H4 inference set: L AI, R AI, dACC (Bonferroni ×3), one-sided.
- H6: 10 mm sphere at MNI (−4, +12, +36) primary; atlas dACC secondary;
  two-sided.

**Inference:** whole-brain cluster-forming p < .001, cluster FWE p < .05,
two-sided (permutation-based cluster inference is acceptable and preferred
if using nilearn); ROI small-volume FWE p < .05 with the stated Bonferroni.
**No uncorrected result is a finding.**

**H5 (equivalence):** extract per-participant ROI estimates of
`confidence_mod` from base vs +PetCO2 (N = 22) and base vs +RVT (N = 26)
models; TOST on the paired difference, bounds dz = ±0.65 / ±0.60; fixed
interpretation grid in prereg §12.4. Report the corrected–uncorrected
correlation alongside. The same comparison for `aware_gt_unaware` (12 / 14
participants) is **descriptive only** — report estimates side by side, draw
no survival/abolition conclusion.

**Registered diagnostic:** run-1 vs run-2 spatial correlation of
`paced_gt_free` and `confidence_mod` contrast maps within participants +
cross-participant ICC of ROI means. Descriptive only. H3 excluded from it.

## 5. Implementation notes (choices left to you, within the spec)

- **nilearn** is the pragmatic choice here (`FirstLevelModel` with
  `hrf_model='spm + derivative'`, `drift_model=None`, smoothing_fwhm=6,
  the sheet columns passed as the design's confounds;
  `non_parametric_inference` for group cluster FWE). SPM12 on MATLAB is
  equally valid if present — the IEAT used SPM12 — but do not mix engines
  across hypotheses.
- Events: build per-condition onsets from the `_desc-model_events.tsv`
  directly; `salience` and `change_direction` columns are present for the H4
  split. Modulator columns are RAW — centre them at design time.
- Anatomical masks: nilearn's Harvard–Oxford fetcher
  (`cort-maxprob-thr50-2mm`) matches the registered 50% threshold.
- Keep every fitted image under a derivatives-style tree, e.g.
  `derivatives/firstlevel/<sub>/`, with contrast maps named by the table
  above; write a `group/` tree per hypothesis. Log software versions.

## 6. What to send back

1. Verification report (step 3) before any fitting.
2. After the freeze: per-hypothesis results tables (peak MNI, cluster size,
   corrected p), the H5 TOST numbers with the interpretation-grid cell,
   the reliability diagnostic, and the design-QC summary.
3. Any deviation, however small, as a `deviations.md` entry rather than a
   silent adjustment.

## 7. Do-not list (each of these was earned the hard way)

- Do not use the Birn RRF on PetCO2 (wrong sign for CO2) and do not
  normalise any kernel by signed area (inverts biphasic kernels).
- Do not fit or select a PetCO2 transit lag; the 4-lag basis IS the model.
- Do not treat `paced_nochange` responses as correct rejections, and do not
  build any detection-stage contrast ("said changed" vs "said no change"
  across stimulus types, i.e. hits vs false alarms) — that split measures
  response bias here, not sensitivity. H3's aware − unaware is different: it
  is restricted to true-change trials and its aware side requires the
  correct direction.
- Do not re-derive run inclusion, motion gating, or subject exclusion —
  the manifest is authoritative.
- Do not put global signal in any model.
- Do not put the confidence modulator in the same design as the awareness
  regressors (9 designs go rank-deficient).
- Do not report any condition-vs-rest contrast (there is no rest baseline).
- Do not quote uncorrected or lag-/ROI-selected statistics as effect sizes.

---

## Addendum (2026-08-21): sub-082 / sub-104 — you are right, and here is the history

> **Status (2026-11-24): resolved.** The verification this addendum asks for was
> run on the Linux side on 2026-08-20 (`VERIFICATION_REPORT.md`): slot mapping
> correct for both, one gate-passing breathing run each. sub-082 and sub-104 were
> reinstated in prereg v0.2 (`firstlevel_manifest_candidates.csv`), and the Ns
> are now those in §4 above (H1 47, not the 45 projected below). Kept for the
> record only.

You asked why these two are excluded when their imaging looks complete. The
exclusion was made against the **stale Drive-side copy** of the derivatives,
whose task labels for these two participants are shifted. Your copy is the
dataset of record and was re-preprocessed after the labels were fixed — so on
your side they SHOULD look complete and correctly named.

**What actually happened in those sessions:** both acquired only THREE BOLD
runs (localizer-1, ONE breathing run, authenticity; their PMU folders have
exactly 3 recordings, and the task logs contain exactly one completed
breathing run each, fully evented: 7 paced + 7 free epochs). Order-based
naming on the stale copy then slid the labels: the file named `task-loc2`
holds 242 volumes (breathing length) and the file named `task-breath1` holds
172 volumes (authenticity length).

**Verification required before re-inclusion** (report the table back):
for each func run of sub-082 and sub-104 on your side, read the BIDS sidecar
`ProtocolName` / `SeriesDescription` (Siemens `bto_ep2d_bold_p2s2_2mm_runN_*`)
and the volume count. Expected if your labels are correct:

| sidecar protocol | expected task label | expected volumes |
|---|---|---|
| ...run1... | loc1 / localizer | 196 |
| ...run3... | breath1 / breathing run 1 | 242 |
| ...run5... | auth1 / authenticity | ~172 |

Any file whose protocol slot and task label disagree, or whose volume count
does not match its slot, blocks re-inclusion and should be reported instead.

**If verification passes:** commit each participant's breathing-run confounds
TSV to `FirstLevel/incoming/` and push. The Windows-side session will rebuild
their covariate sheets, add them to the manifest, and amend the registered Ns
(the registration is not yet frozen, so this is a pre-freeze correction, not
a deviation).

**Already pre-verified on the Drive side** (from the mislabelled-but-identical
confounds, FD column only): events complete for both; motion gate PASSES
(mean FD 0.242 and 0.154, spike fractions 0.004 and 0.000); MAIA complete for
both. If your verification passes, C1/C3/H6 move from N = 43 to **N = 45**,
each contributing one breathing run.

**While you are there:** spot-check two or three already-retained runs by
comparing your (re-run) confounds row counts against the sheets' — fMRIPrep
motion estimates should be effectively identical across re-runs, but if you
see material differences, say so and ship replacements rather than adapting
silently.
