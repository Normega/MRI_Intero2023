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

Participants (N = 43 analysable) breathed with a visual pacer whose rate
increased, decreased, or held steady (step vs ramp salience), interleaved
with free breathing, and reported the change direction with confidence and
arousal ratings. Behaviour and physiology are fully analysed (that
foreknowledge is disclosed in the registration): identification is well above
chance, confidence is metacognitively valid, and "No Change" responses
function as "I don't know". The imaging hypotheses are H1–H6 in the prereg;
two are direct replications of Farb, Zuo & Price (2023, eNeuro) — insula/
opercular DEACTIVATION under exogenous pacing (H1: paced < free, direction
committed), and MAIA-moderated sparing of the ACC (H6: sphere at MNI
−4, +12, +36). The functional images have never been touched beyond fMRIPrep;
keep it that way until step 3.

## 2. Inputs

### On this machine (Linux): the imaging

- BIDS dataset + **fMRIPrep 24.1.1 derivatives** (dataset of record; recently
  re-run). Space: `MNI152NLin2009cAsym`, 2 mm, TR = 2 s.
- **Task naming differs from the covariate sheets.** The BIDS/derivatives tree
  uses `task-breath1` / `task-breath2`; the sheets use
  `task-breathing_run-1` / `task-breathing_run-2`. Mapping:
  `breath1 ↔ task-breathing_run-1`, `breath2 ↔ task-breathing_run-2`.
  (Also `loc1/loc2 ↔ task-localizer_run-{1,2}`, `auth1 ↔ task-authenticity`,
  not needed for this analysis.)

### From the repository (transfer or clone from the Windows/Drive copy)

- `FirstLevel/<sub>_<ent>_desc-model_events.tsv` — onset, duration, condition
  (`paced_aware / paced_unaware / paced_wrongdir / paced_nochange / free /
  rating / subrating`), raw confidence and arousal on breathing epochs.
  Onset 0 = first scanner trigger = volume 1.
- `FirstLevel/<sub>_<ent>_desc-model_confounds.tsv` — one row per volume,
  row-aligned with the fMRIPrep confounds: motion 6+6, FD, CSF, WM, cosine*,
  spike_NN (FD > 0.5 mm), and physio (`rvt_rrf_z`, `rvt_z`, `resp_rate_c`,
  `petco2_hrf_lag00/10/20/30`), demeaned, NA→0. Conventions: `FirstLevel/README.md`.
- `FirstLevel/firstlevel_manifest.csv` — **the run list is this file.**
  `motion_gate_fail == FALSE` → 84 runs, 43 participants. Do not re-derive
  inclusion from anything else.
- `Data/maia_scores.csv` — `maia_ieat_total` (0–5 metric) for H6.
- `Data/prereg_sample.csv`, `Data/qc_ledger.csv` — reference only.

## 3. Verification before modelling (do these, report results)

1. Every manifest-retained run has a `desc-preproc_bold` + `desc-brain_mask`
   in MNI space and a confounds TSV whose **row count equals `n_vol`** in the
   manifest. Flag any mismatch; do not silently drop or trim.
   Watch for sub-065 and sub-076 specifically (late additions to the dataset)
   and confirm sub-082 / sub-104 are NOT in the run list (excluded; shifted
   run labels).
2. Events sanity: per run, `events_end_s ≤ scan_dur_s + 10`; all seven
   condition labels present across the set.
3. Design-matrix QC (allowed pre-freeze: this touches no BOLD data): build
   all design matrices; report per-run VIFs among task regressors and the
   correlation of confidence-modulator with the aware-vs-unaware contrast
   vector. High collinearity is a finding to report, not to fix ad hoc.

## 4. The registered model, exactly

**First level, per run** (SPM-style canonical HRF + temporal derivative;
6 mm FWHM smoothing; drift handled by the cosine columns — do not add another
high-pass):

- 7 condition regressors (boxcars over epoch durations) as listed above.
- Parametric modulators on paced epochs: confidence (primary), arousal
  (secondary model), mean-centred **over the trials actually modelled**.
- All confound columns from the sheet. **Never global signal.**
- Sensitivity variants: base / +RVT (`rvt_rrf_z`) / +PetCO2 (4-lag basis,
  all four columns or none) / +both.

**Combining runs:** fixed effects within participant (average of run
contrasts weighted by design efficiency, or SPM-style FFX contrast across
runs). Run as a random factor is registered OUT.

**Contrasts (first level):**
| name | definition |
|---|---|
| `paced_gt_free` | (aware + unaware + nochange)/3 − free  *(H1: prediction is NEGATIVE in insula/operculum, somatomotor)* |
| `confidence_mod` | confidence parametric on paced epochs (H2) |
| `aware_gt_unaware` | aware − unaware (H3) |
| `salience_high_gt_low` | high − low salience change trials (H4; needs the events split by salience — column present in events file) |
| `arousal_mod` | arousal parametric (secondary) |

**Group level:** one-sample t on participant contrast images (N = 43 for
H1/H2; 21 for H3; 41 for H4). H6: add `maia_ieat_total` (demeaned) as a
covariate on `paced_gt_free`; the MAIA effect is the test, one-sided
(higher MAIA → less deactivation, i.e. positive covariate effect).

**ROIs:**
- Harvard–Oxford (50% threshold): L anterior insula, R anterior insula,
  dACC/MCC, posterior insula, SMA, central operculum.
- H2/H3 inference set: L AI, R AI, dACC (Bonferroni ×3).
- H6: 10 mm sphere at MNI (−4, +12, +36) primary; atlas dACC secondary.

**Inference:** whole-brain cluster-forming p < .001, cluster FWE p < .05
(permutation-based cluster inference is acceptable and preferred if using
nilearn); ROI small-volume FWE p < .05 with the stated Bonferroni. H6
one-sided. **No uncorrected result is a finding.**

**H5 (equivalence):** extract per-participant ROI estimates of
`confidence_mod` from base vs +PetCO2 (N = 22) and base vs +RVT (N = 26)
models; TOST on the paired difference, bounds dz = ±0.65 / ±0.59; fixed
interpretation grid in prereg §12.4. Report the corrected–uncorrected
correlation alongside.

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
  build any contrast on "detected vs not" — the yes/no split measures
  response bias here, not sensitivity.
- Do not re-derive run inclusion, motion gating, or subject exclusion —
  the manifest is authoritative.
- Do not put global signal in any model.
- Do not quote uncorrected or lag-/ROI-selected statistics as effect sizes.
