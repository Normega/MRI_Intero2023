# Localizer stage C2 implementation notes: connectivity (written before any connectivity statistic was computed)

Written 2026-11-27, after stage C1 was committed (3ae7ee9) and **before any
localizer connectivity value, gPPI map or T5 coupling value was computed**.
Stage C2 is **exploratory** throughout (plan §9: "connectivity and gPPI", "T5").
It covers plan §6 (network FC, background FC, gPPI) and T5 (§8). T4, decoding,
RSA, gradients and phase-locked responses are part 3 and are not run here.

Plan of record: `Preregistration/intero2023_localizer_prereg.md` v1.0 (d29136b),
§6, §8 (T5), §9, §11; OSF form `Preregistration/osf_localizer_secondary_prereg.md`
(Transformations 13, Volume-level 7). Readings of ambiguous text and every
analysis not listed in §6/§8 are logged in `Preregistration/localizer_deviations.md`
before running. The swc7p H7 script, the breathing-run pipeline and the stage-A/B
localizer scripts are imported, never modified.

**What was run before these notes (no statistic):** design-level checks only.
- gPPI design check on sub-067 (`Analysis/localizer_gppi_design_check.json`).
- One-participant smoke test of the FC and T5 code (sub-067). Only row counts,
  window volume counts and the H7 reproduction difference were looked at.
  The T5 network pairs reproduce the registered H7 per-subject Fisher z for
  sub-067 to 1.2e-12, with identical volume counts.

## Scripts
| Script | Role |
|---|---|
| `Scripts/intero2023_loc_connectivity.py` | Localizer network FC and background FC (`subject`), T5 in the breathing runs (`t5subject`), checks (`verify`), group tests (`group`), figure (`figure`) |
| `Scripts/intero2023_loc_gppi.py` | gPPI design check (`check`), first levels (`subject`), FLAME 1 (`flames`), clusters and within-network summaries (`tables`) |

Per-participant outputs go to `derivatives/intero2023_glm/loc_conn/` (FC, T5),
`derivatives/intero2023_glm/gppi/<sub>/` (gPPI maps) and
`derivatives/intero2023_glm/group/LOC_gPPI_*` (group maps). The repository gets
group tables only: `Analysis/results/EXPLORATORY_LOC_C2_*.csv`.

## 1. Network FC (plan §6, first bullet), as H7
Every step is the registered H7 connectivity model (swc7p §7 final subsection;
`intero2023_h7_connectivity.py`), applied to the localizer runs.
- **Data:** the 102 usable localizer runs; fMRIPrep MNI152NLin2009cAsym res-2
  preprocessed BOLD, **unsmoothed**.
- **Networks:** CAREN `ROIs/net_{SAL,SMN,ECN,DMN,VIS}.nii`. Network signal =
  **median** over voxels at each volume (mean: sensitivity).
- **Denoising per run** (H7 step 2): zero-phase Butterworth band-pass 0.01–0.1 Hz
  (`nilearn.signal.butterworth`, as H7) applied identically to the voxel series,
  every regressor and the fMRIPrep `global_signal`; OLS on the filtered regressors
  plus an intercept. Nuisance = H7's set from the localizer sheet: 6 motion + 6
  derivatives, FD, CSF, WM. No cosines (the band-pass supersedes them) and no
  spike columns (FD > 0.5 mm volumes are dropped instead; p5ztc Volume-level 7).
- **Two task models** (reading logged):
  - `total` = **network FC**: nuisance only. Condition-locked evoked responses stay
    in the signal.
  - `background` = **background FC on task-model residuals**: nuisance plus the
    localizer task regressors (`intero`, `extero`, `match`, `self`, `instruction`,
    `intero_nc`; Rest implicit), SPM canonical + temporal derivative, as H7's
    Model A removal. This is the analysis that is methodologically identical to
    H7, so comparisons with H7 use it.
- **Windows** (H7 step 3): volume times TR (k + 0.5); an epoch contributes the
  volumes in [onset + 4 s, onset + duration + 4 s), FD > 0.5 mm volumes dropped,
  each network's segment demeaned within the epoch, segments concatenated over
  blocks and runs within a condition. Conditions: `intero` (compliant blocks only;
  `intero_nc` blocks are not used), `extero`, `match`, `self`, `rest`. Instruction
  and fixation periods are not used. In sub-067: 57–60 volumes per condition.
- **Coupling** (H7 step 4, `H7.coupling`): Pearson r per pair, condition and
  participant, Fisher z. All 10 pairs.
- **Contrasts and sets.**
  | Contrast | Conditions | Set | N |
  |---|---|---|---|
  | C1 | Intero − Extero | Intero set | 48 |
  | C2 | Intero − Self | Intero set | 48 |
  | C4 (added) | Match − Extero | all | 54 |
  | C5 (added) | Match − Intero | Intero set | 48 |

  A participant enters a contrast if both conditions have ≥ 2 windows and ≥ 20
  volumes (analyst rule; the registered compliance rule already implies it).
  Rest enters the condition-mean matrices only, never a test.
- **Test:** two-sided paired t on the per-participant z difference; dz = mean / SD
  of the difference.
- **Primary pairs** (plan §6): C1 SAL–SMN and SAL–DMN, Holm ×2 within C1; C2
  SAL–DMN (the plan's "DMN–SAL"; H7's pair label is `SAL-DMN`), one test.
- **All 10 pairs:** BH-FDR over the 10 pairs, within each contrast × analysis ×
  measure × specification (as H7 step 8).
- **Checks reported alongside** (H7 step 6):
  - (i) **global signal:** `global_signal` filtered and regressed identically and
    partialled out of each pair correlation within the concatenated windows
    (`z_gs`). Reading rule borrowed from the H7 gate: a pair difference (p < .05)
    is described as network coupling only if it stays p < .05 in the same
    direction after the partial; otherwise "not separable from the global
    signal". Exploratory; it decides wording only.
  - (ii) five-network mean partialled out (`z_five`).
  - (iii) **amplitude:** each network's temporal SD per condition, paired test per
    contrast.
  - (iv) **co-activation:** the network-mean stage-A joint 6 mm base effect of
    the same contrast, one-sample t.
  - (v) **physiology:** `rvt` (+`rvt_rrf_z`) and `co2` (+4-lag PetCO2 basis)
    added to the regression, in participants with the regressor on every used run
    (RVT: 36 overall, 32 Intero; PetCO2: 26 / 24). Each is compared with the base
    model **in the same participants**.
- **Specifications:** summary (median / mean) × shift (4 / 6 s), for both task
  models, as H7a. H7a's FIR task-removal specification is not run: with two
  blocks per condition per run, a per-run FIR basis would nearly saturate each
  block.
- **CSF regressor (added sensitivity).** A check outside this stage found the
  fMRIPrep CSF regressor in the nuisance set to be task-locked: it falls in every
  task block relative to Rest (dz −0.3 to −0.6). Regressing it removes some
  condition-locked variance from every network and from the global signal. The
  primary analysis keeps it, as H7 and the registered nuisance set do. A `nocsf`
  specification (H7 nuisance set without `csf`; median, 4 s) is reported for both
  task models, mainly to qualify the amplitude and global-signal readings.

## 2. Background FC (plan §6, second bullet)
Section 1 with `task = background`. Every test, check and specification above
is computed for both task models and reported side by side.

## 3. gPPI (plan §6, third bullet)
- **Seeds (10).**
  - Gate-primary fROIs (stage B, 10% rule): the individual C1 fROI in AInsL,
    dACC and PIns; the LOSO group fallback in AInsR and the IEAT sphere (read
    from `LOC_L4_gate.csv`, asserted). They exist for the 44 fROI-eligible
    participants only.
  - A-priori masks: CAREN AInsL, AInsR, dACC, PIns and the IEAT ACC sphere (the
    five constraint parcels).
- **Data:** the stage-A first-level data: 6 mm FWHM within the group mask,
  `mean_scaling` (PSC) per run.
- **Seed signal:** the mean over the seed's voxels, per run, adjusted for the run's
  non-task columns (the nuisance sheet without `cosine_*` incl. spikes, the 200 s
  DCT set, the intercept) by OLS.
- **Deconvolution (registered, p5ztc Transformations 13: "gPPI interaction terms use
  the deconvolved seed signal").** Per run:
  - a microtime grid of 16 bins per TR (0.125 s);
  - the neural signal is represented by a DCT basis with n_vol components on that
    grid (SPM-style);
  - each basis column is convolved with the SPM canonical HRF (nilearn `spm_hrf`
    at TR/16) and sampled at TR (k + 0.5);
  - ridge solution, λ chosen by generalised cross-validation over
    10^[−6, 1] × s_max².

  Check (`Analysis/localizer_gppi_design_check.json`): on a simulated neural signal
  the low-passed recovery is r .997 / .98 / .87 at SNR ∞ / 4 / 1. On sub-067's
  seeds the reconvolved fit is r ≈ .69–.80 with λ ≈ 0.4–1.6 × s_max²: GCV treats
  much of the seed signal as noise. None hits the grid edge. The choice of
  deconvolution is registered, so it is not a free choice here.
- **PPI terms:** for every modelled condition c (`intero`, `extero`, `match`,
  `self`, `instruction`, `intero_nc` where present), the deconvolved signal ×
  the 0/1 boxcar of c on the microtime grid, convolved with the same HRF and
  sampled at TR (k + 0.5). Rest is the implicit baseline and has no PPI term (the
  seed regressor carries it), as in the task model. Check: the same convolution
  of each boxcar reproduces nilearn's task regressor at r ≥ .9997.
- **Model:** the stage-A joint two-run base design (task regressors, nuisance,
  spikes, 200 s DCT, intercept; AR(1), `intero2023_firstlevel.fit`) plus the
  seed regressor (run-specific) and the PPI terms (shared across runs, as the
  task regressors). Same nuisance set as stage A, 6 mm.
- **Contrasts:** gPPI-C1 = `ppi_intero` − `ppi_extero`; gPPI-C4 = `ppi_match` −
  `ppi_extero`. Effect, variance and t are saved.
- **Sets:** a-priori seeds: C1 in the Intero set (48), C4 in all 54. fROI seeds:
  both contrasts in the 44 fROI-eligible participants.
- **Group:** FSL FLAME 1 (one-group intercept; effect, variance, dof image of each
  participant's gPPI joint-model dof; group mask), as stage A. `smoothest`, then
  FSL `cluster` at z > 3.09 on z and −z, GRF cluster p × 2 (two-sided), as the
  stage-A descriptive maps (`intero2023_loc_group.clusters`, Harvard–Oxford and
  cerebellar labels at peaks).
- **Multiplicity:** each of the 20 maps (10 seeds × 2 contrasts) is reported at
  its own two-sided cluster FWE. No further correction is applied. The number of
  clusters that would also survive ×20 is stated.
- **Within-network summaries:** per participant, the mean gPPI contrast effect over
  each CAREN network with the seed's own voxels removed; one-sample t and dz;
  BH-FDR over the 50 seed × network tests within each contrast.
- **Physiology rule, fixed now:** gPPI is fitted without physiology (as the
  primary localizer model). If a seed × contrast map has a surviving cluster,
  that seed is refitted with +`rvt_rrf_z` in the participants with RVT on every
  used run, and the cluster's mean effect is compared with the base model in the
  same participants. Otherwise no +RVT gPPI is run.

## 4. T5: fROI-node coupling in the breathing runs (plan §8)
- **Pipeline:** H7 exactly, on the 90 swc7p breathing runs: unsmoothed BOLD,
  band-pass, H7 nuisance set plus the seven Model A task regressors (canonical +
  derivative), the H7 windows (`H7.windows`: paced = all four paced conditions;
  free; aware; unaware), shift 4 s (6 s sensitivity), median (mean sensitivity).
- **Nodes:** the participant's five gate-primary fROIs (as gPPI) and the five
  CAREN networks.
  - fROI–network pairs use the network with the participant's fROI voxels (union
    of the five) removed, so no pair shares voxels.
  - fROI–fROI pairs: 10. The dACC parcel and the IEAT sphere share 209 voxels.
    For the dACC–sphere pair, voxels in both fROIs are removed from both.
  - Network–network pairs (full networks) are recomputed as well. They must equal
    the registered H7 per-subject values (checked in `verify`).
- **Pairs and family:** 25 fROI–network + 10 fROI–fROI = 35; BH-FDR over the 35
  within each contrast × measure × specification.
- **Contrasts and sets:**
  - T5a paced − free: fROI-eligible ∩ the 47 breathing participants = **41**;
  - T5b aware − unaware: fROI-eligible ∩ S7 (the H7b set) = **21**.
- **Test:** two-sided paired t on the z difference, dz. GS partial (`z_gs`) and
  five-network mean partial (`z_five`) as H7; node temporal SD (amplitude).
  +RVT (S2) and +PetCO2 (S1), as in H7, against the base in the same participants.
- **Comparison with H7:** H7's 10 network pairs, primary specification, recomputed
  in the same 41 / 21 from the registered H7 per-subject files
  (`derivatives/intero2023_glm/h7/`), with FDR over 10.

## 5. Physiology and the global signal (summary of the rules above)
- +RVT and +PetCO2 versions of network FC, background FC and T5, each compared
  with the base model in the same participants.
- The H7 global-signal partial for every FC and T5 test.
- `nocsf` sensitivity for localizer FC (CSF regressor task-locked).
- Breathing differs by condition (Intero deeper than Extero, RVT dz .88; Match
  deepest). Amplitude and global-signal differences between conditions are read
  with that in mind.

## 6. Outputs
- `EXPLORATORY_LOC_C2_netfc.csv`: every localizer FC test (analysis =
  total / background / coactivation; measure = z / z_gs / z_five / sd).
- `EXPLORATORY_LOC_C2_netfc_condition_means.csv`: condition means for the matrices.
- `EXPLORATORY_LOC_C2_T5.csv`, `EXPLORATORY_LOC_C2_T5_H7_same_participants.csv`.
- `EXPLORATORY_LOC_C2_gppi_clusters.csv`, `EXPLORATORY_LOC_C2_gppi_networks.csv`,
  `EXPLORATORY_LOC_C2_gppi_firstlevel_qc.csv`.
- `Analysis/localizer_stageC2_verification.json` (H7 reproduction, volume counts).
- `Figures/localizer/loc_C2_network_fc.png`, and a gPPI map figure only if a
  cluster survives.

## 7. Compute
First levels run in the background with `xargs -P 8`, `OMP_NUM_THREADS=2`,
absolute script paths and `~/miniforge3/envs/neuro/bin/python`. Logs go to
`derivatives/intero2023_glm/logs/loc_{fc,t5,gppi}_<sub>.log`.
