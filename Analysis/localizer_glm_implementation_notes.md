# Localizer GLM implementation notes (written before any localizer statistic was computed)

Written 2026-11-27, after OSF p5ztc was approved and **before any localizer
first-level or group statistic existed**. Stage A only: localizer first levels,
the confirmatory group tests L1, L2 and L3 with the L1 breathing-robustness
co-report, and the registered sensitivity checks for those tests. The plan of
record is `Preregistration/intero2023_localizer_prereg.md` (v1.0, frozen at
d29136b); the OSF form text is `Preregistration/osf_localizer_secondary_prereg.md`.

These are the choices the registration leaves to the analyst. None changes a
registered test. Any later change goes in `Preregistration/localizer_deviations.md`.
The breathing-run pipeline is imported, never modified
(`intero2023_glm_common`, `intero2023_firstlevel.fit`, `intero2023_group`).

## Scripts
| Script | Role |
|---|---|
| `Scripts/intero2023_loc_common.py` | Paths, run list, participant sets, event coding, design builders |
| `Scripts/intero2023_loc_sheets.py` | Nuisance sheets `FirstLevel/localizer/*_desc-model_confounds.tsv` and `localizer_manifest.csv` |
| `Scripts/intero2023_loc_verify.py` | Design-only verification of all 102 runs; one-participant end-to-end check |
| `Scripts/intero2023_loc_firstlevel.py` | First levels, one participant per call |
| `Scripts/intero2023_loc_group.py` | L1–L3, robustness rule, sensitivity analyses, descriptive maps |

Imaging outputs are written to `derivatives/intero2023_glm/firstlevel/<sub>/` and
`derivatives/intero2023_glm/group/LOC_*/`, never into the repository.

## Nuisance sheets (how the breathing sheets were built, and the localizer copy)
The breathing-run sheets were written by three R scripts:
1. `intero2023_covariates_linux.R`: from fMRIPrep's confounds TSV, keep the six
   rigid-body parameters and their `_derivative1`, `framewise_displacement`, `csf`,
   `white_matter`, `cosine*` and `non_steady_state_outlier*`. Non-finite values
   (row 1 of the derivatives and FD) become 0. One `spike_NN` indicator is added
   per volume with FD > 0.5 mm.
2. `intero2023_physio_sheets.R`: `rvt_rrf_z`, `rvt_z` and `resp_rate_c` from
   `Data/regressors/<sub>_<ent>_resp.tsv`. Each is padded or trimmed to n_vol,
   demeaned within run, NA → 0, and omitted if its coverage is below 50%. The
   PMU file is a fallback only when the Biopac `rvt_rrf_z` covers less than 50%.
   The 4-lag PetCO2 basis enters whole or not at all.
3. `intero2023_tissue_confounds.R`: `csf` and `white_matter` re-taken from the
   fMRIPrep derivatives of record.

`intero2023_loc_sheets.py` implements that rule once, in Python.
- **Check:** `validate` rebuilds every committed breathing sheet (103, including
  the gated runs) from the same inputs. All 103 have identical columns and
  identical values (≤ 1.5e-6, the 6-dp rounding)
  (`Analysis/localizer_sheet_validation.csv`).
- `build` writes the 105 localizer sheets and `FirstLevel/localizer/localizer_manifest.csv`.
  The manifest reproduces the registered exclusions: sub-058 run-2 (partial task
  and motion), sub-108 run-1 (motion) and sub-114 run-2 (partial task). That
  leaves 102 runs from 54 participants. The gate and spike counts are asserted
  equal to `Analysis/localizer_qc/motion_gate.csv`.
- **`cosine_*`:** the sheets keep fMRIPrep's `cosine_*` columns, so their
  construction is the breathing construction. The columns are **dropped at model
  time** (`read_confounds`), as registered.
- **Physio:** no PMU localizer files exist, so every localizer RVT is Biopac.
  Usable runs with `rvt_rrf_z`: 75. Usable runs with the PetCO2 basis: 54.

## First level (`Scripts/intero2023_loc_firstlevel.py`)
- **Computation.** `intero2023_firstlevel.fit`: nilearn 0.13.1 `run_glm` (AR(1))
  plus `compute_contrast`.
  - Data: `NiftiMasker` on `ROIs/group_brain_mask.nii`, float32.
  - Scaling: `mean_scaling` per run, giving percent signal change.
  - Each run is masked once per smoothing kernel (6 mm and 4 mm FWHM), and every
    model reuses the same data.
  - Images are written with `g.save_float32`. A variance map with any
    non-positive voxel is an error.
- **Design per run** (`intero2023_loc_common.design`):
  - SPM canonical HRF, no derivative;
  - boxcars over each event's exact duration, clipped to the scan end;
  - frame times TR × (k + 0.5), as the breathing designs and the step-1 QC;
  - nilearn `drift_model="cosine"` at `high_pass = 1/200` Hz, the same DCT set as
    the step-1 efficiency QC;
  - the nuisance sheet without `cosine_*` as `add_regs`;
  - nilearn's intercept.
- **Regressors:** `intero` (compliant blocks, ≥ 4 presses), `intero_nc` (fewer
  than 4 presses, of no interest, present only where such a block exists),
  `extero`, `match`, `self` and `instruction`.
  - `instruction` is every instruction event; its duration already spans the 2 s
    fixation, or only the instruction for v1 Self and Rest.
  - Rest is not modelled: it is the implicit baseline.
  - The coding of `intero_nc` is asserted to agree with the events file's
    `compliant` column.
- **Joint two-run model.** Per-run designs are joined as
  `intero2023_h3_joint_qc.build_joint_design` does, with the localizer task
  columns:
  - task regressors are shared across runs;
  - each run's nuisance, cosine drift and intercept columns are run-specific;
  - each run is PSC-scaled before concatenation, and AR(1) is fitted on the joined
    series, as in H3;
  - residual dof = n − rank(X);
  - for a one-run participant the joint model is that run's model.
- **Per-run models.** The same design fitted to each run alone, saved for the
  split-half analyses of stage B. A per-run contrast is saved only where it is
  estimable: a run whose Intero blocks are all non-compliant has no `intero`
  column, so it has no Intero contrast.
- **Contrasts.**
  - C1 = intero − extero; C2 = intero − self; C4 = match − extero;
    C5 = match − intero; C6 = self.
  - `intero`, `extero` and `match` (each vs Rest, descriptive).
  - For each: effect, variance and t (nilearn's t = effect / √variance).
  - **C3** = voxelwise min(t_C1, t_C2) of the same model. It is saved as a t map
    only, because it has no effect or variance.
- **Fits saved:**

  | Smoothing | Variant | Joint | Per run |
  |---|---|---|---|
  | 6 mm | base (200 s) | yes | yes |
  | 6 mm | `hp128`, `hp320` (DCT cut-off only changes) | yes | no |
  | 6 mm | `rvt` (+`rvt_rrf_z`) | yes, if every used run has it | no |
  | 6 mm | `co2` (+4-lag PetCO2 basis) | yes, if every used run has it | no |
  | 4 mm | base (200 s) | yes | yes (stage B) |

  Names: `<sub>_task-localizer_loc-{joint|run1|run2}_sm-{6|4}_var-<v>_con-<c>_{effect,variance,t}.nii.gz`.
  `<sub>_task-localizer_firstlevel.json` records runs, dofs, columns, clipped
  events and estimable contrasts.
- **Run in background** with `xargs -P 8`, `OMP_NUM_THREADS=2`, absolute script
  paths and `~/miniforge3/envs/neuro/bin/python`. Logs go to
  `derivatives/intero2023_glm/logs/loc_fl_<sub>.log`.

## Participant sets (`participant_sets`)
| Set | Rule | N |
|---|---|---|
| all | ≥ 1 usable run (C4, C6, extero, match) | 54 |
| Intero (L1, L2, C1–C3, C5, intero) | `compliance.csv` `meets_75pct` | 48 |
| L3 | Intero set with `maia_ieat_total` | 46 |
| RVT co-report (L1) and +RVT sensitivity | Intero set, **two usable runs, `rvt_rrf_z` on both** | 29 |
| +PetCO2 sensitivity | Intero set, PetCO2 basis on every used run, as swc7p S1 | 24 (L3 24) |
| No v1 | Intero set minus P51–59 | 40 (L3 subset accordingly) |

- **RVT subset.** The registration says "RVT regressors on both runs (about 29
  of the 44 fROI-eligible)", and it is read literally. Three one-run Intero
  participants have RVT on their only run; they are not in the co-report set.
  - The Intero set's two-run participants are exactly the 44 fROI-eligible, so
    this set is the registered 29.
  - The `rvt` variant is fitted for every participant with RVT on every used
    run (36). The group stage chooses the set.

## Group level (`Scripts/intero2023_loc_group.py`)
- **Engine.** FSL FLAME 1, as H3. Inputs:
  - the joint-model 6 mm effect and variance;
  - a dof image holding each participant's joint-model residual dof;
  - `ROIs/group_brain_mask.nii`.

  Design: a one-group intercept (L1, L2). For L3, an intercept plus demeaned
  `maia_ieat_total`, contrast [0 1].
- **Smoothness.** `smoothest -d <N − p> -r res4d -m <group mask>`. Each mask's
  resels = n_voxels / RESELS.
- **Mask voxel FWE (GRF).** `intero2023_group.grf_p`, the 3-D EC peak p as FSL
  `ptoz -g`, which is one-tailed.
  - **Two-sided** (L1, L3): the mask's peak |z|, with p = 2 × grf_p(|z|).
  - **One-sided** (L2, C5 < 0): the peak of −z, with p = grf_p.
  - Then Bonferroni ×3 for L1 and L2, and ×1 for L3.
  - The voxel count above threshold uses `ptoz` at α / (k × tails).
  - If FLAME fails in a mask (flameo error, or a non-finite z in the mask), the
    permutation result is primary for that mask, as registered.
- **Mask-mean effect.** For each participant, the unweighted mean of the
  joint-model effect over the mask's voxels.
  - L1, L2: group mean and dz (mean / SD).
  - L3: Pearson r and the OLS slope t of the mask means on MAIA.
  - Also, as H3: the REML random-effects meta-analysis of the mask means (with
    the mask-mean variance), and each participant's FLAME weight share
    1 / (varcope_i + σ²_g).
- **Bayes factors** (r = .707).
  - L1, L2: `jzs_bf10` from `intero2023_explore_awareness` on the mask means
    (two-sided). For L2, also the directional BF₋₀ (`bf_plus0` on −means).
  - L3: the JZS regression Bayes factor for one covariate (Liang et al. 2008;
    Rouder & Morey 2012; BayesFactor's `linearReg.R2stat`), from R² of the mask
    means on MAIA, with prior scale .707.
- **Permutation check.** nilearn `non_parametric_inference` within each mask:
  10,000 sign-flips (L3: permutation of MAIA with an intercept),
  `random_state = 20261126`, voxel max-t FWE.
  - Two-sided for L1 and L3.
  - One-sided on negated maps for L2.
  - Same Bonferroni. Peak t from `intero2023_group.ols_t`. Equal-weighted OLS.
- **Breathing-robustness rule (L1).** In the 29, the base and the `rvt` joint C1
  are compared per mask.
  - ratio = group mean of the +RVT mask means / group mean of the base mask means.
  - **breathing-robust** if the sign is the same and ratio ≥ 0.50; otherwise
    **breathing-dependent**.
  - Both FLAME tests in the 29 are reported alongside, descriptively.
  - The rule is applied and reported for each mask, whatever the L1 verdict.
- **Sensitivity analyses** for L1, L2 and L3. Each reruns the same FLAME mask
  test, and reports N, peak z and p ×k:
  - **High-pass:** the `hp128` and `hp320` refits; same participants.
  - **No v1:** P51–59 removed.
  - **Group covariates,** each mean-centred and added as a FLAME column; the test
    stays on the intercept (L1, L2) or on MAIA (L3). Participants without the
    covariate are dropped from that analysis, and N is reported.
    - ΔRVT = the participant mean over runs of `rvt_z_mean` (Intero) − `rvt_z_mean`
      (Extero), from the run-level rows of
      `Analysis/localizer_qc/breathing_by_condition.csv`.
    - Extero entrainment = Extero PLV − null PLV, with null = the mean of the
      participant's Rest, Self and Intero PLV against the virtual pacer. This is
      exactly the step-1 summary, computed per participant.
    - Both use only runs that enter the imaging model (the usable runs), so the
      covariate describes the modelled data.
    - Protocol, effect-coded (v1 +0.5, v2 −0.5) and mean-centred (p5ztc
      Transformations 11; plan §5). It is reported with the drop-v1 analysis.
  - **+PetCO2:** the `co2` refit in the 24 with the basis on every used run, and
    the base model in the same 24.
  - **+RVT:** for L2 and L3, the `rvt` refit and the base model in the same
    two-run RVT set.
- **Whole brain (descriptive and L-test maps).** FLAME 1 on the 6 mm joint
  effects of C1, C2, C5 and intero (N 48), and of C4, C6, extero and match (N 54).
  The L3 MAIA map is included.
  - FSL `cluster` at z > 3.09 on z and on −z.
  - GRF cluster p from DLH and volume, two-sided by Bonferroni over tails
    (p × 2), as H3.
  - **C3 group map:** the conjunction-null minimum of the C1 and C2 FLAME z maps
    (Nichols et al. 2005), in the same 48. Clusters where both z > 3.09 (and
    where both z < −3.09) are tabulated with size and peak. A GRF cluster p
    does not apply to a minimum statistic, so none is given. This map is
    descriptive. The participant-level C3 (min-t) exists for the stage-B fROIs.
  - **Labels:** FSL `atlasquery` at each cluster peak, for Harvard–Oxford
    cortical and subcortical and the cerebellar atlas, as in the breathing
    exploratory tables.
- **Outputs.**
  - `Analysis/results/LOC_<label>_{roi,wholebrain,subjects}.csv`;
  - `Analysis/results/LOC_L1_breathing_robustness.csv`;
  - `Analysis/LOCALIZER_RESULTS.md`.
- **Permutation p-values** have a floor of 1e-4 (10,000 permutations).
