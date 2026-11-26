# GLM implementation notes (written before any group result was computed)

Written 2026-11-26, after the swc7p update was approved and while the first
levels were running, **before any group-level statistic existed**. These are
choices the registration leaves to the analyst (handoff §5). None changes a
registered test. Any later change to them goes in `Preregistration/deviations.md`.

## First level (`Scripts/intero2023_firstlevel.py`)
- **Computation.** nilearn 0.13.1 `run_glm` / `compute_contrast`, done once per
  run on the same smoothed data for every model and variant. On sub-082 run-1
  this reproduces `FirstLevelModel` exactly: identical effect and variance,
  max |difference| 0. Masker: `ROIs/group_brain_mask.nii`, 6 mm FWHM, float32
  data. Then per-voxel percent signal change (`signal_scaling=0`) and AR(1).
- **Designs.** Built by the same functions as the pre-registration design QC
  (`intero2023_twomodel_qc.model_a/b/c`). Frame times are TR × (k + 0.5), as in
  the QC.
- **Combining runs (C1, C3, C5, and the per-run C2 check).** Voxelwise
  inverse-variance (precision-weighted) fixed effects. The handoff allows this
  ("SPM-style FFX").
- **C3.** Only runs whose paced confidence varies enter, so sub-106 run-1 and
  sub-112 run-2 are skipped.
- **H3 joint model.** `build_joint_design()` on the per-run Model A designs.
  Each run is PSC-scaled before concatenation, and AR(1) is fitted on the
  joined series. The residual dof is n − rank(X). For a one-run participant the
  joint model is identical to that run's Model A (checked: difference 0).
- **Physio variants** (+rvt, +co2, +both) are fitted on every run that has the
  regressor. The group stage chooses participants per the registered sets:
  - S1: PetCO2 on every run (25).
  - S2: Biopac RVT on every run (30).
  - S2b: Biopac or PMU RVT on every run (41).
  In the H3 joint model the physio columns are nuisance, so they are
  run-specific.

## Group level: H1, H2, H4, H6 (`Scripts/intero2023_group.py`)
- **Engine.** nilearn `non_parametric_inference`: sign-flip permutation (for H6,
  permutation of the MAIA covariate with an intercept), **10,000
  permutations**, `random_state = 20261126`.
- **ROI tests.** Voxel-level max-t FWE within each registered mask, one-sided in
  the predicted direction, and two-sided for H6. Bonferroni ×3 is applied as
  p_FWE × 3 over the hypothesis's three masks. H1's predicted direction is
  negative, so its maps are negated for the one-sided test. The opposite tail
  is tested the same way, as the registered failed-replication check.
- **Whole brain.**
  - The cluster-forming threshold is **p < .001 one-tailed in each tail**: t at
    p = .001 for the group df, which nilearn gives as `threshold=0.002` with
    `two_sided_test=True`. This is SPM's convention (the IEAT used SPM12) and
    matches H3's z > 3.09.
  - Cluster-size FWE p < .05 comes from the maximum cluster over both tails,
    so it is two-sided.
- **Reported alongside, never as findings:** peak t and its MNI coordinate in
  each mask, the number of voxels surviving, and the a-priori-mask mean
  contrast with its dz (descriptive).

## Group level: H3 (FSL FLAME 1)
- `flameo --runmode=flame1` on the C2 joint-model effect, variance and
  per-participant dof images, one-group design, `ROIs/group_brain_mask.nii`.
- **Smoothness.** `smoothest -d <dof> -r res4d -m group mask`, giving DLH and
  the resel size. Each mask's resel count is n_voxels / resel size.
- **ROI voxel FWE.** The GRF peak p for the maximum z in the mask uses the 3-D
  Euler-characteristic term, as FSL `ptoz -g` does:
  R (4 ln 2)^{3/2} (2π)^{-2} (z² − 1) e^{−z²/2}. It is one-sided positive,
  then Bonferroni ×3.
- **Whole brain.** `cluster` at z > 3.09 on zstat and on −zstat, with GRF cluster
  p from DLH and volume. Two-sided by Bonferroni over the two tails (cluster
  p × 2).
- **Weights.** Each participant's share of the total weight
  1 / (varcope_i + σ²_g), using ROI means of varcope and of FLAME's
  `mean_random_effects_var1`.
- **Robustness.** REML random-effects meta-analysis of each participant's ROI
  mean effect, with the ROI mean variance as its sampling variance.

## Pipeline fix before any result (2026-11-26, 12:50)

The first full first-level pass wrote every map through nilearn's
`inverse_transform`, which copies the **uint8** header of the group mask. Every
saved effect and variance map was therefore quantised to 8 bits, and small
variances were stored as exactly 0. The symptom was found at the first group
step: FLAME estimated only 15,314 of 219,413 voxels, and the permutation tests
crashed on an unrelated return-type error. No group statistic was read: the H3
table was never printed, and the other tests did not run.

**Fix.**
- All images are now written as float32 (`save_float32`).
- The first level refuses to save any variance map with a non-positive voxel.
- The quantised outputs were deleted and every first level refitted.

On sub-082 run-1 the saved maps now equal `FirstLevelModel`'s to within
2e-7 / 5e-7. The 0.019 residual quoted above was this quantisation.
The permutation ROI tests pass `threshold=0.002` only so that nilearn returns
the t map alongside the voxel max-t FWE map.
