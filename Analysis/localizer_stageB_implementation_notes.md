# Localizer stage B implementation notes (written before any stage-B statistic was computed)

Written 2026-11-27, after stage A was committed (87a8be9) and **before any
participant-specific ROI (fROI) was formed, any reliability-gate index was
computed, or any breathing-run contrast was refitted or extracted at 4 mm**.
Stage B covers the participant-specific ROIs, the reliability gate L4, the
transfer tests T1 and T3b, secondary T3, the exploratory C2 and C3 fROI rules,
exploratory T2, and the registered stage-B sensitivity analyses. Connectivity,
gPPI, T4, T5, decoding, RSA and PLV are stage C and are not run here.

The plan of record is `Preregistration/intero2023_localizer_prereg.md` (v1.0,
frozen at d29136b), §7–§9, with the OSF form text
`Preregistration/osf_localizer_secondary_prereg.md` (Indices 3–4, Analyses,
Inference criteria). These notes record the choices the registration leaves to
the analyst. Readings of ambiguous text that could be argued to change a test are
also logged in `Preregistration/localizer_deviations.md`, before the affected
analysis is run. The breathing-run pipeline and the stage-A localizer scripts are
imported, never modified.

## Scripts
| Script | Role |
|---|---|
| `Scripts/intero2023_loc_breath4mm.py` | Breathing-run refit at 4 mm and unsmoothed, one participant per call; `check` mode re-runs one participant at 6 mm and compares with the swc7p maps |
| `Scripts/intero2023_loc_froi.py` | fROI definition, L4 gate, extraction, T1, T3b, T3, T2, exploratory rules, sensitivity → `Analysis/results/LOC_*` |

Imaging outputs (refitted maps, fROI masks) go to
`derivatives/intero2023_glm/firstlevel/<sub>/` and `derivatives/intero2023_glm/froi/`,
never into the repository. Tables in `Analysis/results/` hold numbers only
(counts, means, statistics), no voxel coordinates of participant ROIs.

## Contrast mapping (the naming trap)
The breathing-run files use the swc7p contrast numbering
(`Scripts/intero2023_firstlevel.py`, `intero2023_group.CON`). It is **not** the
localizer numbering.

| Breathing-run file `<sub>_var-base_con-<X>` | Meaning | swc7p test | Used here for |
|---|---|---|---|
| `C1` (per run, then inverse-variance FFX over runs) | paced − free: unweighted mean of the paced conditions present in the run (`paced_aware`, `paced_unaware`, `paced_nochange`) minus `free`; Model A | H1, H6 | T3, T3b |
| `C2joint` (joint two-run Model A) | aware − unaware (`paced_aware` − `paced_unaware`) | H3 | T1 |
| `C3` (Model B per run, FFX over runs whose confidence varies) | confidence parametric modulator | H2 | T2 (exploratory) |
| `E3changejoint` (joint Model A) | change (mean of `paced_aware`, `paced_unaware`, `paced_wrongdir` present) − `paced_nochange` | E3a/E6 | T2 (exploratory) |
| `C2` (per run), `C5` | per-run aware − unaware; step − ramp | H3 per-run check; H4 | not used |

Localizer files (`<sub>_task-localizer_loc-<joint|run1|run2>_sm-4_var-base_con-<c>_<kind>`)
use the localizer numbering: C1 = Intero − Extero, C2 = Intero − Self,
C3 = min(t_C1, t_C2), C5 = Match − Intero. In these notes and the results,
"breathing H3 contrast" or "aware − unaware" always means breathing `C2joint`,
"paced − free" always means breathing `C1`, and "C1/C2/C3/C5" without a qualifier
means the localizer contrast.

## Participant sets
| Set | Rule | N |
|---|---|---|
| fROI-eligible | `compliance.csv` `froi_eligible` (two usable runs, ≥ 3 of 4 compliant Intero blocks, ≥ 1 per run) | 44 |
| LOSO group-map source | localizer Intero set (`meets_75pct`), each target left out | 48 (47 per map) |
| **T1** | fROI ∩ swc7p H3 set (`Data/h3_sets.csv` `in_H3`) | **40** |
| T3 | fROI ∩ swc7p H1 set (all 47 breathing participants) | 41 |
| **T3b** | T3 ∩ swc7p H6 set (`maia_ieat_total` available) | **40** |
| T2 confidence | fROI ∩ swc7p H2 set (44: not sub-066, sub-070, sub-114) | 39 |
| T2 change − No Change | fROI ∩ E3a set (47) | 41 |
| No-v1 sensitivity | T1 / T3b minus protocol-v1 participants (P51–59) | 34 / 34 |

- The four fROI-eligible participants outside T1: sub-056, sub-080 and sub-085
  have no usable breathing run; sub-070 never answered "No Change" on a change
  trial. sub-065 (no MAIA) is in T3 but not T3b.
- Every T1 participant has two gate-passing breathing runs.
- The registered Ns (T1 40; T3b "about 40") are met exactly.
- MAIA is `maia_ieat_total` (IEAT scoring, 0–5), demeaned over the analysed sample.

## fROI definition (plan §7; p5ztc Indices 3)
- **Constraint parcels:** `caren_AInsL` (588 voxels), `caren_AInsR` (527),
  `caren_dACC` (661), `caren_PIns` (285) and `ieat_ACCsphere` (528). All lie fully
  inside the group brain mask.
- **Size:** n = max(ceil(0.10 × parcel voxels), 50): AInsL 59, AInsR 53, dACC 67,
  PIns 50 (10% would be 29, so the minimum applies), sphere 53. No parcel has fewer
  than 50 voxels, so the per-participant fallback for an unformable fROI (p5ztc
  Data inclusion 10) never triggers.
- **Ranking map:** the participant's two-run (`loc-joint`) 4 mm t map. Voxels are
  sorted by t, descending, with a stable sort in mask-voxel order, so ties (none
  expected in float32 t) break deterministically.
- **Primary rule:** the top n voxels by C1 t (Intero > Extero).
  - Group C1 is Intero < Extero almost everywhere in these parcels (stage A, L1).
    The rule is applied as registered, so the top 10% are each person's
    least-deactivated (or breath-preferring) voxels. The distribution of the fROI
    C1 values (mean t and mean effect in the fROI, and the share of participants
    whose fROI mean C1 is > 0) is reported per parcel. It is descriptive and
    circular, because the same data select and describe the fROI.
- **Specificity controls:** the complement (parcel minus fROI) and the bottom n
  voxels by C1 t (Extero > Intero).
- **Exploratory rules:** the top n by the C2 t (Intero > Self) and by the C3 map
  (voxelwise min of the C1 and C2 t of the joint model, as saved in stage A).
- **Secondary rule (T3):** the bottom n by the C5 t (Match < Intero, "negative
  direction").
- **Size sensitivity:** 5% and 20%, with the same 50-voxel minimum, because only
  the percentage is registered as varying. Consequence, stated now: at 5% every
  parcel's fROI is the 50-voxel minimum (AInsL 50, AInsR 50, dACC 50, PIns 50,
  sphere 50). At 20%: 118, 106, 133, 57 and 106.
- fROI masks are written as NIfTI to `derivatives/intero2023_glm/froi/<sub>/` only.

## L4 reliability gate (plan §7; p5ztc Indices 4)
Computed on localizer data only, for each parcel, in the 44 fROI-eligible
participants, from the per-run 4 mm fits (`loc-run1`, `loc-run2`).
1. **Cross-run effect.** For direction j → k, the fROI is defined from run j's C1
   t map with the primary rule and size, and the run-k C1 **effect** (percent
   signal change) is averaged over it. The participant's index is the mean of
   1 → 2 and 2 → 1.
   - **Group test:** a one-sample t-test of the index against 0, one-sided
     (C1 > 0), p < .05. The OSF power table registers this test as one-sided.
   - **Sign agreement:** the share of participants whose index is > 0, which must
     be ≥ .70. The share of participants whose two directions agree in sign is
     reported alongside.
   - (1) holds if both parts hold.
2. **Dice.** 2|A∩B| / (|A| + |B|) between the run-1 and run-2 fROIs. (2) holds if
   the median over participants is ≥ .20. Chance is the fROI share of the parcel
   (about .10; .18 for PIns, where 50 of 285 voxels are selected). That chance
   level is reported.
3. **Spatial similarity.** Pearson r between C1 t maps over the parcel's voxels.
   - Within = the participant's run 1 vs their own run 2.
   - Between = the mean, over the other 43 participants, of r(own run 1, other's run 2).
   - Each r is Fisher-z transformed before averaging and testing.
   - **Test:** a paired t-test of z(within) − mean z(between), one-sided
     (within > between), p < .05. The Wilcoxon signed-rank p is reported alongside.
- **Gate outcome.** The parcel passes L4 if (1), (2) and (3) all hold.
  - **If (3) fails,** the LOSO group fallback becomes the primary fROI for that
    parcel, as registered.
  - **If (1) or (2) fails but (3) holds,** the individual fROI stays primary, as
    the registration triggers the fallback on (3) only. L4 is then "not met" for
    that parcel. The transfer tests are still run and reported, with that caveat
    attached. The registration does not make transfer conditional on L4.
- The gate is also applied to the exploratory C2 and C3 rules ("the same gate"),
  each with its own contrast.
  - **C2 rule:** the C2 effect, C2 fROIs and C2 t maps.
  - **C3 rule:** fROIs from the per-run C3 (min-t) maps. Criterion (1) requires
    both the cross-run C1 effect and the cross-run C2 effect to pass, because C3
    claims selectivity on both. Criterion (3) uses the C3 maps.
- The 5% and 20% size variants are not gated; they are sensitivity analyses of the
  transfer tests.

## LOSO group fallback (plan §7)
- **Source:** for each target participant, the group map is the OLS one-sample t
  of the 4 mm two-run C1 **effect** over the localizer Intero set (48) with the
  target removed (47 maps).
  - The Intero set is the registered set for C1 group inference (L1).
  - Leaving the target out keeps the target's own data out of the selection.
- **fROI:** the top n voxels of that map within the parcel, with the same size rule.
- The same construction is used for the exploratory rules:
  - **C2 rule:** the LOSO t of the C2 effect.
  - **C3 rule:** voxelwise min of the LOSO C1 and C2 t.
- The fallback fROI is computed for every parcel and reported for T1 and T3b
  alongside the primary, whatever the gate says. It replaces the individual fROI
  as primary only where criterion (3) fails.

## Breathing-run refit at 4 mm, and unsmoothed (plan §3; p5ztc Transformations 3)
`Scripts/intero2023_loc_breath4mm.py <sub>` refits the registered swc7p models.
The only change is the smoothing kernel. Everything else is `intero2023_firstlevel.main`:
- the same run list (`g.run_list()`: 90 runs, 47 participants);
- the same designs (`q.model_a`, `q.model_b` via `g.read_run_physio(..., [])`,
  i.e. the base nuisance set): SPM + derivative HRF, fMRIPrep 128 s cosines in
  the confounds, `drift_model=None`;
- the same masker (`ROIs/group_brain_mask.nii`, float32), `mean_scaling` per run,
  and `fit` (`run_glm` AR(1) + `compute_contrast`);
- per-run Model A gives C1 (`g.c1_spec`) and Model B gives C3 (runs whose
  confidence varies, `g.conf_varies`), each combined over runs with
  `intero2023_firstlevel.ffx`;
- the joint Model A (`hj.build_joint_design`, runs PSC-scaled then concatenated)
  gives C2joint (aware − unaware) and E3changejoint (the E3a spec of
  `intero2023_explore_awareness.e3subject`);
- only the base variant is refitted (see the deviations log for the physiological
  variants).

**Kernels:**
- `smoothing_fwhm = 4`;
- `smoothing_fwhm = None` (unsmoothed sensitivity, labelled `sm-0`).

**Outputs:** `<sub>_sm-<4|0>_var-base_con-<C1|C3|C2joint|E3changejoint>_{effect,variance}.nii.gz`
and `<sub>_breath_sm4_sm0.json` (runs, dofs, contrasts).
- These names never collide with the swc7p 6 mm files, which have no `sm-` entity.
- A variance map with a non-positive voxel is an error, as in swc7p.

**Check before the batch.** `check sub-060` runs the same code at 6 mm into the
scratch directory and compares C1, C3, C2joint and E3changejoint with the swc7p
files. The result goes to `Analysis/localizer_stageB_refit_check.json`.

**Run:** all 47 breathing participants, `xargs -P 8`, `OMP_NUM_THREADS=2`,
absolute script paths, `~/miniforge3/envs/neuro/bin/python`. Logs go to
`derivatives/intero2023_glm/logs/locB_breath_<sub>.log`. After the batch every log
is checked for "done" and every expected file for existence.

## Extraction (p5ztc Indices 3)
For each participant, rule, size, parcel and breathing contrast:
- **fROI index:** the unweighted mean of the contrast effect over the fROI voxels.
- **Complement index:** the same mean over the rest of the parcel.
- **Parcel index:** the mean over the whole parcel (the anatomical-mask control).
- **Bottom-10% index:** the same mean over the bottom-10% voxels (specificity
  control, descriptive).
- **Variance:** each index's variance (FLAME input) is the mean of the voxel
  variances over the same voxels. This is the H3 convention (`flame_tests`,
  `reml_meta`).
- **Smoothing:** extraction is from the 4 mm refit (primary), the unsmoothed refit
  and the swc7p 6 mm maps (sensitivities).

## T1: aware > unaware in the C1 fROIs (plan §8–§9; p5ztc Inference criteria)
- **Masks:** AInsL, AInsR and dACC fROIs (primary rule, or the fallback where the
  gate's (3) fails), N 40, breathing `C2joint` at 4 mm.
- **FLAME-weighted, as H3.** FSL `flameo --runmode=flame1` with a one-group design.
  - Cope, varcope and dof enter exactly as in H3 (dof = the joint-model residual
    dof from the refit json).
  - Every "voxel" of the input image is one index: the fROI, complement and parcel
    of each mask, for every rule and variant in the same call. flameo treats voxels
    independently, so each index gets its own FLAME 1 fit.
  - **p** = one-sided upper tail of FLAME's z, × 3 (Bonferroni over the three masks).
  - GRF does not apply to a single value.
  - Reported alongside: the FLAME weight shares, the REML meta-analysis
    (`G.reml_meta`), the unweighted mean, dz and one-sample t.
- **Paired fROI − complement test.** d_i = fROI index − complement index, in the
  same 40.
  - **Statistic:** dz = mean(d) / SD(d), with a 95% percentile bootstrap CI
    (10,000 resamples of participants, seed 20261127; p5ztc Indices 3 and
    Transformations 14).
  - **Secondary:** the difference dz(fROI) − dz(complement), with its own
    bootstrap CI from the same resamples (plan §8's "dz difference" wording).
  - The paired t, one-sided, is reported descriptively.
- **Verdict per mask:** T1 is **supported** in a mask if:
  - the FLAME one-sided p × 3 < .05, **and**
  - the 95% CI of dz(fROI − complement) lies above 0.

  T1 is supported overall if it is supported in at least one mask, the H3 rule.
  "More sensitive than the within-run analyses" is claimed only for a mask where
  T1 is supported.
- **Bayes factors.** `jzs_bf10` (two-sided) and `bf_plus0` (directional, effect
  > 0) at r = .707 on the unweighted fROI means.
  - They are reported next to E5 (BF01 5.1 / 2.2 / 5.8; BF0+ 3.4 / 1.2 / 8.1 in
    the 45; 6 mm parcel means), read from `EXPLORATORY_E5_bayes_factors.csv`, and
    next to the same BFs for the 4 mm parcel means in the 40.
  - BF01 > 3 in an fROI is read as moderate evidence against the awareness effect
    there.
- **Anatomical-mask comparison (registered: same participants, same estimator).**
  The parcel index and the complement index go through the same FLAME,
  BF and dz computation in the same 40. The registered H3 voxelwise result (6 mm,
  N 45) is quoted from `H3_var-base_set-all_roi.csv` for context, never as the
  comparison.

## T3b: MAIA × paced − free > 0 in the ACC-sphere fROI (one-sided)
- **Data:** the sphere's primary C1 fROI (or the fallback, if the gate's (3)
  fails), breathing `C1` at 4 mm, N 40.
- **Model, as H6:** OLS of the fROI index on demeaned MAIA with an intercept. The
  test statistic is the MAIA slope t.
- **p:** the one-sided (positive) permutation p over 10,000 permutations of MAIA,
  seed 20261126. This mirrors H6's engine, MAIA permutation with an intercept, on a
  single value. The parametric one-sided p is reported alongside.
- **Verdict:** supported if the permutation p < .05.
- Reported with it: Pearson r, the JZS regression BF10 (`jzs_reg_bf10`, r = .707),
  the same model for the complement and the whole sphere, and fROI − complement
  (slope difference, with a bootstrap CI).

## T3 (secondary): paced − free in the C5 fROIs
Read in the deviations log. C5 fROIs (bottom 10% by C5 t) are formed in the five
registered constraint parcels.
- **Test, as H1:** paced − free < 0, one-sided, N 41.
  - Statistic: the one-sample t of the fROI index.
  - p: from 10,000 sign-flips (seed 20261126), Bonferroni × 5 over the parcels.
- **H6 analogue:** MAIA × paced − free in the sphere C5 fROI, as T3b, N 40.
- Reported for each: fROI − complement (dz and bootstrap CI), the parcel index and
  the BFs.
- This is a secondary analysis; it cannot support or overturn T3b.

## Exploratory, labelled as such
- **C2 and C3 fROI rules** (registered exploratory). They pass through the same
  L4 gate (above) and the same transfer tests: T1 with its paired test and BFs,
  and T3b.
- **T2** (registered exploratory) uses the primary C1 fROIs in AInsL, AInsR and
  dACC.
  - Contrasts: breathing `C3` (confidence; N 39) and `E3changejoint` (change − No
    Change; N 41), at 4 mm.
  - FLAME on the index as for T1, but two-sided. Uncorrected p and p × 3 are both
    reported.
  - Also reported: dz, BF10 and fROI − complement.

## Sensitivity analyses (stage B)
The analyses below apply to T1 (each mask, FLAME p × 3, paired dz and CI, BF) and
to T3b, with the gate-determined primary fROI. Each is reported; none replaces
the primary.
1. **Unsmoothed extraction:** the `sm-0` refit (registered).
2. **fROI size 5% and 20%** (registered), with the 50-voxel minimum, as above.
3. **6 mm extraction:** the swc7p maps themselves. Added at the PI's instruction;
   logged in the deviations.
4. **Excluding protocol v1** (registered; N 34 / 34).
5. **The individual fROI and the LOSO fallback, side by side,** whichever is
   primary.

## Bonferroni, tails and floors
- T1: one-sided, ×3. T3b: one-sided, ×1. T3: one-sided, ×5. T2: two-sided, ×3.
- L4: each criterion at α = .05, no correction (as registered).
- Permutation p has a floor of 1e-4.

## Outputs
- `Analysis/results/LOC_froi_definition.csv`: per participant × rule × size × parcel,
  the number of voxels, the mean C1 t and C1 effect in the fROI and complement.
  Plus `LOC_froi_definition_summary.csv`.
- `Analysis/results/LOC_L4_gate.csv` (per parcel × rule) and `LOC_L4_gate_subjects.csv`.
- `Analysis/results/LOC_froi_extract.csv`: long-format indices.
- `Analysis/results/LOC_T1.csv`, `LOC_T1_weights.csv`, `LOC_T3b.csv`, `LOC_T3.csv`,
  `LOC_T2.csv`, `LOC_exploratory_rules.csv`, `LOC_stageB_sensitivity.csv`.
- `Analysis/LOCALIZER_RESULTS.md`, "Stage B" section.
