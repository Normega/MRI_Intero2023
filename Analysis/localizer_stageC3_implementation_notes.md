# Localizer stage C3 implementation notes: pattern expression, decoding, RSA, gradients, respiratory phase (written before any stage-C3 statistic was computed)

Written 2026-11-27, after stage C2 was committed (cec9af3) and **before any
stage-C3 value was computed**: no LSS beta, index, classifier accuracy,
distance, slope or phase effect existed when these notes were committed.
Stage C3 is **exploratory** throughout (plan §9: "RSA, decoding and gradients",
"T4", "phase-locked responses"; §10 ranks the techniques). Nothing here tests a
registered hypothesis, and nothing changes an earlier result.

Plan of record: `Preregistration/intero2023_localizer_prereg.md` v1.0 (d29136b),
§3, §7, §8 (T4), §9, §10, §11. Readings of ambiguous text and analyst choices
are logged in `Preregistration/localizer_deviations.md` (entry "Stage C3")
before running. The registered pipeline (`intero2023_glm_common`,
`intero2023_twomodel_qc`, `intero2023_firstlevel.fit`, `intero2023_loc_common`)
is imported, never modified.

**What was run before these notes (no statistic, no BOLD):**
- LSS design check on sub-060 (`Analysis/localizer_stageC3_t4_design_check.json`):
  all 28 LSS designs are full rank, and the target column plus the remaining
  same-condition column reproduce the Model A column exactly (max |Δ| 0).
- Respiratory-phase design check, design matrices and belt only
  (`Analysis/results/EXPLORATORY_LOC_C3_respphase_design_check.csv`): see §6.
- Counts of epochs, maps and runs (below).

## Scripts
| Script | Role |
|---|---|
| `Scripts/intero2023_loc_t4.py` | T4: LSS betas (`lss`), design check (`verify`), indices (`index`), reliability (`reliability`) |
| `Scripts/intero2023_loc_t4_mixed.R` | T4 mixed models (lme4 / lmerTest / emmeans) |
| `Scripts/intero2023_loc_physio_relations.py` | §10.2 between-participant relations |
| `Scripts/intero2023_loc_decoding.py` | LOSO decoding and cross-decoding |
| `Scripts/intero2023_loc_rsa.py` | Cross-run crossnobis RSA (`subject`, `group`) |
| `Scripts/intero2023_loc_gradient.py` | Insula posterior → anterior C1 gradient |
| `Scripts/intero2023_loc_respphase.py` | Inhale vs exhale by condition (`segments`, `check`, `subject`, `group`) |
| `Scripts/intero2023_loc_c3_figures.py` | Figures |

Imaging intermediates go to `derivatives/intero2023_glm/loc_c3/` only. The
repository gets tables (`Analysis/results/EXPLORATORY_LOC_C3_*.csv`) and figures.

## 1. T4: per-epoch expression of the localizer C1 pattern (plan §8, §10.1)
- **Set.** Localizer Intero set (48; each has their own C1 map) ∩ breathing-run
  set (swc7p, 90 runs / 47): **44 participants**. Out: sub-068, sub-074, sub-104
  (breathing but no localizer Intero contrast) and the localizer-only participants.
- **Epochs.** Every paced epoch (aware 291, unaware 142, wrong direction 77,
  No Change 85; 595) and every free epoch (594) of the 44 participants' swc7p runs.
  "About 14" in the plan is the paced count per participant (median 14 over
  two runs).
- **LSS** (Mumford et al. 2012, by-condition form). For each target epoch, the
  registered swc7p Model A design of that run (`q.model_a`: the seven condition
  boxcars, SPM + derivative, nuisance sheet, no drift) with the target moved to
  its own `target` regressor (+ derivative); the other epochs of its condition
  stay in their regressor. AR(1) through `intero2023_firstlevel.fit`. The
  target's canonical beta is kept (PSC).
  - **Data:** breathing BOLD at **4 mm** (plan §3: pattern analyses), PSC per run,
    voxels of CAREN **SAL ∪ SMN** (`net_SAL` 4,156 + `net_SMN` 9,033 = 13,189
    voxels, no overlap). Smoothing is of the whole image before masking, so the
    values equal group-mask 4 mm data in those voxels.
- **Templates** (same voxels):
  - **own** (primary): the participant's joint two-run 4 mm localizer C1 effect
    (Intero − Extero, PSC). For one-run participants, that run's model.
  - **loso** (registered alternative): OLS one-sample t of the joint 4 mm C1
    effect over the Intero set (48) without the participant, as the stage-B
    fallback map.
  - **run1 / run2**: the participant's per-run 4 mm C1 effect (reliability only).
- **Index** (the plan says "dot product (or correlation)"; both are computed):
  - **dot (primary):** b · w / (‖w‖ √V), the projection of the epoch beta pattern b
    on the unit-length template, per voxel. If b = w it equals RMS(w). Templates
    are normalised so that participants with stronger C1 maps do not get larger
    indices for the same epoch pattern.
  - **corr (secondary):** Pearson r(b, w) over voxels: the spatial shape only, with
    both means removed.
  - Also per network (SAL alone, SMN alone), descriptive.
- **How the sign reads.** C1 is mostly negative (Intero < Extero, stage A). A
  **positive dot index** means the epoch is lower where Intero is lower than
  Extero: it looks more like Intero than Extero, **including the global
  deactivation**. Because the breathing runs have no clean implicit baseline
  (events cover about 94% of each run), the mean level of b is defined only
  relative to the leftover baseline and the run intercept, so the dot index mixes
  pattern and level; only within-participant (and within-run) comparisons of
  epochs are interpretable. A **positive corr index** means the epoch's spatial
  shape resembles C1 regardless of its mean level.
- **Mixed model (registered form, plan §8):** `index ~ confidence + outcome +
  (1 | participant)`, paced epochs, `outcome` = aware / unaware / No Change as the
  swc7p events define them (`paced_aware`, `paced_unaware` = reported "No Change"
  on a change epoch, `paced_nochange` = no-change epoch, any report). Confidence
  is the raw 1–7 rating. **Wrong-direction epochs** are not one of the three
  outcomes and are left out (sensitivity: as a fourth level). lme4 REML,
  lmerTest Satterthwaite df, Wald 95% CI; outcome contrasts aware − unaware,
  aware − No Change, unaware − No Change from emmeans (Satterthwaite), and the
  omnibus F for outcome. Fitted for dot and corr with the own and loso templates
  (primary: dot, own), and for SAL and SMN separately (descriptive).
- **Sensitivity (dot and corr, own template):** run intercepts nested in
  participant (`+ (1 | participant:run)`); confidence split into within- and
  between-participant parts; wrong direction as a fourth outcome.
- **Added descriptive check:** paced vs free epochs, `index ~ epoch type +
  (1 | participant) + (1 | participant:run)`. H1 found paced < free activity,
  so with a negative template a higher paced dot index would be expected from the
  level alone; corr separates the shape.
- **Reliability** (`EXPLORATORY_LOC_C3_T4_reliability.csv`):
  1. template split-half: the index computed with the run-1 vs the run-2
     localizer C1 map, across paced epochs within participant (participant means
     removed), Spearman–Brown corrected; and the between-participant r of the
     participant means;
  2. breathing run 1 vs run 2: between-participant r of the participant mean
     paced index (Spearman–Brown), and of the participant aware − unaware
     difference where ≥ 10 participants have both cells in both runs (otherwise
     reported as not estimable);
  3. the ICC (participant-intercept share) of the index in the paced epochs.
  These serve as the noise ceiling: an outcome effect cannot be detected if the
  index's within-participant variation is mostly noise.
- No correction across indices; dot / own is primary and everything else is
  labelled.

## 2. Physiology-informed between-participant models (plan §10.2; brief)
Stage C1 already related tracking PLV to C1 in the fROIs and masks, Match
adherence to breathing-run adherence, and adherence to H1. What is new after
C2 and T4 is a small family (BH-FDR over 4), using C1's measures and its
`corr_row` (Spearman, permutation p, bootstrap CI, Pearson, JZS BF10):
1. tracking PLV × C1 SAL–DMN background-FC difference (Intero − Extero, base,
   median, 4 s shift; the C2 primary finding);
2. tracking PLV × participant mean paced-epoch T4 index (dot, own);
3. tracking PLV × participant aware − unaware T4 index (participants with ≥ 1
   epoch of each);
4. breathing-run adherence (`adhere_r`) × participant mean paced-epoch T4 index.
Anything else is declared covered by C1.

## 3. LOSO decoding of Intero vs Extero (plan §10.3)
- **Samples:** the per-run 4 mm `intero` and `extero` (each vs Rest) effect maps
  of the Intero set (48 participants, 92 runs, 184 maps; every run has both).
- **Masks:** whole brain (group mask, 219,413 voxels); CAREN SAL, SMN, DMN; insula
  = CAREN AInsL ∪ AInsR ∪ PIns.
- **Features:** *raw* PSC (primary), and *pattern*: each map centred and divided
  by its SD over the mask's voxels, which removes the map's mean level (the
  global deactivation) and scale.
- **Classifier:** linear SVM (`sklearn.svm.SVC`, precomputed linear kernel), C = 1,
  no tuning (with far more voxels than maps the solution is effectively the
  hard-margin one); leave one participant out (48 folds).
- **Score:** balanced accuracy over all held-out maps (each participant
  contributes equal numbers of each class); a participant-bootstrap 95% CI
  (10,000).
- **Permutation p:** 5,000 permutations, seed 20261129. In each, the two labels
  of every run's Intero/Extero pair are swapped with probability .5 (training
  and test alike), and the whole LOSO is rerun. p = (1 + #null ≥ observed) /
  5,001.
- **Cross-decoding (optional in the plan; done):** the classifier trained on the
  other participants' Intero vs Extero maps is applied to the held-out
  participant's **Match** maps (scored as the Intero class) and Extero maps.
  Balanced accuracy and the share of Match maps called Intero. Permutation as
  above, permuting the training labels only. Descriptive addition: the share of
  the held-out **Self** maps called Intero (no p), to show whether the
  classifier separates "breath" or just "not Extero".
- 5 masks × 2 feature sets = 10 tests; the count surviving Bonferroni ×10 is
  stated.

## 4. Cross-run crossnobis RSA (plan §10.4; descriptive)
- **Set:** fROI-eligible participants (44): two usable runs, all four modelled
  conditions estimable in each.
- **Patterns:** the registered per-run 4 mm localizer model (`L.run_designs`,
  base, 200 s DCT, AR(1)) is refitted in the ROI voxels to obtain the betas of
  `intero`, `extero`, `match`, `self` and the whitened residuals. The betas are
  checked against the saved per-run 4 mm effect maps. **Rest** is the implicit
  baseline, so its pattern is zero and all **five** conditions enter.
- **Noise normalisation:** per run, the covariance of the whitened residuals with
  Ledoit–Wolf shrinkage (`sklearn.covariance.ledoit_wolf_shrinkage`), applied as
  Σ^(−1/2) (computed exactly through the residual SVD). Crossnobis distance
  d(i, j) = (b_i¹ − b_j¹)ᵀ(b_i² − b_j²) / P over the two runs. It is unbiased:
  its expectation is 0 when two conditions do not differ.
- **ROIs:** the five a-priori parcels (AInsL, AInsR, dACC, PIns, IEAT sphere);
  the networks SAL, SMN, DMN, VIS; the gate-primary fROIs (individual C1 in AInsL,
  dACC, PIns; LOSO in AInsR, sphere). The individual fROIs were selected on the
  two-run C1 (both runs), so their Intero–Extero distance is biased upwards; the
  a-priori masks and the LOSO fROIs are not.
- **Model RDMs** (conditions I, E, M, S, R; distance 0 / 1):
  - *breath*: attention to or control of the breath, {I, M} vs {E, S, R};
  - *visual*: moving circle {I, E, M} vs words {S} vs fixation {R};
  - *self*: S vs every other condition.
  They are collinear (visual and self share every pair that involves S); their
  correlations and VIFs are reported.
- **Statistics (descriptive):** group mean RDM; each distance vs 0 (one-sample t);
  per-participant Pearson and Spearman r between the 10 distances and each model,
  group mean, one-sample t and Wilcoxon. No correction; the count with p < .05 is
  stated.

## 5. Insula posterior → anterior gradient of C1 (plan §10.5)
- **Set:** localizer Intero set (48). **Map:** the joint 4 mm C1 effect (PSC).
- **Region:** CAREN insula, AInsL ∪ AInsR ∪ PIns (1,400 voxels; PIns is split by
  the sign of x for the per-hemisphere slopes).
- **Slope:** per participant, OLS of the voxel C1 effect on MNI y, in PSC per
  10 mm; pooled over hemispheres (primary) and per hemisphere. A positive slope
  means C1 is larger (less negative) anteriorly.
- **Group test:** one-sample t (two-sided), dz, 95% CI, and a sign-flip
  permutation p (10,000, seed 20261126).
- **Sensitivity:** x and z as covariates; Harvard–Oxford insula (`ho_AInsL` ∪
  `ho_AInsR` ∪ `ho_PIns`); the 6 mm C1 map. Descriptive: the same slope for
  Intero vs Rest and Extero vs Rest, to show which condition carries it; a y-bin
  profile for the figure.
- **Principal cortical gradient:** no gradient map (e.g. Margulies et al. 2016)
  is available locally (a file-system search found none), and none is
  downloaded. The projection is **not done**.

## 6. Inhale- vs exhale-locked responses by condition (plan §10.6; descriptive)
- **Runs:** usable localizer runs whose belt is good or degraded, with trigger
  timing exact or good, pulses = volumes, and a resolved participant polarity
  (stage-C1 rules, `intero2023_localizer_qc.physio_segments`): **69 runs, 40
  participants**. Intero rows are used only for Intero-set participants.
- **Breaths and phase:** the step-1 `breath_table` (inhale onset = trough, exhale
  onset = peak after polarity) and the Hilbert phase of the polarity-corrected
  band-passed belt (0 at the end of inhalation), as stage C1. The first 3 s of
  each block are dropped. Volume k is sampled at k·TR + 0.958 s (fMRIPrep
  slice-timing reference, `StartTime`).
- **Two candidate models**, each added to the registered per-run design (6 mm,
  mask-mean PSC series of AInsL, AInsR, dACC, PIns, COper, SM and the IEAT sphere;
  AR(1)):
  - *onset model*: impulses at the inhale and at the exhale onsets within each
    condition's task window, SPM canonical HRF; contrast inhale − exhale;
  - *phase model*: cos φ and sin φ at each volume within each condition's task
    window, not convolved (first-order RETROICOR form, per condition).
    Inhale − exhale = mean over the inhale half-cycle minus the exhale half
    = −(4/π)·b_sin. Amplitude = √(b_cos² + b_sin²).
- **Feasibility rule and check (design only).** The onset model is fitted only if
  the median VIF of every inhale and exhale column is ≤ 10. The threshold was
  written into the script before the check ran; the script was not committed
  before the check, so this cannot be verified from the history. Check result
  (69 runs): r(inhale, exhale) median .91–.92 in every condition, VIF 9.2–29.5.
  At about 12 breaths/min the inhale − exhale difference is a 0.2 Hz ripple that
  the HRF all but removes. **The onset model is therefore not fitted.** The phase
  model has VIF 1.2–1.6 and is the analysis run.
- **Participant value:** inverse-variance fixed effect over runs.
- **Tests (descriptive):** inhale − exhale per condition vs 0; condition
  differences in inhale − exhale (Intero − Extero, Match − Extero, Intero − Rest,
  Match − Intero); amplitude differences (Intero − Extero, Match − Extero,
  Intero − Rest). One-sample t; BH-FDR within each family over masks × tests.
- **How to read it.** Unconvolved respiratory-phase BOLD is dominated by
  physiological artefact (B0 changes from chest motion, head motion, CSF
  pulsation), not neural responses, and it grows with breathing depth, which
  differs by condition (Intero and Match are deeper). The amplitude is also
  positively biased by noise. The results are therefore described, not
  interpreted as neural inhale- or exhale-locked responses.
