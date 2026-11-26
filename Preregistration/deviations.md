# Deviations from the OSF registration

**Registration of record:** OSF Preregistration *Neural correlates of awareness
of respiratory change during paced breathing*, https://osf.io/swc7p, registered
**2026-08-20 12:23 UTC** (08:23 EDT), public, no embargo. **Updated
2026-11-26** (OSF schema response `6b0863b9dab3517393f851af`, submitted 16:59
UTC and approved the same day), carrying D1–D26 below with the justification
text from `osf_template_answers.md`. The
registered form text is the repository's `Preregistration/osf_template_answers.md`
at commit `7472481` (2026-08-20 08:09 EDT); two fields (foreknowledge and data
collection procedures) were edited on the form before submission. The companion
plan was v0.1 (`Preregistration/intero2023_imaging_prereg.md` at commit
`aa0752b`). The form says the plan is attached, but no file is attached to the
registration.

**What follows the registration.** Plan versions v0.2 (2026-08-24) to v0.8
(2026-11-25) were all written after the registration. Those documents call
their changes "pre-freeze"; that wording is wrong with respect to swc7p (see the
correction note in each document). Every change below is therefore a deviation
from the registration. They were filed together as the 2026-11-26 OSF update.

**Imaging results seen before any deviation: none.** No first-level model has
been fitted at any point, no contrast or statistical image exists, and no BOLD
voxel intensity has been related to any task regressor. The work behind these
deviations used:
- timing sheets, behaviour and physiology;
- fMRIPrep confound columns (framewise displacement, the tissue signals, and the
  global signal);
- fMRIPrep brain masks, for the ROI grid and coverage.

One item used a BOLD-derived quantity in deciding a rule. The PMU fallback gate
(D11) was revised after the fallback runs' RVT–global-signal correlations had
been seen. That fallback feeds robustness analysis S2b only.

Columns: **Registered** is swc7p; **Now** is plan v0.8; **When** is the plan
version and date the change was decided. Rationale in full is in
`methods_decisions.md` (decision numbers given).

## Sample and exclusions

| # | Registered (swc7p) | Now (v0.8) | When | Why |
|---|---|---|---|---|
| D1 | sub-082 and sub-104 excluded for shifted run labels; motion gate fails 11 of 95 runs; 43 analysable participants | sub-065, sub-076, sub-082 and sub-104 reinstated; sub-051 (no MRI) and sub-056 (truncated run-1, run-2 fails the gate) named; gate fails 12 of 103 runs; **90 runs, 47 participants** | v0.2, 2026-08-20 / 08-24 | The registered inventory was built from an obsolete copy of the data; in the dataset of record the runs map correctly (`VERIFICATION_REPORT.md`; decision 4) |
| D2 | One-run participants enter with that run's estimate | Stated as an inclusion rule: a participant enters with at least one gate-passing breathing run; S6 tests two-run participants only | v0.2, 2026-08-24 | Clarification of the registered rule (decision 4) |
| D3 | Not stated | sub-066, sub-070 and sub-114 (constant confidence on every paced trial) excluded from H2 only; C3 taken from the other run for sub-106 and sub-112 | v0.3, 2026-11-24 | A constant modulator is not estimable (decision 9) |
| D4 | MAIA complete for all 43 | MAIA complete for 46 of 47; sub-065 has none and does not enter H6 | v0.2, 2026-08-24 | Follows from D1 |

## First-level models

| # | Registered (swc7p) | Now (v0.8) | When | Why |
|---|---|---|---|---|
| D5 | One model: awareness conditions plus confidence and arousal modulators | Three models: **A** awareness, **B** confidence, **C** salience; arousal exploratory | v0.2, 2026-08-24 | The single model was rank-deficient in 9 of 91 designs, and awareness and confidence are collinear (decision 2) |
| D6 | Paced-vs-free contrast not defined in detail | C1 = free vs the unweighted mean of the paced conditions present in the run; trailing events clipped to the scan end (3 runs) | v0.2, 2026-08-24 | The fixed 1/3 weighting is undefined in 18 of 91 runs (decision 3) |
| D7 | Confidence and arousal mean-centred within participant | Confidence mean-centred within run, over the modelled trials | v0.2, 2026-08-24 | Each run is a separate first-level model |
| D8 | Not stated | No implicit baseline; only between-condition contrasts are interpreted | v0.3, 2026-11-24 | Events tile about 94% of each run (decision 6) |
| D9 | CSF and white-matter regressors as in the sheets | Regenerated from the fMRIPrep derivatives of record; no other column changed | v0.6, 2026-11-25 | They did not reproduce from the stated inputs (median r 0.985 / 0.996; decision 15) |
| D10 | Runs combined by fixed effects for every contrast | H3 is fitted as one joint model of each participant's runs (shared task regressors, run-specific nuisance and intercept) | v0.7, 2026-11-25 | Per-run fitting discards 56 of 290 aware trials (decision 17) |
| D11 | S1 PetCO2 N = 22, S2 RVT N = 26 | Three physio pipeline bugs fixed: PetCO2 25, Biopac RVT 30. New robustness analysis S2b with RVT from the Siemens PMU cushion where the belt failed (N = 41). **One S2b gate was revised after the fallback runs' quality metrics had been seen** | v0.4, 2026-11-25 | Valid recordings had been withheld by the bugs (decisions 10–11; plan §6) |

## Hypotheses

| # | Registered (swc7p) | Now (v0.8) | When | Why |
|---|---|---|---|---|
| D12 | **H2 primary**, H3 secondary. A null H3 beside a positive H2 "will be attributed to power" | **H3 authoritative**, H2 supportive. A null H3 beside a positive H2 is not evidence that confidence rather than awareness drives the effect, because the comparison is confounded with per-participant precision | v0.2, 2026-08-24 | Confidence covaries with awareness (median within-person r = .316), so the confidence effect cannot arbitrate between them (decision 1) |
| D13 | H3: at least 3 trials per cell, N = 21; group one-sample t | H3: at least 1 aware and 1 unaware trial, **N = 45**; FSL FLAME 1, random-field-theory inference. The 3-trial rule is kept as S7 (N = 23 in the corrected sample) | v0.6, 2026-11-25 | The 3-trial rule kept 23 of 47 and over-represented infrequent detectors (decision 14) |
| D14 | H6 one-sided "replication": higher MAIA predicts greater (spared) ACC activity during paced breathing; N = 43; Harvard–Oxford dACC secondary | H6 **two-sided conceptual extension** on paced vs free; N = 46; CAREN dACC reported descriptively | v0.3, 2026-11-24 | The IEAT's MAIA effect was on interoception vs exteroception, a contrast this task lacks, so the sign is not fixed (decision 5) |
| D15 | H1 regions named in prose; no inference set listed | H1 inference set: CAREN posterior insula, central operculum and sensorimotor cortex, one-sided, Bonferroni ×3; the opposite tail is tested at the same correction | v0.5, 2026-11-25 | Fixes the test before any map exists (decision 12) |
| D16 | Regions from the Harvard–Oxford atlas | Regions from the CAREN network × AAL-3 atlas; Harvard–Oxford as atlas sensitivity only | v0.5, 2026-11-25 | Network-defined regions match the hypotheses (decision 12) |
| D17 | H2 falsification and H4 support criteria not spelled out; H4 N = 41 | H2 falsified if no effect survives in any of the three masks, dACC included; H4 supported if at least one mask survives; H4 N = 43 | v0.3, 2026-11-24 | Clarification; N follows D1 |
| D18 | H5: TOST on the **H2** effect; bounds dz = ±0.65 (PetCO2, N = 22) and ±0.59 (RVT, N = 26) | H5: TOST on the **H3** effect if H3 is supported, ±0.61 (PetCO2, N = 25) and ±0.58 (Biopac RVT, N = 28); H2 supportive, ±0.61 / ±0.59 (N = 25 / 27); H3 governs on disagreement; if H3 is not supported, H5 is read from H2 and stated about confidence only | v0.3–v0.8; final 2026-11-25 | Follows D12–D13 and D11; H3 is the authoritative awareness test (decision 19) |
| D19 | Not registered | **H7 added as a confirmatory hypothesis**: salience-network coupling with SMN and ECN, H7a paced > free (N = 47), H7b aware > unaware (N = 23); Holm ×2; global-signal interpretation gate | v0.5, 2026-11-25 | Added before any imaging analysis (decision 13). It is a hypothesis added after registration, and is reported with that history |

## Inference, sensitivity and exploratory analyses

| # | Registered (swc7p) | Now (v0.8) | When | Why |
|---|---|---|---|---|
| D20 | Sidedness not stated for whole-brain tests | Whole-brain tests two-sided; ROI tests one-sided in the predicted direction, except H6 | v0.2, 2026-08-24 | Clarification |
| D21 | S1 N = 22; S2 N = 26; S4 H2 only, N = 24 | S1 25; S2 30 (Biopac belt only); S2b 41; S4 H2 and H3, N = 26; new S6 (two-run participants, N = 43) and S7 (H3 under the 3-trial rule, N = 23) | v0.2–v0.6 | Follows D1, D11 and D13 |
| D22 | Exploratory list as registered | Added: graded-awareness model (+confidence aware / −confidence unaware); all 10 network pairs, FDR-corrected | v0.5–v0.6, 2026-11-25 | Labelled exploratory regardless of outcome |

## Registration form fields

| # | Registered (swc7p) | Now | Why |
|---|---|---|---|
| D23 | Study design text cut off after "7 cycles of [paced-breathing epoch." | Full description of the run structure | Text was truncated on the form |
| D24 | Data collection: single site, community adult sample, MRI safety screening | Adds University of Toronto REB protocol 45079, 2-hour session, $65 compensation, and the MRI contraindications from the consent form | Completes the field |
| D25 | Plan described as attached; none attached | Plan v0.8 attached; v0.1 identified by commit `aa0752b` | Correction |
| D26 | Foreknowledge option: "limited observation of the data could not influence analysis decisions" | "Authors have observed the data, but have not performed the proposed analyses" (chosen in the update) | The hypotheses are informed by analysed behaviour, and D11's gate revision used a BOLD-derived quantity |

## What the approved update contains

Checked against the OSF API on 2026-11-26. **14 fields were updated, and each
matches `osf_template_answers.md` exactly:** hypotheses, foreknowledge and its
explanation, data collection, sample size and rationale, manipulated and
measured variables, statistical models, transformations, inference criteria,
inclusion and exclusion, missing data, and other planned analyses. The
indices field was already identical.

**Three fields still carry the original 2026-08-20 text.** Where they
conflict with the updated fields, the updated fields and plan v0.8 govern:
- *Study design:* still truncated after "7 cycles of [paced-breathing epoch." (D23 not applied). The full run structure is in plan §5.
- *Additional blinding:* minor wording only ("post-registration deviations will be logged …"); no conflict.
- *Context and additional information:* still describes H1 and H6 as "replications with direction and coordinates fixed by that paper". **This conflicts with the updated hypotheses field, where H6 is a two-sided conceptual extension.** The updated H6 governs.

No file is attached to the registration (D25). The plan of record is v0.8 in
the linked repository.

## After the update

Any change after the swc7p update is approved is logged below with its date,
what changed, why, and whether the affected result had been seen.

| Date | Change | Why | Result seen? |
|---|---|---|---|
| 2026-11-26 | **Added an EXPLORATORY two-sided analysis of aware vs unaware (E1).** It uses the registered H3 FLAME fit (base model, N = 45), with no refit and no change to the model, masks or sample. Both tails are tested in the three H3 masks (left and right anterior insula, dACC), with GRF voxel FWE per tail as in the registered test and Bonferroni ×6 (3 masks × 2 tails). The REML meta-analysis is reported for both directions. The registered whole-brain test is already two-sided and is reported as it stands. Script: `Scripts/intero2023_explore_h3.py`. | Requested by the PI after the registered one-sided H3 test was null. Aware < unaware was not registered. | **Yes.** The registered H3 result (aware > unaware, no mask surviving) had been seen. This analysis is exploratory and cannot support or overturn H3. |
| 2026-11-26 | **Added EXPLORATORY whole-brain maps of the physiological regressors themselves (E2).** In the registered +RVT and +PetCO2 Model A fits, the RVT ⊗ RRF beta and the sum of the four PetCO2-lag betas (the net CO2 response) are taken as first-level contrasts. They are combined over runs by the same precision-weighted fixed effects, then tested at group level. The test is a one-sample sign-flip permutation (10,000 permutations, seed 20261126), two-sided, with cluster-forming p < .001 per tail and cluster FWE p < .05, as for the registered whole-brain tests. Samples are S2b (RVT, N = 41), S2 (Biopac RVT, N = 30) and S1 (PetCO2, N = 25). The four single-lag maps are kept, not tested. Script: `Scripts/intero2023_explore_physio.py`. | Requested by the PI: what are the whole-brain correlates of physiological variation? The registration uses these regressors as nuisance terms only. | **Yes.** All registered results had been seen, including the S1/S2/S2b comparisons. These maps were not. Exploratory. |
| 2026-11-26 | **EXPLORATORY E3: change vs No Change, then awareness within change-responsive regions.** (a) In a joint two-run Model A, fitted as for H3, the contrast is the unweighted mean of the change conditions present (aware, unaware, wrong-direction) minus No Change. It is computed in every participant with at least one No Change trial. The group test is a sign-flip permutation, 10,000 permutations, two-sided, with cluster-forming p < .001 per tail and cluster FWE. (b) The voxels in surviving clusters (either sign) form a mask. Within it, the registered H3 FLAME z map is tested in both tails by GRF voxel FWE, using the mask's resels and the H3 smoothness, Bonferroni ×2. The change contrast weights aware and unaware equally, so its weights are orthogonal to aware − unaware. Script: `Scripts/intero2023_explore_awareness.py`. | Requested by the PI after the registered awareness tests were null. | **Yes.** All registered results had been seen. Exploratory. |
| 2026-11-26 | **EXPLORATORY E4: TFCE whole-brain test of aware vs unaware.** The H3 joint-model C2 effects (N = 45) are tested by sign-flip permutation with TFCE, two-sided, **5,000 permutations** (fewer than the registered analyses, for runtime), seed 20261126. It is equal-weighted: a precision-weighted model is not exchangeable under sign-flipping, which is why the registered H3 used GRF. Same script. | As above. TFCE is more sensitive than a p < .001 cluster threshold to weak, spatially extended effects. | **Yes.** Exploratory. |
| 2026-11-26 | **EXPLORATORY E5: Bayes factors for the null.** These are JZS one-sample Bayes factors (Rouder et al., 2009) on each participant's ROI-mean effect in the three H3 masks: H3 C2 joint (N = 45), with H2 C3 (N = 44) alongside. They are reported as BF10 (two-sided) and BF+0 (directional, as registered). The Cauchy prior scale is r = 0.707, with r = 0.5 and 1.0 as robustness. Equal-weighted. Same script. | To show whether the H3 null is evidence of absence or merely inconclusive. | **Yes.** Exploratory. |
