# OSF swc7p: update 2, the localizer addendum (paste-ready text)

The plan of record is `Preregistration/intero2023_localizer_prereg.md` (v1.0) in the linked GitHub repository. Each answer below points into it by section (§).

## Justification for the update
This update adds an analysis plan for the two localizer runs (`task-localizer`). They were collected in every session but were not covered by the registered plan, except as a breathing-at-rest reference for PetCO2.

The plan was written before any localizer imaging data were analysed. No localizer BOLD map, ROI value or statistic has been computed.

What was examined before registration:
- the localizer task logs;
- framewise displacement;
- NIfTI headers;
- design efficiency;
- the raw respiration-belt recordings, for quality, per-condition breathing, and the distribution and reliability of breath-tracking accuracy.

No localizer behavioural or physiological measure has been related to any outcome. Three plan choices were made after these checks:
- the high-pass cutoff;
- the breath-tracking measure;
- the exclusion of two partial runs.

These are disclosed in §2–§4. The breathing-run results registered in this project (OSF swc7p, plan v0.8) and the exploratory analyses that followed them are known to the authors. Tests that use breathing-run data are therefore labelled "pre-specified, target data known".

## Hypotheses
- **L1.** Attending to the breath vs attending to a visual circle (Intero − Extero; matched buttons and circle motion) changes activity in the salience-network masks: left and right anterior insula, dACC. Two-sided.
- **L2.** Paced breathing lowers activity relative to free breathing in the somatomotor masks (posterior insula, central operculum, sensorimotor), in runs independent of the registered H1 test (Match − Intero < 0).
- **L3.** Interoceptive awareness (MAIA) moderates the breath-vs-visual attention effect in the IEAT ACC sphere. Two-sided. This is a direct replication of the IEAT finding.
- **L4.** Participant-specific interoceptive-attention ROIs are reliable across the two localizer runs (§7 gate).
- **T1.** In participant-specific ROIs defined by Intero > Extero within the anterior insula and dACC masks, aware > unaware change detection in the breathing runs, and more so than in the remainder of the same mask.
- **T3b.** MAIA × (paced − free) is positive in the participant-specific ACC-sphere ROI.

## Design
Two localizer runs per participant. Each has 10 blocks of 30 s: two shuffled sets of five tasks.
- Track the breath with button presses (Intero).
- Track a pacer circle with button presses (Extero).
- Breathe with the circle while pressing (Match).
- Trait-word self-judgements (Self).
- Fixation (Rest).

Each block is preceded by a 6 s instruction and a 2 s fixation. Details, including protocol v1 timing: §1.

## Sampling and existing data
- **Existing data.** All data were collected before this update; 54 participants have localizer runs.
- **After exclusions** (§2: two partial runs; registered motion gate): 102 runs, 54 participants.
- **Ns:**
  - Intero contrasts: N = 48.
  - Participant-specific ROIs: N = 44.
  - Transfer tests: N = 40, minimum detectable dz .45.

## Variables
- Contrasts C1–C6 (§3).
- MAIA (as registered).
- Breathing-run contrasts as registered (C2 joint aware − unaware, C1 paced − free), refitted at 4 mm.
- Exploratory physiological measures (§4):
  - breath-tracking phase-locking;
  - pacer entrainment and adherence in both tasks;
  - breathing rate, depth and RVT.

## Analysis plan
- **First level** (§3).
  - Joint two-run GLM as registered, with Rest as the implicit baseline.
  - **200 s high-pass.** The pre-set rule replaces 128 s when the 128 s efficiency loss exceeds 20%; it was 49–55%. 128 s and 320 s are sensitivity analyses.
  - 6 mm smoothing for group maps; 4 mm for ROIs.
- **Group** (§5). FLAME 1; GRF voxel FWE in frozen masks; whole-brain cluster FWE at z > 3.09.
- **Participant-specific ROIs** (§7).
  - Each participant's top 10% of voxels by Intero > Extero within each mask (minimum 50), from the two-run model.
  - Reliability gate; leave-one-subject-out group fallback.
  - Intero > Self, and the conjunction of both, are exploratory.
- **Transfer** (§8). The same participants, estimator and N for ROI, mask and complement.

## Inference criteria
- L1: two-sided, Bonferroni ×3.
- L2: one-sided, ×3.
- L3: two-sided.
- T1: one-sided, ×3, plus a paired ROI-minus-complement test (bootstrap CI on the dz difference). JZS Bayes factors (r = .707) are reported throughout.
- **L1 breathing robustness.** In participants with RVT on both runs, the +RVT effect must keep its sign and ≥ 50% of the base effect in the same participants. Otherwise L1 is reported as breathing-dependent.
- **"More sensitive than the within-run analyses."** T1 must be supported **and** the ROI-minus-complement difference must be > 0. "Significant in the ROI, not in the mask" is not sufficient.

## Data exclusion
- Two runs whose task stopped mid-run (sub-058 run-2, sub-114 run-2).
- The registered run-level motion gate.
- Intero compliance: a block needs ≥ 4 presses; a participant needs ≥ 3 of 4 compliant blocks (2 of 2 with one run).
- ROI eligibility needs two usable runs (§2).

## Missing data
- Participants without MAIA are excluded from L3 and T3b only.
- Physiological covariates and the +RVT co-report use available runs, reported as subset results.

## Other (exploratory)
- Network connectivity and gPPI (§6).
- Per-epoch pattern expression in the breathing runs (§8 T4).
- Decoding, RSA and gradients (§10).
- Breath-tracking phase-locking and pacer adherence across both tasks, and their relations to MAIA, breathing-run behaviour and brain measures (§4).
