# Localizer analysis plan: DRAFT (v0.2, 2026-11-27)

**Status.** This draft is for the PI's review. Nothing has been fitted to localizer imaging data. It merges two plans written independently and blind to localizer BOLD: one by Claude, one by a separate review agent. Where they differed, the choice taken here is marked **[decision]**.

**What has been seen** (to be declared in any addendum):
- the localizer task logs, for timing, block order, button presses and the key-5 overwrite (see §1);
- one localizer-derived physiology number already in the registered plan: PetCO2–global-signal coupling of +0.126 (§6 of `intero2023_imaging_prereg.md`);
- all breathing-run results, registered and exploratory E1–E10.

No localizer BOLD map, ROI value or statistic has been computed.

## 1. Design and behaviour (from `InteroPilot.4.51.py`, instruction images and logs)
Two runs of about 196 volumes each (TR 2 s, about 382 s of task). Each run has 10 blocks: two shuffled sets of five tasks. Each block is a 6 s instruction, 2 s fixation and 30 s of task, so there are 4 blocks per condition per participant.

| Condition | Instruction | Task |
|---|---|---|
| Intero | "Track the Breath": 1 inhale, 2 exhale | Free breathing; presses follow the breath; the circle follows the presses |
| Extero | "Track the Circle": 1 growing, 2 shrinking | Presses follow a pacer circle; breathing not instructed |
| Match | "Match Breath to Circle" | Paced breathing with presses |
| Self | "Does the word describe you?" (yes / no) | Six trait words, 5 s each |
| Rest | "Focus on the fixation cross" | Fixation |

**Pacer rate.** Protocol v2 (P60 onwards) uses a fixed 12/min. Protocol v1 (P51–59) used each participant's calibrated rate, about 8–16/min.

**Log quality.**
- Intero presses are clean: median 10 per block, press-interval CV 0.12. 16 of 176 blocks have fewer than 4 presses.
- In Extero, Match and Self the logged response is the last key pressed, and the scanner trigger (key 5) usually overwrites it. Accuracy and yes/no answers therefore cannot be scored.
- 103 of 105 scanned runs have complete logs. sub-058 run-2 and sub-114 run-2 have 5 blocks; check them against scan length.
- Respiration belt data exist for 97 runs, 72 of them with exact trigger timing.

## 2. First level
- **Pipeline.** As in the registered breathing runs: joint two-run model, PSC per run, AR(1), the same nuisance set and FD spike rule.
- **Regressors.** Five 30 s boxcars and one 8 s instruction-plus-fixation regressor. Intero blocks with fewer than 4 presses go into a separate no-interest regressor. No press regressors: at 10–12 per minute they convolve to the block regressor.
- **Smoothing.** 6 mm for group maps.
  - **[decision]** fROIs and pattern analyses use 4 mm, a compromise between the unsmoothed version one plan proposed and the 6 mm of the group maps.
  - Breathing-run contrasts are re-extracted at the same 4 mm; an unsmoothed version is the sensitivity check.
- **High-pass filter.** Each condition occurs only twice per run, so its contrast energy is low-frequency.
  - Before any fit, compute each contrast's efficiency under the 128 s cosine set, from the design matrices alone.
  - If the loss exceeds 20%, use a 200 s cutoff.
- **Physio.** The primary model has none, so it runs at full N.
  - Sensitivity: + RVT⊗RRF, and + PetCO2 where available, in runs with exact triggers. Report as subset-dependence, as for S1/S2.
  - E2 showed RVT loads on about 30% of the brain, so this sensitivity analysis is required.
- **Contrasts.**

| # | Contrast | Question |
|---|---|---|
| C1 | Intero − Extero | Target of attention: breath vs visual; motor and visual demands roughly matched |
| C2 | Intero − Self | Breath vs conceptual self |
| C3 | min-t(C1, C2) | Selective for interoceptive attention |
| C4 | Match − Extero | Same pacer and presses, plus coupling of the breath to it |
| C5 | Match − Intero | Paced vs free breathing, presses matched: **independent-run H1 analogue** |
| C6 | Self − Rest | Self-reference (with reading) |
| — | each condition − Rest | Descriptive |

Per-run versions of all contrasts are saved for split-half analyses.

- **Exclusions.**
  - The registered run-level motion gate.
  - Intero contrasts need at least 3 of 4 compliant Intero blocks.
  - fROI definition needs both runs; group maps accept one run.

## 3. Behavioural and physiological measures (no BOLD; can run now)
- **Breath-tracking accuracy (Intero).** Belt-detected inhale and exhale onsets are matched to presses. Per participant, measure:
  - the proportion of breaths with a correct 1-then-2 pair;
  - the press–belt lag (mean and SD).

  This is the study's only objective interoceptive accuracy measure. Relate it to MAIA, to breathing-run detection rate and confidence, and to C1.
- **Entrainment and adherence.** Belt–pacer phase-locking in Extero (unintended) and in Match (intended).
- **Rate and depth by condition.** A manipulation check, and a per-participant covariate for C1, C4 and C5.

## 4. Activation (group)
- **Inference.**
  - FLAME 1 with GRF voxel FWE in masks, and whole-brain cluster FWE at z > 3.09, as for H3.
  - Permutation with voxel max-t as a check.
  - Both tails are reported.
  - **The IEAT and H1 suggest Intero < Extero (deactivation) may dominate.** L1 is therefore two-sided.
- **Masks (frozen).** CAREN AInsL, AInsR, dACC, PIns, COper, SM and SMA; the IEAT ACC sphere; CAREN networks for connectivity.
- **Covariates (secondary).**
  - Breath-tracking accuracy; entrainment.
  - Protocol v1/v2. Sensitivity: drop P51–59.

## 5. Connectivity
- **Network FC, as in H7.** Unsmoothed, the same denoising, median CAREN network signals, block windows shifted 4 s and concatenated per condition. That gives about 60 volumes per condition (per-participant SE of r about .13).
  - All 10 pairs, FDR.
  - Primary pairs for C1: SAL–SMN and SAL–DMN. For C2: DMN–SAL.
  - Report amplitude (temporal SD) alongside, and run the H7 global-signal check.
- **Background connectivity.** The same analysis on task-model residuals, to remove co-activation.
- **gPPI.** Seeds are the §6 fROIs and the a-priori masks; targets are whole brain; contrasts C1 and C4.
- **Not done.** Beta-series (4 blocks per condition) and dynamic FC.

## 6. Participant-specific ROIs (fROIs)
- **Constraint parcels.** The frozen a-priori masks: AInsL, AInsR, dACC, PIns, the IEAT sphere. The anterior insula and dACC masks are broad (y down to −20.5 and −44.5), which is exactly the case where sub-selection can help.
- **Primary rule (PI decision, 2026-11-27).** Within each parcel, each participant's top 10% of voxels by C1 t (Intero > Extero; two-run model, 4 mm), with a floor of 50 voxels.
- **Exploratory rules (PI decision, 2026-11-27),** reported alongside the primary and never substituted for it.
  - Top 10% by C2 alone (Intero > Self).
  - Top 10% by C3, the conjunction min-t(C1, C2).

  Each goes through the same reliability gate and the same transfer tests (T1–T3). Their results are labelled exploratory.
- **Other secondary rules (pre-specified).**
  - Top 10% by C5, negative direction (for H1/H6 transfer).
  - 5% and 20% sensitivity for the primary rule.
- **GSS-style group fallback.** Top 10% of the leave-one-subject-out group C1 map within the parcel.
- **Specificity control.** The parcel complement (the remaining 90%), and the bottom 10% (Extero > Intero voxels).
- **Reliability gate.** Computed on localizer data only, before any transfer. A parcel passes only if all three hold:
  1. Defined in one run and estimated in the other (both directions, averaged), C1 is above 0 at group level, p < .05, with at least 70% sign agreement.
  2. Median run-1/run-2 Dice of at least .20 (chance is .10).
  3. Within-participant cross-run spatial r exceeds between-participant r (paired test).

  If (3) fails, the group-fallback fROI becomes primary for that parcel.
- **Independence.** Selection uses localizer data only. Breathing data, awareness, confidence and MAIA never enter it.

## 7. Transfer to the breathing runs
Contrasts are re-extracted from the registered first levels (refit at 4 mm). The same participants (localizer ∩ the target set), the same estimator and the same N are used for fROI, parcel and complement.

| | Test | Inference |
|---|---|---|
| **T1** | aware − unaware (C2 joint) in the AIns and dACC fROIs, one-sided | ×3; JZS BF (r = .707); FLAME-weighted, as H3 |
| T2 | confidence modulator (C3); change − No Change (E3) | Two-sided |
| T3 | paced − free (C1) in C5 fROIs; MAIA × paced − free in the ACC-sphere fROI (H6 transfer) | As H1/H6 |
| T4 | **Per-epoch pattern expression.** LSS betas for about 14 paced epochs; index = dot product with the participant's own C1 map (leave-one-subject-out group map as an alternative) within SAL ∪ SMN; mixed model: index ~ confidence + outcome + (1 \| participant) | Exploratory; confidence is the better-powered target |
| T5 | fROI-node coupling, aware vs unaware and paced vs free | Exploratory |

**"More sensitive than the within-run analyses" is defined in advance.**
- Required: T1 is supported in the fROI, **and** a paired test of fROI minus parcel complement on the same contrast is above 0 (bootstrap CI on the dz difference).
- Not sufficient: "significant in the fROI, not in the mask".
- BFs are reported next to E5 (BF01 5.1 / 2.2 / 5.8). A BF01 above 3 in the fROIs is itself a result: it rules out spatial dilution as the reason the H3 null occurred.
- Power: with N = 45, the minimum detectable dz is about .43.

## 8. Confirmatory vs exploratory
**Proposed pre-registered addendum,** to swc7p or a linked registration, frozen before any localizer BOLD is read:
- **L1:** C1 ≠ 0 in AInsL, AInsR and dACC (two-sided, ×3).
- **L2:** C5 < 0 in the H1 somatomotor masks (replicates H1 in independent runs).
- **L3:** MAIA × C1 in the IEAT ACC sphere (two-sided). This is the direct IEAT replication; H6 was its conceptual extension.
- **L4:** the fROI reliability gate (§6) is met.
- **T1:** aware > unaware in the C1 fROIs, and greater than in the parcel complement.
- **T3b:** MAIA × paced − free is positive in the ACC-sphere fROI.

The transfer tests are "pre-specified, target data known". The breathing-run awareness results are already known (null). What is new is the localizer-defined measurement.

**Exploratory:** everything else, including the C2 and C3 fROI rules, connectivity, gPPI, RSA, decoding, gradients, T2, T4, T5 and phase-locked responses.

**Decided 2026-11-27:** the PI chose to add this plan to the pre-registration (swc7p).

## 9. Other techniques (ranked; exploratory)
1. Per-epoch pattern expression (T4).
2. Physiology-informed between-subject models (tracking accuracy, entrainment, MAIA).
3. Leave-one-subject-out decoding of Intero vs Extero from per-run contrast maps. Within-subject classification is not viable: 4 blocks per class, 2 folds.
4. Cross-run crossnobis RSA of the 5 conditions in fROIs. Model RDMs for attention target, breathing mode, motor and visual motion are partly collinear, so this is descriptive.
5. Insula posterior → anterior gradient of C1 (a per-participant slope), and projection onto cortical gradients.
6. Inhale- vs exhale-locked responses by condition (physio subset).

## 10. Risks
- **C1 is not a pure attention contrast.** It also differs in:
  - agency (self-driven vs external circle);
  - visual tracking demand;
  - press predictability;
  - breathing rate and depth.

  Interpret by convergence across C1, C3 and C4.
- **C2 and C3 are confounded by press rate** (about 20 presses vs 6), visual motion and reading. Self has no semantic control and no usable responses.
- **Individual maps are noisy** (4 blocks per condition), hence the reliability gate and the fallback.
- **fROIs cannot fix the target contrast.** It has a median of 2 misses and ICC .2–.45, so expect a modest gain at most.
- **Breathing leaks into BOLD.** Extero entrainment varies by person, and v1 participants were paced near their own rate.
- **Block order is shuffled and unbalanced,** so carry-over is possible (E8/E10 showed carry-over in the breathing runs).
- **CAREN masks are on the FSL MNI grid** (already declared).
- **Multiplicity across fROI rules.** One primary rule only.
