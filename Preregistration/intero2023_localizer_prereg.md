# Localizer analysis plan: linked secondary-data preregistration (v1.0, 2026-11-27)

**Status.** This is the plan of record for the two localizer runs (`task-localizer`). **Registered at OSF: https://osf.io/p5ztc (Secondary Data Preregistration, approved 2026-11-27), linked to swc7p (same sample).** Any later change is logged in `Preregistration/localizer_deviations.md`. It is not an update of swc7p. The paste-ready form answers are in `osf_localizer_secondary_prereg.md`. Nothing has been fitted to localizer imaging data. History:
- **v0.1:** merged from two plans, written independently and blind to localizer BOLD (one by Claude, one by a separate review agent).
- **v0.2:** the PI's choice of fROI rule.
- **v1.0:** adds the step-1 QC results and the PI's decisions of 2026-11-27: 200 s high-pass; phase-locking as the exploratory tracking measure; the breathing-robustness rule for L1.

**Seen before registration.** Every item below is declared.
1. **Localizer task logs:**
   - timing and block order;
   - attempt selection;
   - the v1 timing structure;
   - the two partial runs;
   - Intero press counts and compliance.
2. **Localizer data other than BOLD:**
   - fMRIPrep framewise displacement only (the motion gate);
   - NIfTI headers (volume count, TR);
   - design efficiency and VIFs from the design matrices.
3. **Raw localizer belt (Biopac):**
   - belt and trigger-timing grades;
   - belt polarity;
   - breathing rate, depth, RVT and pacer phase-locking per condition, with the paired condition tests;
   - breath-tracking accuracy: distribution, chance level, missingness, the press lead, and run-1 vs run-2 reliability.
4. **One localizer physiology number already in the registered plan:** PetCO2–global-signal coupling of +0.126 (`intero2023_imaging_prereg.md` §6).
5. **All breathing-run results,** registered and exploratory E1–E10 (`Analysis/RESULTS_SUMMARY.md`).

**Not seen:** any localizer BOLD map, ROI value or statistic; any relation of a localizer behavioural or physiological measure to MAIA, breathing-run behaviour, demographics or brain data. QC details: `Analysis/localizer_qc/README.md`; script `Scripts/intero2023_localizer_qc.py` (commit 334d0c1).

## 1. Design
Two runs of about 196 volumes (TR 2 s). Each run has 10 blocks: two shuffled sets of five tasks, so there are 4 blocks per condition per participant. Each block is an instruction (6 s), a fixation (2 s) and 30 s of task.

**Protocol v1 (P51–59).** Self and Rest have no fixation, so their task starts about 6.2 s after the instruction instead of about 8.2 s. Extero and Match run whole pacer cycles (30.1–34.4 s).

| Condition | Instruction | Task |
|---|---|---|
| Intero | "Track the Breath": 1 inhale, 2 exhale | Free breathing; presses follow the breath; the circle follows the presses |
| Extero | "Track the Circle": 1 growing, 2 shrinking | Presses follow a pacer circle; breathing not instructed |
| Match | "Match Breath to Circle" | Paced breathing, with presses |
| Self | "Does the word describe you?" (yes / no) | Six trait words, 5 s each |
| Rest | "Focus on the fixation cross" | Fixation |

**Pacer.** v2 (P60 onwards): 12/min. v1: each participant's calibrated rate (7–16/min).

**Responses.** Intero presses are clean: median 11 per block. In Extero, Match and Self the scanner trigger (key 5) overwrites the logged response, so accuracy and yes/no answers cannot be scored. This is a limitation.

**Events.** `Data/events/<sub>_task-localizer_run-<k>_events.tsv` (105 runs) uses exact block durations. Six v1 runs end before the task and are clipped to the scan end, as for the breathing runs.

## 2. Data and exclusions
- **Partial runs:** sub-058 run-2 and sub-114 run-2 are excluded. The task stopped in the 5th block while the scanner ran on, and one condition is missing from each.
- **Motion:** the registered run-level gate (mean FD > 0.5 mm, or > 25% of volumes with FD > 0.5 mm) excludes 2 runs.
- **Result:** 102 runs, 54 participants (48 with two runs).
- **Intero compliance:** a block is compliant with ≥ 4 presses. A participant enters Intero contrasts with ≥ 3 of 4 compliant Intero blocks, or 2 of 2 if they have only one usable run.
  - Non-compliant Intero blocks are modelled by a separate regressor of no interest.
  - N for Intero contrasts: **48**.
  - Contrasts without Intero (C4, C6): every usable participant (**54**).
- **fROI eligibility:** two usable runs, ≥ 3 of 4 compliant blocks, and ≥ 1 compliant Intero block in each run: **N = 44**.
- **Transfer sets:** fROI ∩ the registered H3 set gives **N = 40**. At N = 40 the minimum detectable dz is **.45** (two-sided .05, 80% power).

## 3. First level
- **Pipeline.** As the registered breathing runs: joint two-run model, PSC per run, AR(1), the same nuisance set and FD spike rule.
  - **Exception:** fMRIPrep's `cosine_*` columns are 128 s, so they are dropped.
- **Regressors.**
  - Four task boxcars (Intero, Extero, Match, Self) and one instruction-plus-fixation regressor.
  - **Rest is the implicit baseline.** With Rest modelled explicitly the regressors tile the run (VIF 9–12); with it implicit, VIF is about 2. Differential contrasts change by ≤ 1%.
  - Non-compliant Intero blocks go into a separate regressor.
  - No press regressors.
  - HRF: SPM canonical.
- **High-pass: 200 s cosine set (primary).** At 128 s the median efficiency loss is 49–55%, which fires the pre-set trigger (> 20%). At 200 s the loss is still 22–27% median, and this is stated as a limitation.
  - **Sensitivity:** 128 s and 320 s (8–11% loss).
- **Smoothing.**
  - 6 mm for group maps.
  - 4 mm for fROI definition and pattern analyses.
  - Breathing-run contrasts used in transfer are refitted at 4 mm (unsmoothed as sensitivity), with no other change to the registered model.

**Contrasts** (per-run versions are also saved for split-half analyses):

| # | Contrast | Question |
|---|---|---|
| C1 | Intero − Extero | Target of attention: breath vs visual circle, with buttons and circle motion roughly matched |
| C2 | Intero − Self | Breath vs conceptual self |
| C3 | min-t(C1, C2) | Selective for interoceptive attention |
| C4 | Match − Extero | Same pacer and presses, plus breath coupling |
| C5 | Match − Intero | Paced vs free breathing (presses matched): **independent-run H1 analogue** |
| C6 | Self (vs Rest) | Self-reference, with reading |
| — | Intero, Extero, Match (each vs Rest) | Descriptive |

- **Physiology.** The primary model has none, so it runs at full N.
  - The +RVT⊗RRF model is a **required co-report** for L1 (§5). It runs in participants with RVT on both runs (about 29 of the 44 fROI-eligible).
  - +PetCO2 is a sensitivity analysis where available.

## 4. Behavioural and physiological measures (exploratory)
Purpose: to know whether participants were breathing along with the pacer, and how well they tracked their own breath, in **both** the localizer and the change-detection runs.

Belt polarity is set per participant from the Match blocks (sign change between sub-068 and sub-069). sub-071 is unresolved and is excluded from phase measures.

- **Breath-tracking accuracy (Intero).**
  - **Primary measure:** phase-locking value (PLV) of the presses to the belt's breath phase, in compliant blocks. It is reported against a circular-shift chance level.
  - **Secondary measure:** a press-pair score centred on each run's press lead, estimated from the *other* Intero block (cross-validated).
  - **Disclosed:** the pre-specified fixed-window score was at chance because presses lead the belt extremum by about 0.85 s, and the centred score was devised after that.
  - Measured in 38 participants (29 with both runs). Run-to-run reliability: ICC .43 (PLV), .56 (pair score).
- **Pacer coupling, localizer.**
  - Extero entrainment (unintended): PLV minus a same-participant virtual-pacer null.
  - Match adherence (intended): PLV, the same way.
- **Pacer coupling, change-detection runs.**
  - The same PLV, computed on the paced epochs, alongside the registered adherence measure (`adhere_r`, `Data/resp_runs.csv`). This puts adherence in both tasks on one scale.
- **Rate, depth and RVT** per condition. These serve as the manipulation check, and as covariates (§5).
- **Exploratory relations, all reported as exploratory:**
  - tracking PLV with MAIA;
  - tracking PLV with breathing-run detection rate and mean confidence;
  - tracking PLV with breathing-run adherence;
  - tracking PLV with C1 in the fROIs;
  - localizer Match adherence with breathing-run adherence;
  - breathing-run adherence as a moderator of H1 (paced − free).

## 5. Activation (group)
- **Inference.**
  - FLAME 1. GRF voxel FWE within masks; whole-brain cluster FWE at z > 3.09, as for H3.
  - Permutation voxel max-t is reported as a check.
  - Both tails are reported.
- **Masks (frozen).** CAREN AInsL, AInsR, dACC, PIns, COper, SM and SMA; the IEAT ACC sphere; CAREN networks (§6).
- **Breathing-robustness rule for L1.** The primary L1 result comes from the model without physiology. In the participants with RVT on both runs, the +RVT model is fitted and compared with the base model *in the same participants*.
  - L1 is **breathing-robust** if the +RVT mask-mean effect has the same sign and is ≥ 50% of the base effect.
  - Otherwise L1 is reported as **breathing-dependent**.
- **Group covariates (secondary).**
  - ΔRVT (Intero − Extero, participant mean).
  - Extero entrainment (PLV minus null).
  - Protocol v1/v2. Sensitivity: drop P51–59.

## 6. Connectivity (exploratory)
- **Network FC, as in H7.**
  - Unsmoothed data; median CAREN network signals; block windows shifted 4 s; concatenated per condition (about 60 volumes per condition).
  - All 10 pairs, FDR.
  - Primary pairs for C1: SAL–SMN and SAL–DMN. For C2: DMN–SAL.
  - Amplitude (temporal SD) and the H7 global-signal check are reported alongside.
- **Background connectivity:** the same, on task-model residuals.
- **gPPI:** seeds are the §7 fROIs and the a-priori masks; targets are whole brain; contrasts C1 and C4.
- **Not done:** beta-series (4 blocks per condition) and dynamic FC.

## 7. Participant-specific ROIs (fROIs)
- **Constraint parcels:** the frozen AInsL, AInsR, dACC, PIns and the IEAT sphere.
- **Primary rule (PI decision).** Within each parcel, each participant's top 10% of voxels by C1 t (Intero > Extero; two-run model, 4 mm), with a minimum of 50 voxels.
- **Exploratory rules (PI decision),** reported alongside the primary and never substituted for it:
  - top 10% by C2 alone (Intero > Self);
  - top 10% by C3, the conjunction min-t(C1, C2).

  Each goes through the same gate and the same transfer tests.
- **Other secondary rules:**
  - top 10% by C5, negative direction (for H1/H6 transfer);
  - 5% and 20% sensitivity for the primary rule.
- **GSS-style group fallback:** top 10% of the leave-one-subject-out group C1 map within the parcel.
- **Specificity control:** the parcel complement (the remaining 90%), and the bottom 10% (Extero > Intero).
- **Reliability gate (L4),** on localizer data only, before any transfer. A parcel passes if all three hold:
  1. Defined in one run and estimated in the other (both directions, averaged), C1 is above 0 at group level, p < .05, with ≥ 70% sign agreement.
  2. Median run-1/run-2 Dice ≥ .20 (chance is .10).
  3. Within-participant cross-run spatial r exceeds between-participant r (paired test).

  If (3) fails, the group fallback becomes primary for that parcel.
- **Independence.** Selection uses localizer data only. Breathing data, awareness, confidence and MAIA never enter it.

## 8. Transfer to the breathing runs
Same participants (fROI ∩ target set), same estimator and same N for fROI, parcel and complement.

| | Test | Inference |
|---|---|---|
| **T1** | aware − unaware (C2 joint) in the AIns and dACC fROIs, one-sided | ×3; FLAME-weighted, as H3; JZS BF (r = .707) |
| T2 | confidence modulator (C3); change − No Change (E3) | Two-sided; exploratory |
| T3 | paced − free (C1) in C5 fROIs; **T3b**: MAIA × paced − free in the ACC-sphere fROI | As H1/H6 |
| T4 | Per-epoch pattern expression: LSS betas for about 14 paced epochs; index = dot product with the participant's own C1 map (leave-one-subject-out group map as an alternative) within SAL ∪ SMN; index ~ confidence + outcome + (1 \| participant) | Exploratory |
| T5 | fROI-node coupling, aware vs unaware and paced vs free | Exploratory |

**"More sensitive than the within-run analyses."** Both must hold:
- T1 is supported in the fROI;
- a paired test of fROI minus parcel complement, same contrast and same participants, is above 0 (bootstrap CI on the dz difference).

"Significant in the fROI, not in the mask" is not sufficient. BFs are reported next to E5 (BF01 5.1 / 2.2 / 5.8). A BF01 > 3 in the fROIs is itself a result: it rules out spatial dilution as the reason for the H3 null.

## 9. Confirmatory hypotheses
| | Hypothesis | Test |
|---|---|---|
| **L1** | Interoceptive vs exteroceptive attention differs in the salience masks: C1 ≠ 0 in AInsL, AInsR and dACC | Two-sided, ×3; breathing-robustness rule (§5) |
| **L2** | Paced breathing lowers activity relative to free breathing in independent runs: C5 < 0 in the H1 masks (PIns, COper, SM) | One-sided, ×3 (replicates H1) |
| **L3** | MAIA × C1 in the IEAT ACC sphere (direct IEAT replication) | Two-sided |
| **L4** | The fROI reliability gate is met | §7 |
| **T1** | Aware > unaware in the C1 fROIs, and greater than in the parcel complement | §8 |
| **T3b** | MAIA × paced − free is positive in the ACC-sphere fROI | One-sided |

T1 and T3b are "pre-specified, target data known": the breathing-run results are known (H3 null; H6 supported in the a-priori sphere). What is new is the localizer-defined measurement. **Everything else is exploratory**, including:
- the C2 and C3 fROI rules;
- connectivity and gPPI;
- RSA, decoding and gradients;
- T2, T4 and T5;
- the §4 behavioural relations;
- phase-locked responses.

## 10. Other techniques (exploratory, ranked)
1. Per-epoch pattern expression (T4).
2. Physiology-informed between-subject models (§4 measures).
3. Leave-one-subject-out decoding of Intero vs Extero from per-run contrast maps.
4. Cross-run crossnobis RSA of the 5 conditions in fROIs (descriptive; the model RDMs are partly collinear).
5. Insula posterior → anterior gradient of C1 (a per-participant slope), and projection onto cortical gradients.
6. Inhale- vs exhale-locked responses by condition (physio subset, per-participant polarity).

## 11. Risks
- **C1 is not a pure attention contrast.** It also differs in:
  - agency (self-driven vs external circle);
  - press predictability;
  - **breathing.** Intero breathing is deeper than Extero (RVT dz .88), and Extero entrains breathing to the pacer.

  Hence the breathing-robustness rule, the covariates, and interpretation by convergence across C1, C3 and C4.
- **C2 and C3** are confounded by press rate (median 11 logged Intero presses vs 6 words per block), visual motion and reading. Self has no semantic control and no usable responses.
- **Four blocks per condition.** Even at 200 s the high-pass costs about 25% efficiency, so individual maps are noisy. That is why the reliability gate and the fallback exist.
- **fROIs cannot fix the target contrast:** a median of 2 misses and ICC .2–.45. Expect a modest gain at most.
- **v1 participants** (8) had a personalised pacer and different timing. Only 4 have a usable belt.
- **Block order** is shuffled and unbalanced, so carry-over is possible.
- **CAREN masks** are on the FSL MNI grid (already declared).
- **Multiplicity:** one primary fROI rule; everything else is labelled.
