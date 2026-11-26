# Deviations from the OSF registration

**Registration of record:** OSF Preregistration *Neural correlates of awareness
of respiratory change during paced breathing*, https://osf.io/swc7p, registered
**2026-08-20 12:23 UTC** (08:23 EDT), public, no embargo, no later update. The
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
from the registration. They are filed together as an OSF registration update.

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
| D26 | Foreknowledge option: "limited observation of the data could not influence analysis decisions" | Author to choose; see the note in the update text | The hypotheses are informed by analysed behaviour, and D11's gate revision used a BOLD-derived quantity |

## After the update

Any change after the swc7p update is approved is logged below with its date,
what changed, why, and whether the affected result had been seen.

| Date | Change | Why | Result seen? |
|---|---|---|---|
| | | | |
