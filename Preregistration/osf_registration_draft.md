# OSF Registration Draft — Neural correlates of respiratory change awareness

**Registration type:** Preregistration of a secondary data analysis (template after Mertens & Krypotos, 2019)
**Companion document:** `intero2023_imaging_prereg.md` **v0.6** (full analysis plan; to be attached to the registration). Section references (§) below point into it.
**Status:** DRAFT. To be frozen on OSF before any first-level model is fitted. Synced to plan v0.6 on 2026-11-25.

---

## Q1. What are the hypotheses under investigation?

There are seven hypotheses, each stated in full with falsification criteria in §3 of the attached plan. H1–H6 are each tested in one of three first-level models: **A** awareness, **B** confidence, **C** salience. H7 uses a separate connectivity model. Awareness and the confidence modulator are never estimated in the same design (§7).

- **H1 (Model A, N = 47).** Paced (exogenously controlled) breathing produces relative *deactivation* of the insula/operculum pathway and somatomotor cortex compared with free breathing. This replicates the Active Matching result of Farb, Zuo & Price (2023, eNeuro; "IEAT"). The direction is committed to the replication: the volitional-control account predicts the opposite, and a reversed result will be reported as a failed replication. The design has no rest baseline, so "deactivation" means relative to free breathing. Inference is in three somatomotor-network masks (posterior insula, central operculum, sensorimotor cortex), Bonferroni ×3.
- **H3 (authoritative awareness test; Model A, N = 45).** Aware (correct-direction) > unaware (responded "No Change") change trials, in bilateral anterior insula and dACC. Inference is ROI-restricted, Bonferroni ×3. Every participant with an aware and an unaware trial in the same run enters. Unaware counts run from 1 to 8, so the group model weights participants by precision (FSL FLAME 1). The earlier ≥ 3-trial rule kept 23 participants and over-represented less frequent detectors. It is retained as sensitivity analysis S7.
- **H2 (supportive awareness test; Model B, N = 44).** Trial-by-trial confidence on paced trials correlates positively with activity in the same three ROIs. Model B's confidence effect absorbs the awareness effect, so a null H3 alongside a positive H2 is not evidence that confidence rather than awareness drives the effect. That comparison is confounded with per-participant precision: H3 rests on about 3–4 trials per cell.
- **H4 (Model C, N = 43).** Salience (step vs ramp) effects fall in the salience-network ROIs. This is consistent with the behavioural finding that salience affects noticing, not direction discrimination.
- **H5.** Physiological correction attenuates but does not abolish the awareness-related effects. It is adjudicated by a preregistered TOST equivalence criterion on the H2 effect (bounds dz = ±0.61 for PetCO2 at H2 N = 25 and ±0.59 for Biopac-belt RVT at H2 N = 27; interpretation grid fixed in §12.4). The same comparison for H3 is descriptive: 25 and 28 of its participants have the two physiological regressors. A robustness version (S2b) adds RVT from the scanner's PMU respiration cushion where the Biopac belt failed (H2 N = 38; H3 39 of 45).
- **H6 (conceptual extension; Model A, N = 46; two-sided).** MAIA interoceptive sensibility moderates paced-vs-free activity in a 10 mm sphere at the IEAT's ACC coordinate (MNI −4, +12, +36). The IEAT's MAIA effect was on interoception vs exteroception, a contrast this task does not contain. The sign on this study's contrast is therefore not fixed by prior work.

- **H7 (connectivity; both parts confirmatory).** The salience network couples more strongly with the somatomotor and executive networks (CAREN SAL–SMN, SAL–ECN), as follows. **H7a** (N = 47): during paced than during free epochs. **H7b** (N = 23, the participants with ≥ 3 epochs per cell): on aware than on unaware change epochs. The test is a two-sided paired t on Fisher-z coupling, Holm-corrected across the two pairs. A supported effect counts as network coupling only if it survives partialling the global signal, because pacing changes CO2 by design. H7b is powered only for large effects (dz ≥ 0.61).

## Q2. What data will be used, and how were they generated?

MRI_Intero2023 is a single-session fMRI study (3T Siemens, TR = 2 s, fMRIPrep 24.1.1 derivatives). Participants performed a paced-breathing change-detection task with free-breathing periods and rating windows. The paced trials varied change direction (faster / slower / no change) and salience (step / ramp). Respiration belt and expired CO2 (capnograph) were recorded concurrently, along with scanner PMU physiology. Data collection is complete (2023-09 – 2024-06). Nothing further will be collected.

## Q3. Have the data been accessed or analysed? (full disclosure)

**Yes — everything except the functional imaging.** This is the core disclosure and is stated identically in §0 of the attached plan:

*Fully analysed before this registration:* task logs and event timing; all trial-level behaviour (three-category identification, ordinal SDT, confidence/metacognition); respiration-belt processing, including pacer adherence; capnography and PetCO2 extraction; RVT; response-function validation; and all QC and exclusion decisions.

*Not examined in any form:* the functional images beyond automated preprocessing (fMRIPrep, and MRIQC image-quality metrics used for run QC). No first-level model has been fitted, no contrast image exists, and no statistical map has been computed or viewed.

*Computed before freezing, from timing and confound sheets only:* the design matrices themselves. These are rank, collinearity, contrast efficiency, event coverage, trial counts and estimability. No BOLD voxel was read.

*Physiological re-examination before freezing (prereg §6, v0.4):* three pipeline bugs that had withheld valid recordings were fixed. A PMU-belt fallback for RVT was developed, using the global signal as its validity check. One gate was revised after the fallback runs' quality metrics had been seen; that sequence is disclosed in full in §6, and the fallback is used only in robustness analysis S2b.

*Before freezing, v0.6:* the CSF and white-matter regressors were regenerated from the fMRIPrep derivatives of record. They had come from an earlier anatomical pass; no other column changed. The H3 inclusion rule was revised on trial counts alone.

*ROI masks (prereg §8, v0.5):* built before freezing from the CAREN atlas and fMRIPrep's run-level brain masks, which fix the grid and the coverage mask. No BOLD intensity was read. The masks are frozen by SHA-256 checksum. The generating code and checksums are published; the mask files themselves are not, because they derive from a third-party atlas.

*Documented exceptions:* fMRIPrep's `global_signal` confound column was used to verify that the physiological regressors carry the physiologically correct sign (§6). Its `framewise_displacement` column was used to apply the preregistered motion gate. No task regressor has been related to any imaging measure.

*Consequence, stated plainly:* H2–H4 are informed by known behavioural results. They are untested predictions about brain data, but they are not blind hypotheses. H1 is a directional replication of an independent published sample (IEAT). H6 takes its region from that sample but registers no direction.

## Q4. How was the sample determined, and what are the exclusion rules?

The sample is the full collected cohort minus preregistered exclusions, all already applied and enumerated in §4:
- **Identifiers:** pilot/QA IDs are excluded.
- **sub-051 and sub-056:** sub-051 has no MRI. sub-056's run-1 is truncated and its run-2 fails the motion gate.
- **Motion gate:** a run fails if mean FD > 0.5 mm or > 25% of volumes exceed 0.5 mm. It was applied 2026-08-14 from the FD column only. It fails 12 of 103 runs and removes 5 participants entirely.
- **One-run rule:** a participant enters when at least one breathing run passes. Run count is not an inclusion criterion.
- **Estimability rule:** a contrast counts for a participant when each of its cells has ≥ 3 trials, pooled across runs. H3 is the exception: it needs one aware and one unaware trial in the same run, and its group model weights participants by precision.

**Post-gate sample: 90 runs, 47 participants.** Final Ns: H1 47, H2 44, H3 45, H4 43, H6 46, H7a 47, H7b 23. Sensitivity samples: S1 (PetCO2) 25, S2 (Biopac RVT) 30, S2b (RVT incl. PMU fallback) 41, S4 (adherence) 26, S6 (two-run only) 43, S7 (H3 under the ≥ 3-trial rule) 23. A registered check shows the motion gate does not select against pacer-adherent breathers (fd–adherence r = −0.003). sub-082 and sub-104, excluded in an earlier draft for "shifted run labels", are included. That problem affected only an obsolete copy of the data (§4).

## Q5. Which variables enter the analysis?

First level uses three models (§7), all sharing one nuisance set:
- **Model A** has paced_aware, paced_unaware, paced_wrongdir, paced_nochange, free, rating and subrating regressors.
- **Model B** collapses the paced epochs into one regressor with a mean-centred confidence modulator.
- **Model C** splits the paced change epochs by salience.

Arousal is exploratory only. The nuisance set is the fMRIPrep motion parameters and derivatives, FD, CSF, WM (both from the derivatives of record), cosines and FD spikes. Sensitivity models add RVT⊗RRF and a 4-lag HRF-convolved PetCO2 basis (§6–§7). Only between-condition contrasts are interpreted; there is no rest baseline (§5). All inputs exist and are frozen in the repository (`FirstLevel/`, `Data/regressors/`, `Data/maia_scores.csv`), and the manifest carries the gates.

Group level is a one-sample t on fixed-effects-combined participant contrasts, with MAIA total as a covariate on the paced-vs-free contrast for H6. H3 instead uses FSL FLAME 1 on each participant's contrast and variance, combined over the runs that contain both cells, with random-field-theory inference. An exploratory graded-awareness model codes aware trials +confidence and unaware trials −confidence, so a low-confidence aware trial counts as less aware.

Connectivity (H7, §7) uses the median signal of the CAREN SAL, SMN and ECN networks from the unsmoothed preprocessed BOLD. Each run is band-passed at 0.01–0.1 Hz, and the data are regressed on the nuisance set and the Model A task regressors. Windows are epoch-locked, shifted 4 s and demeaned per epoch. Coupling is the Fisher-z correlation per condition.

## Q6. What are the statistical models and inference criteria?

§7–§8 in full. Whole-brain tests use cluster-forming p < .001 and cluster FWE p < .05, two-sided. ROI tests use small-volume FWE p < .05, one-sided in the predicted direction, Bonferroni-corrected across the ROIs tested per hypothesis. For H3 the FWE correction is by random field theory on the FLAME result. The masks come from the CAREN consensus network atlas (network × AAL-3 region): anterior insula and dACC from the salience network, and posterior insula, central operculum and sensorimotor cortex from the somatomotor network. A Harvard–Oxford version is reported as atlas sensitivity only. H6 is two-sided. H7 uses two-sided paired tests, Holm-corrected across two network pairs, with a global-signal interpretation gate. H5 is decided by the TOST grid. No uncorrected result will be reported as a finding.

Power for every registered N is tabulated in §4/§8. This includes the statement that H6 has 80% power only for |r| ≥ 0.40, and ≈ 44–59% at IEAT-sized correlations.

## Q7. Prior publications and involvement with this dataset

No publication uses this dataset. The design descends from the IEAT study (Farb, Zuo & Price, 2023), whose senior author leads this project. H1 is an explicit replication of that paper's pacing result in a new sample. H6 is a conceptual extension of its MAIA result. The behavioural and physiological analyses summarised in Q3 are unpublished and will be reported alongside the imaging results.

## Q8. Blinding, deviations, and timeline

No blinding is possible for the behavioural results (Q3). The imaging analyst has not seen any statistical image. Deviations after freezing will be logged in `Preregistration/deviations.md` (git-versioned) with date, rationale, and whether the affected result had been seen. Post-hoc analyses will be labelled exploratory without exception (§10–§11). The registration will be frozen before the first GLM is fitted, and the analysis code repository will be linked as an OSF component.

---

## Author sign-off

| Name | Role | Approves freezing this registration | Date |
|---|---|---|---|
| Norman A. S. Farb | PI (sole author) | ☐ | |

*Sole-author project. Freezing before first-level modelling is the commitment.*
