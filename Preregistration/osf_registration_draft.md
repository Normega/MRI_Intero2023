# OSF Registration Draft — Neural correlates of respiratory change awareness

**Registration type:** Preregistration of a secondary data analysis (template after Mertens & Krypotos, 2019)
**Companion document:** `intero2023_imaging_prereg.md` (full analysis plan; to be attached to the registration). Section references (§) below point into it.
**Status:** DRAFT — awaiting co-author sign-off. Freeze on OSF before any first-level model is fitted.

---

## Q1. What are the hypotheses under investigation?

Six hypotheses, stated in full with falsification criteria in §3 of the attached plan:

- **H1** — Paced (exogenously controlled) breathing produces relative *deactivation* of the insula/operculum pathway and somatomotor cortex compared with free breathing, replicating the Active Matching result of Farb, Zuo & Price (2023, eNeuro; "IEAT"). Direction committed to the replication; the volitional-control account predicts the opposite and a reversed result will be reported as a failed replication.
- **H2 (primary)** — Trial-by-trial confidence on change trials correlates positively with activity in bilateral anterior insula and dACC (N = 43).
- **H3 (secondary)** — Aware (correct-direction) > unaware (responded "No Change") change trials in the same three ROIs (N = 21; ROI-restricted; a null alongside a positive H2 is attributed to power in advance).
- **H4** — Salience (step vs ramp) effects fall on salience-network regions, consistent with the behavioural finding that salience affects noticing, not direction discrimination.
- **H5** — Physiological correction attenuates but does not abolish the awareness effects; adjudicated by a preregistered TOST equivalence criterion (bounds dz = 0.65 / 0.59; interpretation grid fixed in §12.4).
- **H6** — Higher MAIA interoceptive sensibility predicts spared ACC activity during paced breathing, tested one-sided in a 10 mm sphere at the IEAT peak (MNI −4, +12, +36); N = 43.

## Q2. What data will be used, and how were they generated?

MRI_Intero2023: a single-session fMRI study (3T Siemens, TR = 2 s, fMRIPrep 24.1.1 derivatives) in which participants performed a paced-breathing change-detection task (faster / slower / no change × step / ramp salience), free-breathing periods, and rating windows, with concurrent respiration-belt and expired-CO2 (capnograph) recording plus scanner PMU physiology. Data collection is complete (2023-09 – 2024-06). Nothing further will be collected.

## Q3. Have the data been accessed or analysed? (full disclosure)

**Yes — everything except the functional imaging.** This is the core disclosure and is stated identically in §0 of the attached plan:

*Fully analysed before this registration:* task logs and event timing; all trial-level behaviour (three-category identification, ordinal SDT, confidence/metacognition); respiration-belt processing including pacer adherence; capnography and PetCO2 extraction; RVT; response-function validation; all QC and exclusion decisions.

*Not examined in any form:* the functional images beyond fMRIPrep's automated preprocessing. No first-level model has been fitted, no contrast image exists, no statistical map has been computed or viewed.

*Single documented exception:* fMRIPrep's `global_signal` confound column was used to verify that the physiological regressors carry the physiologically correct sign (§6). fMRIPrep's `framewise_displacement` column was used to apply the preregistered motion gate. No task regressor has been related to any imaging measure.

*Consequence, stated plainly:* H2–H4 are informed by known behavioural results. They are untested predictions about brain data, but they are not blind hypotheses. H1 and H6 are directional replications of an independent published sample (IEAT).

## Q4. How was the sample determined, and what are the exclusion rules?

The sample is the full collected cohort minus preregistered exclusions, all already applied and enumerated in §4: pilot/QA IDs; two participants with shifted run labels in the derivatives (detected by volume-count reconciliation, 242/244 exact); runs failing the motion gate (mean FD > 0.5 mm or > 25% volumes above 0.5 mm — applied 2026-08-14 from the FD column only, removing 11/95 runs and 5 participants). **Final Ns: 43 (H1/H2/H6), 41 (H4), 21 (H3), 22/26 (S1/S2).** A registered check shows the motion gate does not select against pacer-adherent breathers (fd–adherence r = −0.003).

## Q5. Which variables enter the analysis?

First-level: seven event regressors (aware / unaware / wrong-direction / no-change / free / rating / sub-rating), confidence and arousal as parametric modulators, fMRIPrep nuisance set plus FD spikes; sensitivity models add RVT⊗RRF and a 4-lag HRF-convolved PetCO2 basis (§6–§7). All inputs exist and are frozen in the repository (`FirstLevel/`, `Data/regressors/`, `Data/maia_scores.csv`); the manifest carries the gates. Group level: one-sample t on fixed-effects-combined participant contrasts; MAIA total as covariate for H6.

## Q6. What are the statistical models and inference criteria?

§7–§8 in full. Whole brain: cluster-forming p < .001, cluster FWE p < .05. ROI: small-volume FWE p < .05, Bonferroni across the ROIs tested per hypothesis. H6 one-sided. H5 by the TOST grid. No uncorrected result will be reported as a finding. Power for every registered N is tabulated in §4/§8, including the statement that H6 sits at ≈ 67% power for the IEAT-sized effect.

## Q7. Prior publications and involvement with this dataset

No publication uses this dataset. The design descends from the IEAT study (Farb, Zuo & Price, 2023), whose senior author leads this project; H1 and H6 are explicit replications of that paper's results in a new sample. The behavioural and physiological analyses summarised in Q3 are unpublished and will be reported alongside the imaging results.

## Q8. Blinding, deviations, and timeline

No blinding is possible for the behavioural results (Q3). The imaging analyst has not seen any statistical image. Deviations after freezing will be logged in `Preregistration/deviations.md` (git-versioned) with date, rationale, and whether the affected result had been seen; post-hoc analyses will be labelled exploratory without exception (§10–§11). Registration will be frozen before the first GLM is fitted; the analysis code repository will be linked as an OSF component.

---

## Co-author sign-off

| Name | Role | Approves freezing this registration | Date |
|---|---|---|---|
| Norman A. S. Farb | PI (sole author) | ☐ | |

*Sole-author project. Freezing before first-level modelling is the commitment.*
