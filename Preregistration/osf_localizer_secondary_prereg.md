# OSF Secondary Data Preregistration: Intero2023 localizer (paste-ready answers)

**What this is.** These are answers to OSF's *Secondary Data Preregistration* template (van den Akker et al., 2021), in the template's field order. Paste each answer under the field of the same name. Required fields are marked *.

**How it relates to swc7p.**
- This is a **new registration linked to OSF swc7p**, not an update of it. It covers the same sample.
- swc7p registered the breathing change-detection runs; this registration covers the localizer runs, and their use to define regions tested in the breathing runs.
- File both registrations in the same OSF project.
- The plan of record is `Preregistration/intero2023_localizer_prereg.md` (v1.0) in https://github.com/Normega/MRI_Intero2023. Section references (§) below point into it.

**Before submitting.**
1. Fill in the item marked **[PI to confirm]**, under Prior publication.
2. The plan is frozen at commit d29136b (https://github.com/Normega/MRI_Intero2023/blob/d29136b/Preregistration/intero2023_localizer_prereg.md); cite that link in the registration.

---

## Metadata

### Title*
Neural correlates of attention to the breath: an fMRI localizer for interoceptive attention, and its use to test awareness of breathing changes (Intero2023)

### Description*
A secondary-data preregistration for the two localizer runs of the Intero2023 paced-breathing fMRI study.
- **Same sample as OSF swc7p.** swc7p registered the breathing change-detection runs.
- **The localizer contrasts** attention to the breath with attention to a matched visual circle, self-referential judgement, paced breath-matching and fixation.
- **Three aims:**
  - map the neural correlates of these attentional targets (activation and connectivity);
  - replicate the paced < free breathing effect (swc7p H1) in independent runs, and the IEAT interoceptive-awareness effect;
  - use participant-specific interoceptive-attention ROIs to test whether awareness of breathing changes in the change-detection runs is detectable in localizer-defined regions when it was not in anatomical masks.
- **Status:** no localizer imaging data have been analysed.

### Contributors*
Norman Farb (University of Toronto Mississauga)

### License*
CC-By Attribution 4.0 International

### Subject*
Social and Behavioral Sciences → Psychology → Cognitive Psychology; Neuroscience and Neurobiology → Cognitive Neuroscience

### Tags
interoception; breathing; attention; fMRI; insula; functional localizer; participant-specific ROI; MAIA; secondary data

---

## Study information

### Research questions*
- **RQ1.** Does attending to the breath, compared with attending to a matched external visual stimulus, change activity in the salience-network regions (anterior insula, dorsal anterior cingulate)?
- **RQ2.** Does paced breathing lower activity relative to free breathing in somatomotor regions, in runs independent of those used for swc7p H1?
- **RQ3.** Is the breath-vs-visual attention effect in anterior cingulate moderated by self-reported interoceptive awareness (MAIA)?
- **RQ4.** Can participant-specific interoceptive-attention regions be defined reliably from two localizer runs?
- **RQ5.** In those participant-specific regions, do aware breathing-change trials evoke more activity than unaware trials? Is that difference larger than in the rest of the same anatomical region?
- **RQ6.** In the participant-specific ACC region, does MAIA predict the response to paced vs free breathing?

### Hypotheses*
- **L1 (RQ1).** Intero − Extero ≠ 0 in the left anterior insula, right anterior insula and dACC masks. Non-directional.
  - Rationale: the two tasks share buttons, a moving circle and response rhythm, and differ in the attended target. Prior work (the IEAT) and swc7p H1 point to both increases and decreases.
- **L2 (RQ2).** Match − Intero < 0 in the posterior insula, central operculum and sensorimotor masks. Directional: paced breathing lowers activity. This replicates swc7p H1 in independent runs.
- **L3 (RQ3).** MAIA × (Intero − Extero) ≠ 0 in the IEAT ACC sphere. Non-directional. This is the direct form of the IEAT finding; swc7p H6 tested a conceptual extension.
- **L4 (RQ4).** Participant-specific ROIs pass a pre-specified reliability gate (§7):
  - the cross-run effect is above 0;
  - Dice ≥ .20;
  - within-participant cross-run spatial similarity exceeds between-participant similarity.
- **T1 (RQ5).** In participant-specific Intero > Extero ROIs within the anterior insula and dACC masks: aware > unaware (directional), **and** the aware − unaware effect is larger in the ROI than in the remainder of the same mask (directional).
- **T3b (RQ6).** MAIA × (paced − free) > 0 in the participant-specific ACC-sphere ROI. Directional, following swc7p H6.

---

## Data description

### Datasets used*
**Source.** The Intero2023 fMRI dataset (single site, University of Toronto, 2023–2024). Each session has five functional runs:
- 2 localizer runs (this registration);
- 2 breathing change-detection runs (swc7p);
- 1 authenticity run (not used).

Acquisition: 3 T, TR 2 s; about 196 volumes per localizer run. Preprocessing: fMRIPrep 24.1.1, MNI152NLin2009cAsym.

**Localizer.** Each run has 10 blocks of 30 s: two shuffled sets of five tasks. Each block is preceded by a 6 s instruction and a 2 s fixation.

| Task | What participants do |
|---|---|
| Intero | Track their breath with inhale/exhale button presses |
| Extero | Track a pacer circle's growing and shrinking with the same buttons |
| Match | Breathe with the circle while pressing |
| Self | Yes/no on trait words |
| Rest | Fixation |

**Other data used:**
- the respiration belt (Biopac), and PetCO2 where available;
- MAIA-2;
- the breathing-run first-level contrasts registered in swc7p.

### Data availability*
The dataset is not publicly available.

### Data access
Raw imaging and physiological data are held by the PI under the ethics approval (University of Toronto REB protocol 45079) and are not publicly shared. Analysis code, task event files, and derived non-identifying tables are public at https://github.com/Normega/MRI_Intero2023.

### Data identifiers
https://github.com/Normega/MRI_Intero2023 (code and derived tables); linked registration https://osf.io/swc7p.

### Access date*
- **Imaging and physiology:** collected 2023–2024; held by the PI since collection.
- **Localizer task logs and raw respiration files:** copied to the analysis server on 2026-11-27 for the quality checks listed under Prior knowledge.
- **Localizer BOLD:** preprocessed (fMRIPrep) in 2026, but not analysed; only headers and motion estimates have been read.

### Data collection procedures*
- **Recruitment and ethics.**
  - Community adult sample.
  - MRI safety screening by questionnaire before the session.
  - Approved by the University of Toronto Research Ethics Board (protocol 45079).
  - Sessions lasted 2 hours (30 min training, about 1 h in the scanner, 30 min questionnaires); participants were compensated $65.
  - Eligibility was MRI safety.
- **Task program.** The localizer was run with the PsychoPy task program `InteroPilot.4.51.py`.
- **Protocol versions.**
  - v1 (8 participants): the pacer was set to each participant's calibrated breathing rate (7–16 /min), and Self and Rest had no pre-task fixation.
  - v2: the pacer was fixed at 12 /min.
- **Known task-logging limitation.** In Extero, Match and Self the logged response is overwritten by the scanner trigger key, so accuracy and yes/no answers cannot be scored. Intero presses are recorded cleanly.

Full documentation: §1 and `Analysis/localizer_qc/README.md`.

### Codebook
The events files (`Data/events/*_task-localizer_*_events.tsv`) have these columns:
- onset, duration, trial_type (intero, extero, match, self, rest, instruction);
- block_index, protocol_version, n_presses, compliant.

The QC tables in `Analysis/localizer_qc/` are documented in its README.

---

## Variables

### Manipulated variables
Attention target / task, within participants: Intero, Extero, Match, Self and Rest (implicit baseline); 4 blocks each over two runs. The first-level contrasts are:

| # | Contrast | Meaning |
|---|---|---|
| C1 | Intero − Extero | Breath vs visual |
| C2 | Intero − Self | Breath vs self |
| C3 | min-t(C1, C2) | Conjunction |
| C4 | Match − Extero | Breath coupling |
| C5 | Match − Intero | Paced vs free |
| C6 | Self vs Rest | Self-reference |

From the breathing runs (swc7p):
- awareness of the pacer change: aware vs unaware, the C2 joint contrast;
- paced vs free breathing: C1.

### Measured variables*
- **Outcomes.** Voxelwise BOLD effects (percent signal change):
  - mask-level, in frozen a-priori masks: CAREN left and right anterior insula, dACC, posterior insula, central operculum, sensorimotor cortex; the IEAT ACC sphere;
  - fROI-level, in participant-specific ROIs (§7).
- **Predictors and covariates.**
  - MAIA-2 total (the IEAT-defined total, 0–5 metric, as in swc7p).
  - ΔRVT (Intero − Extero breathing volume, participant mean).
  - Extero entrainment: phase-locking of breathing to the pacer, minus a same-participant virtual-pacer null.
  - Protocol version (v1/v2).
- **Exploratory physiological measures** (§4):
  - breath-tracking accuracy, as the phase-locking of Intero presses to belt breath phase;
  - pacer adherence (phase-locking) in Match and in the paced epochs of the change-detection runs;
  - breathing rate, depth and RVT per condition.

### Missing data
- **MAIA** is missing for some participants: 46 of the 48 Intero-eligible; 42 of the 44 fROI-eligible. Those participants are excluded from L3 and T3b only.
- **Respiration:**
  - usable for rate and depth in 71 runs (41 participants);
  - usable for breath-level timing in 70 runs (40 participants);
  - RVT on both runs in about 29 of the 44 fROI-eligible.

  Physiological covariates and the +RVT co-report use available data and are reported as subset results.
- No imputation.

### Unit of analysis
**Participants (fixed effects over runs, random effects at group level).**

**Run exclusions.**
- Two runs in which the task program stopped mid-run (sub-058 run-2, sub-114 run-2).
- The swc7p run-level motion gate: mean FD > 0.5 mm, or > 25% of volumes with FD > 0.5 mm. This excludes 2 runs.
- Result: **102 runs, 54 participants.**

**Intero compliance.**
- A block is compliant with ≥ 4 presses.
- A participant needs ≥ 3 of 4 compliant blocks (2 of 2 with one usable run).

**Expected N.**

| Analysis | N |
|---|---|
| Contrasts without Intero | 54 |
| Intero contrasts (L1) | 48 |
| L3 | 46 |
| fROI-eligible (two usable runs; ≥ 1 compliant Intero block per run) | 44 |
| T1 (fROI ∩ the swc7p H3 set) | 40 |
| T3b | about 40 (fROI ∩ swc7p H6 set) |

### Statistical outliers
No participant is removed as a statistical outlier. FLAME 1 down-weights high-variance participants. Each participant's share of the total weight is reported, as in swc7p.

### Sampling weights
Not applicable.

---

## Knowledge of data

### Prior publication/dissemination
- OSF swc7p (registration, and its approved update) for the breathing change-detection runs of this dataset.
- The breathing-run results (registered tests and exploratory analyses E1–E10) are summarised in the linked repository (`Analysis/RESULTS_SUMMARY.md`); not yet published.
- **[PI to confirm]** whether the companion behavioural manuscript (the BCAT preprint) includes any of this sample.
- No localizer analysis has been presented or published.

### Prior knowledge*
**Norman Farb** designed the study and knows the following.

*Breathing-run results* (swc7p):
- **Supported:** H1, paced < free in somatomotor and posterior insula masks; H6, MAIA × paced − free in the IEAT ACC sphere.
- **Not supported:** awareness (aware vs unaware), confidence, step vs ramp, and connectivity.
- **Exploratory:** further awareness tests in the anterior insula and dACC masks were null, mostly with BF01 about 5–6.

*Localizer, examined before this registration*, for quality control and to set analysis parameters:
1. **Task logs:** timing, block order, attempt selection, the v1 timing structure, two runs where the task stopped, Intero press counts and compliance.
2. **Motion and headers:** framewise displacement (motion gate) and NIfTI headers.
3. **Design efficiency** from the design matrices. This set the 200 s high-pass: the 128 s cutoff lost 49–55% of contrast efficiency, beyond the pre-set 20% trigger.
4. **Raw respiration belt:**
   - quality and polarity;
   - per-condition breathing: Intero breathing is deeper than Extero (RVT dz .88); breathing entrains to the pacer in Extero (phase-locking .57 vs null .32);
   - breath-tracking accuracy: its distribution, chance level and run-to-run reliability. A fixed-window score was at chance because presses lead the belt by about 0.85 s, so the phase-locking measure was adopted.
5. **Registered plan:** one localizer-derived number already appears in the swc7p plan: PetCO2–global-signal coupling of +0.126.

*Not seen:* no localizer BOLD map, ROI value or statistic has been computed. No localizer behavioural or physiological measure has been related to MAIA, breathing-run behaviour, demographics or brain data.

*Analysis process:* analysis plans were drafted and checked with AI-assisted analysis tools, under the PI's direction. Two independent plans were written blind to localizer BOLD and merged.

---

## Analyses

### Statistical models*
**First level** (§3).
- Joint two-run GLM, the same as swc7p:
  - PSC per run, AR(1), swc7p nuisance set and FD spike regressors;
  - SPM canonical HRF.
- Regressors: four task boxcars (Rest is the implicit baseline), an instruction-plus-fixation regressor, and a regressor of no interest for non-compliant Intero blocks.
- **High-pass: 200 s** cosine set; fMRIPrep's 128 s cosines are dropped.
- Smoothing: 6 mm for group maps; 4 mm for participant-specific ROIs.

**Group** (§5).
- FLAME 1, one-sample (L1, L2) or with demeaned MAIA (L3).
- GRF voxel-level FWE within each frozen mask.
- Whole-brain cluster FWE at z > 3.09, both tails.

**Participant-specific ROIs** (§7).
- Within each mask, each participant's top 10% of voxels by C1 t (two-run model, 4 mm); minimum 50 voxels.
- Reliability gate (L4), on localizer data only.
- If the within- vs between-participant similarity test fails for a mask, the fallback is a group-constrained ROI: the top 10% of the leave-one-subject-out group C1 map.

**Transfer** (§8).
- Breathing-run contrasts are refitted at 4 mm, with no other change to the swc7p model.
- Extract the ROI mean, the mask complement mean and the mask mean, in the same participants.
- T1: FLAME-weighted one-sample test, as swc7p H3, plus a paired ROI − complement test (bootstrap CI on the dz difference).
- T3b: MAIA covariate, as swc7p H6.

### Effect size
- No effect size is predicted for L1–L3.
- For T1, the smallest effect of interest is the minimum detectable effect at N = 40: dz ≈ .45.

### Statistical power
- **T1:** N = 40 gives 80% power for dz ≥ .45 at α = .05 two-sided (equivalently .025 one-sided); with the ×3 correction, dz ≥ .48. The swc7p awareness effect in anatomical masks was dz ≤ .22, so T1 is powered only for a substantial gain.
- **L1:** N = 48 gives 80% power for dz ≥ .41 (two-sided).

### Inference criteria
- **L1:** two-sided, Bonferroni ×3.
- **L2:** one-sided, ×3.
- **L3:** two-sided.
- **T1:** one-sided, ×3, and the ROI − complement CI excludes 0.
- **T3b:** one-sided.
- **Breathing robustness for L1.** In participants with RVT on both runs, the +RVT model is fitted. L1 is "breathing-robust" if the +RVT effect keeps its sign and ≥ 50% of the base effect in the same participants; otherwise it is reported as breathing-dependent.
- **"More sensitive than the swc7p analyses"** requires T1 support **and** ROI − complement > 0. Significance in the ROI alone is not enough.
- JZS Bayes factors (r = .707) are reported throughout. BF01 > 3 in the participant ROIs is interpreted as evidence against the awareness effect there.

### Assumption violation / model non-convergence
- FLAME is robust to unequal participant variance. Permutation voxel max-t (10,000 sign-flips) is reported as a check.
- If FLAME fails in a mask, the permutation result is primary for that mask.
- If an fROI cannot be formed (fewer than 50 voxels in the mask), the group-constrained fallback is used for that participant.

### Reliability and robustness testing
- The reliability gate (L4).
- High-pass at 128 s and 320 s.
- +RVT, and +PetCO2 where available.
- Unsmoothed fROI extraction.
- fROI size 5% and 20%.
- Excluding protocol v1 participants.
- Covariates ΔRVT and Extero entrainment.

### Exploratory analysis
All of the following are reported as exploratory:
- participant-specific ROIs from Intero > Self, and from the conjunction with Intero > Extero;
- network connectivity (CAREN networks, as swc7p H7) and gPPI by condition (§6);
- transfer of confidence and change − No Change;
- per-epoch pattern expression of each participant's C1 map in the breathing runs (LSS betas; mixed model with confidence and awareness);
- decoding, RSA and insula gradients (§10);
- breath-tracking phase-locking and pacer adherence in both tasks, and their relations to MAIA, breathing-run behaviour and brain measures (§4).
