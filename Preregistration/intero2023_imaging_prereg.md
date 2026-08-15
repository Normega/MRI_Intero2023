# Preregistration — Neural correlates of respiratory change awareness (MRI_Intero2023)

**Status:** draft for discussion
**Version:** 0.1
**Date:** 2026-08-13
**Analysis code:** `Repo/R/`, `Repo/Scripts/` (all preprocessing and behavioural code already written and versioned)

---

## 0. Transparency: what is already known

This document is written **after** the behavioural and physiological data were fully processed and analysed, and **before** any functional imaging analysis has been run. That asymmetry is the point of registering now, but it has to be stated plainly rather than glossed:

**Already analysed, and used to shape the hypotheses below:**

- All task logs, event timing, and trial-level behaviour (`intero2023_behaviour.R`)
- Respiration belt processing, pacer adherence, per-volume RVT and rate
- Capnography, PetCO2 extraction, response-function convolution
- All quality control, exclusions, and the analysable sample reported in §4

**Not yet examined in any form:**

- Every functional image beyond fMRIPrep's automated preprocessing and its confound time series
- No first-level model has been fitted; no contrast image exists; no group map has been computed or viewed

The global signal *has* been used, for one narrow purpose: validating that the physiological regressors have the physiologically correct sign (§6). No task regressor has been related to any imaging measure.

Because the behavioural results are known, hypotheses H2–H4 below are **informed by observed behaviour**. They are predictions about *brain* data, which remain untested, but they are not blind to the behavioural pattern. Anything genuinely exploratory is marked as such in §10 and will be reported as exploratory regardless of outcome.

---

## 1. Study aims

The task manipulates the rate of paced breathing (faster / slower / no change, delivered as either an abrupt step or a gradual ramp) and asks participants to report what their breathing did. It also includes unpaced free-breathing periods. The aims:

1. Identify the neural response to **voluntary paced breathing** relative to free breathing.
2. Identify the correlates of **awareness of a respiratory change** — the central question.
3. Determine whether those correlates survive **correction for physiological (vascular) confounds**, which in this paradigm are partly manipulated by the task itself.
4. **Replicate key results of the direct predecessor study**, the IEAT (Farb, Zuo & Price, 2023, eNeuro 10(6):ENEURO.0088-23.2023), in a larger sample with an added awareness manipulation.

### 1b. Relation to the predecessor study (IEAT)

The IEAT (N = 22, 44 scans) compared attention to the breath against a difficulty-matched visual task, and included an **Active Matching** condition — breathing paced to an external stimulus — that is the direct ancestor of this study's paced epochs. Its results constrain several predictions here:

- **Interoceptive attention produced widespread cortical deactivation**, not activation of interoceptive cortex, with preserved performance.
- **Exogenous pacing (Active Matching vs Active Interoception) drove additional deactivation along the insula/operculum pathway, somatomotor cortex, ventral occipital cortex, and cerebellum.** The right anterior insula was implicated *only* by this endogenous→exogenous contrast, not by interoception vs exteroception.
- **Self-reported interoceptive sensibility (MAIA) predicted sparing from deactivation in the anterior cingulate** (largest cluster k = 1205, peak MNI −4, +12, +36, BA 24) and left-lateralised language regions; ACC–DAN connectivity increased during interoception.
- Its discussion draws a distinction directly relevant here: **anterior insula tracks liminal detection accuracy** (heartbeat-style paradigms), while **attention to salient respiratory signals implicates the ACC**. The present task's respiratory changes are salient and supraliminal (direction accuracy 0.78), which strengthens the a priori case for the ACC alongside the insula in the awareness hypotheses.

Correspondence between designs: this study's **paced epochs ≈ Active Matching** (exogenous respiratory control); **free epochs ≈ Active Interoception** (endogenous rhythm under interoceptive monitoring). The present design adds what the IEAT could not ask: trial-by-trial *awareness* of an imposed respiratory change.

MAIA-2 items are present in this study's screener battery, so the sensibility moderation is replicable (H6).

---

## 2. Key behavioural findings that constrain the imaging hypotheses

These are established, not predicted:

| Finding | Value |
|---|---|
| Three-category identification accuracy (chance = 1/3) | **0.525**, κ = 0.271, χ²(4) = 141.7, p = 1e-29 |
| Ordinal SDT sensitivity | **0.698** probit units per stimulus step, z = 13.2 |
| Direction accuracy on true-change trials | **0.778**; per-participant median 0.833, 85% above chance |
| Accuracy on No Change trials | 0.301 — **at chance** (p = 0.79) |
| P(respond "changed") | ~0.70 regardless of stimulus (strong liberal bias) |
| Confidence: "No Change" response | **lowest** confidence in the set (5.03) |
| Confidence–accuracy (type-2 AUC) | **0.661**, 74% of participants above 0.5, p = 0.0001 |
| Salience (step vs ramp) effect on *noticing* | p = 0.043 |
| Salience effect on *direction accuracy* | p = 0.26 (n.s.) |
| Arousal across conditions | flat (4.15–4.27) |

Two consequences carry directly into the imaging model:

- **Sensitivity is real and individually reliable**, so an awareness contrast is a split on signal rather than on noise.
- **"No Change" functions as "I don't know"**, not as a positive judgement that nothing happened. It is the lowest-confidence response and is at chance when correct. It will therefore be modelled as a low-awareness state, never as a correct rejection.

---

## 3. Hypotheses

Regions are specified a priori as ROIs in §8. Directional predictions are one-sided in the sense stated; whole-brain tests remain two-sided.

**H1 — Exogenous respiratory control (replication of IEAT Active Matching).**
Two competing directions exist for the paced vs free contrast, and this document commits to the replication direction. The naive volitional-control account predicts paced > free activation in SMA, sensorimotor cortex, and insula. The IEAT found the opposite: exogenous pacing produced greater **deactivation** along the insula/operculum pathway, somatomotor cortex, ventral occipital cortex, and cerebellum relative to endogenous breathing.
*Registered prediction:* **paced < free** (relative deactivation) in the insula/operculum and somatomotor ROIs, replicating the IEAT.
*The test is diagnostic either way:* a paced > free effect in these regions would contradict the predecessor study in a larger sample and will be reported as a failed replication, not reframed. This hypothesis was revised from the volitional-control direction on 2026-08-14, on reading the IEAT results and **before any imaging analysis**.

**H2 — Awareness, parametric (PRIMARY).**
Trial-by-trial confidence on change trials will correlate positively with activity in **bilateral anterior insula** and **dorsal anterior cingulate**.
*Rationale:* confidence is metacognitively valid here (type-2 AUC 0.661), and these regions are the canonical substrate for interoceptive salience and confidence. This is the best-powered awareness test (N = 48).
*Falsified if* no positive confidence effect survives small-volume correction in either anterior insula ROI.

**H3 — Awareness, categorical (SECONDARY, N = 22).**
Change trials answered with the correct direction ("aware") > change trials answered "No Change" ("unaware") in bilateral anterior insula and dACC.
This is the more interpretable contrast but is available in fewer participants (§4). It is registered as a confirmatory test of the same construct as H2, with ROI-restricted inference only. A null here with H2 positive will be attributed to power, not taken as disconfirming.

**H4 — Salience acts on detection.**
High-salience (step) > low-salience (ramp) change trials will elicit greater activity in the salience network (anterior insula, dACC).
*Rationale:* behaviourally, salience affected *noticing* (p = 0.043) but not *direction discrimination* (p = 0.26), so its neural effect is predicted to fall on detection-related rather than discrimination-related regions.
*Falsified if* the effect appears instead in sensory/discrimination regions with no salience-network involvement.

**H5 — Physiological correction attenuates but does not abolish awareness effects.**
Adding physiological regressors (§6) will reduce cluster extent in H1–H4 but the anterior insula and dACC effects in H2 will survive.
This is directional and consequential: PetCO2 is a **mediator** in this design, not a pure nuisance — the pacing manipulation changes end-tidal CO2 by construction. Regressing it out therefore removes part of the effect of interest for H1 in particular. Both corrected and uncorrected models will be reported side by side, and neither is designated "the" result for H1. Precedent: the IEAT found RVT associated with widespread deactivation, and that controlling respiration rate and RVT at both levels attenuated but did not eliminate its interoception deactivation effect.

**H6 — MAIA moderation of the ACC (direct replication, N = 47).**
Higher MAIA interoceptive sensibility will predict **greater (spared) ACC activity** in the paced-related contrasts, replicating the IEAT's central individual-differences result (peak −4, +12, +36, k = 1205).
*Test:* MAIA total (scored as in the IEAT: mean of the six subscales omitting Not-Worrying and Not-Distracting; scores in `Repo/Data/maia_scores.csv`, scoring code `intero2023_maia.R`) as a group-level covariate on the H1 contrast; small-volume inference in a 10 mm sphere at the IEAT peak (−4, +12, +36), plus the atlas dACC ROI as a secondary mask. **One-sided** in the replication direction, as the direction is fully specified by the prior study.
*Coverage and power:* complete MAIA for 47 of 48 imaging participants (M = 3.31, SD = 0.92 on the canonical 0–5 metric; all 37 items answered by every scored respondent). At N = 47, one-sided α = .05: 80% power for r ≥ 0.35; at the IEAT's observed total-score correlation (r ≈ 0.32), power ≈ 72%. Adequately powered for the replication-sized effect; this is stated so a marginal result is read against it.
*Rationale:* the strongest form of replication available — same construct, same region, a priori coordinate from an independent sample. MAIA subscale analyses (Attention Regulation vs Noticing dissociation reported in the IEAT) are exploratory (§10).
*Administration note (registered):* this study administered the 37-item MAIA-2 with response boxes coded 1–6; scores are converted to the canonical 0–5 metric before use. The IEAT used the original MAIA; the six retained subscales are common to both.

---

## 4. Sample, exclusions, and analysable N

Exclusions are already applied and are not contingent on imaging results.

**Excluded a priori:**

- Pilot/QA identifiers P00, P01, P50, P999
- **sub-082 and sub-104** — run labels are shifted in the derivatives tree (the file named `task-breath1` holds 172 volumes where the breathing slot expects 242; 172 is the authenticity length). Every other run agrees with `inventory.tsv` exactly (242/244). Excluded rather than guessed.
- Runs with no usable event timing or no trigger-locked physiology

**Analysable sample (established, not projected):**

| | runs | participants |
|---|---|---|
| Breathing runs with events | 105 | — |
| + usable imaging | **95** | **48** |
| + RVT regressor | 64 | 30 (all runs) |
| + PetCO2 regressor | 51 | 24 (all runs) |
| + both physiological regressors | 43 | 20 (all runs) |
| Verified pacer adherence (≥1 run) | — | 24 |

**Estimability per contrast** (≥3 trials per cell per participant):

| Contrast | N |
|---|---|
| C1 paced vs free | **48** |
| C2 aware vs unaware (categorical, H3) | **22** |
| C3 confidence parametric (H2) | **48** |
| C4 wrong-direction vs aware | 13 — **not analysed as a planned contrast** |
| C5 salience high vs low (H4) | **47** |

**Additional motion exclusion, applied before modelling:** a run is dropped if mean framewise displacement > 0.5 mm or if > 25% of its volumes exceed FD 0.5 mm. A participant is dropped if this leaves fewer than one usable breathing run. This criterion is fixed now and will not be adjusted after seeing results.

**The gate has now been applied** (2026-08-14, using fMRIPrep's FD column only — no BOLD signal or contrast has been examined). It removes 11 of 95 runs and 5 participants entirely (sub-058, -062, -080, -085, -094); two more lose one of two runs. Gated runs are unambiguous (median 81 spike volumes vs 3 in retained runs; median FD 0.52 vs 0.18 mm). **Post-gate: 84 runs, 43 participants** (`Repo/FirstLevel/firstlevel_manifest.csv`).

*Registered check on a paradigm-specific selection risk:* deep paced breathing could inflate FD and the gate could then preferentially remove adherent breathers. It does not: fd–adherence correlation across belt-usable runs r = −0.003 (n = 81), and gated runs show *lower* pacer adherence (median 0.230 vs 0.335).

**Final registered Ns (post-gate):**

| Contrast | N |
|---|---|
| C1 paced vs free / C3 confidence parametric | **43** |
| C2 aware vs unaware (H3) | **21** |
| C5 salience (H4) | **41** |
| H6 MAIA moderation (all 43 have complete MAIA) | **43** |
| S1 ±PetCO2 (basis on all retained runs) | **22** — TOST bound dz = 0.65 (80%) |
| S2 ±RVT (regressor on all retained runs) | **26** — TOST bound dz = 0.59 (80%) |

H6 power updates to: 80% for r ≥ 0.37 one-sided; ≈ 67% at the IEAT's r ≈ 0.32.

---

## 5. Design and event model

Per breathing run, the following events are modelled (already extracted and validated against a hand-checked participant, 80 assertions in `Repo/Tests/`):

- **Paced breathing epochs** (~20 s), one per change trial, split by condition
- **Free breathing epochs** (~20 s), modelled as a separate trial type
- **Rating window**, modelled explicitly rather than left in the implicit baseline
- **Valence/sub-rating events** as nuisance

Trial counts per participant (median): 7 aware-correct, 2 unaware, 1 wrong-direction, 2 No Change, 14 free. The design is short — 14 paced trials per participant — which is why the parametric model (H2) is primary and the categorical split (H3) is secondary.

**Wrong-direction trials** (median 1 per participant) are modelled as a **separate regressor of no interest** so they do not contaminate the aware or unaware regressors. They are not contrasted.

**No Change trials** (median 2) are modelled as their own regressor and are not treated as correct rejections.

---

## 6. Physiological covariates

Two independent physiological measures are available, both already built and validated:

**RVT** (respiration volume per time, Birn et al. 2006), convolved with the respiration response function (Birn et al. 2008), z-scored within run.
Validated: RVT ⊗ RRF predicts the global signal at **r = +0.223, positive in 86% of runs**, at zero lag with no lag selection. Available for 153 runs / 41 participants overall.

**PetCO2**, as an HRF-convolved basis at lags 0/10/20/30 s.
The lag is **not identifiable** in these data — band-passed group lag profiles are correctly signed and cross zero at ~9 s but plateau from 18 s to 40 s without turning over, so a single fitted lag would be false precision. The 4-lag basis gives median multiple R² = 0.125 against the band-passed global signal. Available for 112 runs / 28 participants.

**Registered note on kernel sign.** The RRF is the correct kernel for RVT and the *wrong* kernel for PetCO2: it was derived as the impulse response from RVT to BOLD, and RVT moves opposite to end-tidal CO2. Its signed area is ≈ −16. PetCO2 will therefore never be convolved with the RRF in any reported model.

**Registered note on mediation.** In the breathing runs specifically, fixed-lag PetCO2–global coupling is ~0 (r = −0.004 at 20 s) while spontaneous-breathing localizer runs give +0.126. This is expected: pacing manipulates PetCO2 by design, so in breathing runs PetCO2 carries task variance. This is the reason H5 is framed as attenuation rather than as a purity check.

---

## 7. First-level model

One GLM per run, in MNI152NLin2009cAsym space (fMRIPrep 24.1.1 derivatives, already generated).

**Regressors of interest** (boxcar over epoch duration, convolved with canonical HRF + temporal derivative):

1. Paced breathing — aware (correct direction)
2. Paced breathing — unaware (responded "No Change")
3. Paced breathing — wrong direction *(no interest)*
4. Paced breathing — No Change trials
5. Free breathing
6. Rating window
7. Valence/sub-rating events *(no interest)*

**Parametric modulators**, mean-centred, applied to paced breathing epochs:

- **Confidence** (1–7) — primary, for H2
- **Arousal** (1–7) — secondary, reported but not a primary test given its flat behavioural profile

**Nuisance regressors (base model, all participants):**

- 6 rigid-body motion parameters + temporal derivatives
- Framewise displacement
- CSF and white matter mean signals
- Cosine drift terms (fMRIPrep, 128 s high-pass)
- Spike regressors for volumes with FD > 0.5 mm

**Physiological models (added in sensitivity analyses only, §9):**

- `+RVT`: RVT ⊗ RRF (1 regressor)
- `+PetCO2`: 4-lag HRF basis (4 regressors)
- `+both`

The global signal is **not** included as a nuisance regressor in any model.

---

## 8. Group model, ROIs, and inference

**Group model:** one-sample t-test on first-level contrast images, participant as the unit. Where a participant contributes two runs, run-level contrasts are combined at a fixed-effects intermediate level before entering the group model.

**Why fixed-effects combination (decided 2026-08-14, before any imaging analysis):**

- *Run-as-replication is not available for the awareness contrast.* Within a single run the cells collapse (median 3 aware / 1 unaware per run); requiring ≥3 trials per cell in both runs leaves exactly 1 participant, versus 25 pooled. Pooling is the condition for H3's estimability, not a trade-off against replication.
- *The runs are exchangeable.* Direction accuracy 0.800 vs 0.754 (p = 0.22), confidence 5.51 vs 5.31, P(say change) 0.726 vs 0.713 — a slight non-significant run-2 softening, no heterogeneity that would make an average misleading. Each run keeps its own GLM and baseline, so only the contrast is averaged.
- *Pooling materially improves per-participant reliability.* Cross-run individual differences are stable (direction accuracy r = 0.688, confidence r = 0.736, n = 42/51); Spearman-Brown gives the pooled estimate reliability ≈ 0.82 vs ≈ 0.69 for one run.
- *Run as a random factor is degenerate with two levels* and is not used.

**Registered reliability diagnostic (descriptive, not a gate):** for C1 (paced vs free) and C3 (confidence parametric), run-1 vs run-2 agreement of first-level contrast estimates — spatial correlation of contrast maps within participants, and cross-participant ICC of ROI means. Reported alongside the behavioural benchmark (r ≈ 0.7). H3 is excluded from this diagnostic because it is estimable in both runs for only 1 participant; this is stated now so its absence is not read as an omission later.

**A priori ROIs** (Harvard–Oxford probabilistic atlas, thresholded at 50%, bilateral where applicable):

1. Anterior insula (primary for H2/H3, with dACC)
2. Posterior insula
3. Dorsal anterior cingulate / mid-cingulate (primary for H2/H3, with anterior insula)
4. Supplementary motor area
5. Central operculum / sensorimotor cortex (H1)
6. **IEAT ACC replication sphere:** 10 mm radius at MNI (−4, +12, +36), the peak of the IEAT's MAIA-moderated ACC cluster (H6)

**ROI restriction for H3** (resolved 2026-08-14): H3 is tested in three masks only — left anterior insula, right anterior insula, and dACC (Bonferroni ×3). The dACC is retained alongside the insula on two grounds: the IEAT found the ACC, not the insula, to be the region tied to breath-focused interoceptive awareness; and its liminal-vs-salient distinction predicts ACC involvement for exactly the kind of salient respiratory signal this task uses. Narrowing to anterior insula alone would have encoded the heartbeat-detection literature's prediction while discarding the direct predecessor's.

**Smoothing** (resolved 2026-08-14): 6 mm FWHM, matching the IEAT for replication comparability.

**Thresholds:**

- **Whole brain:** cluster-forming p < 0.001 uncorrected, cluster-level FWE p < 0.05
- **ROI:** small-volume FWE p < 0.05 within each ROI, Bonferroni-corrected across the number of ROIs tested for that hypothesis
- No uncorrected results will be reported as findings; any uncorrected map shown will be labelled descriptive

**Power.** Smallest detectable Cohen's d at 80% power (one-sample t):

| N | α = 0.05 (ROI) | α = 0.001 (cluster-forming) |
|---|---|---|
| 48 (H1, H2, H4) | 0.41 | 0.63 |
| 24 (±PetCO2) | 0.60 | 0.96 |
| 22 (H3) | 0.63 | 1.01 |

H3 at N = 22 requires a large effect (d ≈ 1.0) at whole-brain thresholds, which is why it is ROI-restricted and secondary. This is stated in advance so a null there is not reinterpreted afterwards.

---

## 9. Registered sensitivity analyses

Each is fully specified now and will be reported whatever the outcome.

**S1 — ±PetCO2 (the primary methodological question, N = 24).**
Fit H1–H4 with and without the PetCO2 basis in the 24 participants who have it on all runs. Compare (a) peak z and cluster extent in each ROI, (b) whether each hypothesis's conclusion changes. Because this is a within-participant model comparison, it is better powered than the between-group comparison originally considered.

**S2 — ±RVT (N = 30).** As S1, with the RVT regressor.

**S3 — CO2-available vs CO2-unavailable participants.** Test whether the 24 participants with usable capnography differ from the remainder on the uncorrected H2 contrast. **Registered caveat:** CO2 availability is near-collinear with cohort and protocol version (the capnograph outage was an unbroken date block from 2023-10-12 to 2024-01-18), so this comparison **cannot isolate CO2** and will be reported as a cohort check only, not as evidence about physiological correction.

**S4 — Adherence conditioning (N = 24).** Repeat H2 restricted to participants with verified pacer adherence. Adherence is a run-level property (per-epoch classification is underpowered at AUC 0.66–0.69), and this is registered as robustness, not as the primary analysis.

**S5 — Protocol version.** v1 (n = 221 trials) and v2 (n = 1218 trials) differ in rating structure. Include protocol version as a group-level covariate; report whether H2 conclusions change.

---

## 10. Explicitly exploratory (will be labelled as such)

- Any brainstem or thalamic effect (coverage and physiological noise are unfavourable at 2 mm with this protocol)
- Free-breathing epochs related to measured respiratory change — the per-epoch respiration estimate has a noise floor (sd of log period ratio = 0.399) comparable to the imposed signal, so this is underpowered by construction
- Arousal as a parametric modulator
- Any relationship between neural measures and questionnaire/screener variables beyond the registered MAIA test (H6); MAIA subscale dissociations (Attention Regulation vs Noticing, per the IEAT) are exploratory
- Direction-specific effects (Faster vs Slower), given the behavioural asymmetry (0.614 vs 0.511)
- ACC-seeded PPI connectivity with the dorsal attention network during paced vs free breathing (conceptual extension of the IEAT's H4)

---

## 11. Deviations protocol

Any departure from this document will be recorded in `Repo/Preregistration/deviations.md` with the date, the reason, and whether it was made before or after the affected result was seen. Analyses added after seeing imaging results will be reported as exploratory without exception.

---

## 12. Open decisions for co-author input

1. ~~**Combining runs.**~~ **Resolved 2026-08-14:** fixed-effects within participant, with a registered descriptive run-reliability diagnostic. Rationale in §8.
2. ~~**Smoothing kernel.**~~ **Resolved 2026-08-14:** 6 mm FWHM, matching the IEAT.
3. ~~**H3 ROI restriction.**~~ **Resolved 2026-08-14:** bilateral anterior insula + dACC (three masks, Bonferroni ×3). The dACC is retained on the IEAT's evidence; see §8.
4. ~~**Formal equivalence test for H5.**~~ **Resolved 2026-08-14.** Registered criterion:
   - *Test:* TOST (two one-sided tests, α = .05 each) on the **paired within-participant difference** between the corrected and uncorrected ROI contrast estimates (the H2 confidence effect, extracted per participant from each of the three primary ROIs).
   - *Bounds, set by what the samples can support at ≥80% power:* **dz = ±0.65** for ±PetCO2 (N = 24; 82% power) and **dz = ±0.60** for ±RVT (N = 30; 84% power). Tighter bounds are not honestly supportable: dz = 0.3 would require N ≈ 90.
   - *Interpretation grid, fixed in advance:* TOST passes → correction leaves the effect **statistically equivalent**; TOST fails but the corrected effect remains significant in the ROI → **attenuated but present**; TOST fails and the corrected effect is not significant → **indeterminate** (not "abolished" — the design cannot distinguish abolition from power loss at these bounds).
   - *Context for the bound's meaning:* the paired SD is the SD of the difference between two estimates from nearly identical models. If the two estimates correlate at r = 0.95, dz = 0.65 on the difference corresponds to ≈ 0.21 between-participant SD of the estimate itself — a materially tight bound despite the modest N. The observed corrected–uncorrected correlation will be reported alongside the TOST so the bound can be read in raw units.
5. ~~**MAIA scoring check.**~~ **Resolved 2026-08-14:** MAIA-2 complete (all 37 items) for all 58 scored participants; 47 of 48 imaging participants covered. Scores in `Repo/Data/maia_scores.csv`. H6 finalised at N = 47.
