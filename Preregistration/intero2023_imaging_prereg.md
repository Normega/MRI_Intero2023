# Preregistration — Neural correlates of respiratory change awareness (MRI_Intero2023)

**Status:** draft for discussion — **pre-freeze**
**Version:** 0.2
**Date:** 2026-08-24 (v0.1: 2026-08-13)
**Analysis code:** `Repo/R/`, `Repo/Scripts/` (all preprocessing and behavioural code already written and versioned)

**Changes in v0.2, all made before the registration was frozen and before any model was fitted.** Each is evidenced in `VERIFICATION_REPORT.md` and recorded with its rationale in `Preregistration/methods_decisions.md`:

1. **§4** — the imaging inventory is corrected: sub-065, sub-076, sub-082 and sub-104 are reinstated, sub-051 and sub-056 are named as excluded, and the sample cascade and N table are rebuilt. Post-gate sample: **90 runs, 47 participants** (was 84 / 43).
2. **§7** — the single first-level model is split into three (**A** awareness, **B** confidence, **C** salience). Awareness and the confidence modulator are never estimated in the same design. This removes all 9 rank-deficient designs the single model produced.
3. **§3** — the awareness primacy is reversed: the categorical contrast (**H3**) is authoritative and the parametric confidence effect (**H2**) is supportive, with the cost of that choice registered in advance.
4. **§5** — C1's paced side averages over the paced conditions present in each run; trailing events are clipped to the scan end, counted, and never dropped.
5. **§9** — sensitivity analysis **S6** (two-run participants only) is added.

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

**Computed before the freeze, and used to write v0.2:** the design matrices themselves — rank, per-contrast SE multipliers, collinearity among task regressors, trial counts and estimability (`Scripts/intero2023_designqc.py`, `Scripts/intero2023_twomodel_qc.py`, outputs in `analysis/`). These are properties of the timing and the confound sheets. No BOLD voxel is read by either script, and no quantity below is a function of any contrast estimate. Run-level image quality (framewise displacement, tSNR) is used for exclusion decisions, as §4's registered motion gate already required.

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

Regions are specified a priori as ROIs in §8. Directional predictions are one-sided in the sense stated; whole-brain tests remain two-sided. Each hypothesis is assigned to one of the three first-level models defined in §7; **awareness and confidence are never estimated in the same model** (§7, and `Preregistration/methods_decisions.md` decision 2).

**H1 — Exogenous respiratory control (replication of IEAT Active Matching). Model A, N = 47.**
Two competing directions exist for the paced vs free contrast, and this document commits to the replication direction. The naive volitional-control account predicts paced > free activation in SMA, sensorimotor cortex, and insula. The IEAT found the opposite: exogenous pacing produced greater **deactivation** along the insula/operculum pathway, somatomotor cortex, ventral occipital cortex, and cerebellum relative to endogenous breathing.
*Registered prediction:* **paced < free** (relative deactivation) in the insula/operculum and somatomotor ROIs, replicating the IEAT.
*The test is diagnostic either way:* a paced > free effect in these regions would contradict the predecessor study in a larger sample and will be reported as a failed replication, not reframed. This hypothesis was revised from the volitional-control direction on 2026-08-14, on reading the IEAT results and **before any imaging analysis**.

**H3 — Awareness, categorical (AUTHORITATIVE). Model A, N = 23.**
Change trials answered with the correct direction ("aware") > change trials answered "No Change" ("unaware") in bilateral anterior insula and dACC. ROI-restricted inference only (three masks, Bonferroni ×3, §8).

*Why this is the authoritative awareness test (revised 2026-08-24, before any imaging analysis).* H3 is a within-participant comparison of two classes of trial. H2 is a parametric effect that is **confounded with that same comparison**: confidence and awareness covary in the behaviour, with a within-participant point-biserial r(confidence, aware) of median 0.316 (IQR 0.126–0.570, positive in 37 of 43 participants; aware trials rated a median 0.50 points higher). A measure that carries the confound cannot arbitrate it, so the categorical contrast is designated the confirmatory test and the parametric one is designated supportive. The earlier version of this document had the primacy the other way round, on power grounds alone.

*Registered limitation — the H3 sample is selected, and it is selected in a direction that matters.* Under the registered ≥3-trials-per-cell rule, H3 is estimable in 23 of the 47 participants. Those 23 have a pooled median of 4 unaware trials against 2 for the 24 who do not qualify (median unaware share of change trials 0.36 vs 0.19). The H3 sample therefore **over-represents participants who detected changes less often**, and its estimate does not generalise to the full sample. The subsample's detection rate is reported next to every H3 result. The 23 participants are fixed now, before any model is fitted, and listed in `analysis/design_qc_twomodel.csv`:

> sub-057, sub-060, sub-061, sub-064, sub-065, sub-066, sub-067, sub-068, sub-071, sub-073, sub-076, sub-077, sub-078, sub-084, sub-100, sub-102, sub-106, sub-108, sub-112, sub-120, sub-130, sub-132, sub-134

*Falsified if* no aware > unaware effect survives small-volume correction in either anterior insula ROI or the dACC.

**H2 — Awareness, parametric (SUPPORTIVE). Model B, N = 44.**
Trial-by-trial confidence on paced trials will correlate positively with activity in **bilateral anterior insula** and **dorsal anterior cingulate**.
*Rationale:* confidence is metacognitively valid here (type-2 AUC 0.661), and these regions are the canonical substrate for interoceptive salience and confidence.

*Registered statement of what H2 can and cannot show.* Model B contains no awareness regressors, so the confidence modulator is estimable in every run — but that is a property of the design, not of the constructs. **Model B's confidence effect absorbs the awareness effect rather than being adjusted for it**, and no claim of independence is made from it. A confidence effect in Model B is consistent with an awareness effect, a graded-confidence effect, or both, and Model B cannot distinguish them. Interpretation of H2 is conditioned on H3: where the two overlap spatially, the parsimonious reading is one effect, not two.

*Two inferences are ruled out in advance*, so the N asymmetry between the authoritative and supportive tests cannot be read as a result:

1. **A null H3 alongside a positive H2 is not evidence that confidence rather than awareness drives the effect.** At N = 23 against N = 44 that comparison is confounded with power.
2. **H3's estimate is not to be generalised beyond its subsample**, per the selection stated above, whichever way it comes out.

*Falsified if* no positive confidence effect survives small-volume correction in either anterior insula ROI.

**H4 — Salience acts on detection. Model C, N = 43.**
High-salience (step) > low-salience (ramp) change trials will elicit greater activity in the salience network (anterior insula, dACC).
*Rationale:* behaviourally, salience affected *noticing* (p = 0.043) but not *direction discrimination* (p = 0.26), so its neural effect is predicted to fall on detection-related rather than discrimination-related regions.
*Registered note (added 2026-08-24):* salience is **not crossed with awareness**. The median participant contributes 4 high-salience and 8 low-salience change trials, so a 2 × 2 split would leave cells of about one trial. Model C splits the change epochs by salience and by nothing else; the awareness × salience interaction is not tested, and is not exploratory either — it is not estimable.
*Falsified if* the effect appears instead in sensory/discrimination regions with no salience-network involvement.

**H5 — Physiological correction attenuates but does not abolish awareness effects.**
Adding physiological regressors (§6) will reduce cluster extent in H1–H4 but the anterior insula and dACC effects **in H3, and secondarily in H2,** will survive.
This is directional and consequential: PetCO2 is a **mediator** in this design, not a pure nuisance — the pacing manipulation changes end-tidal CO2 by construction. Regressing it out therefore removes part of the effect of interest for H1 in particular. Both corrected and uncorrected models will be reported side by side, and neither is designated "the" result for H1. Precedent: the IEAT found RVT associated with widespread deactivation, and that controlling respiration rate and RVT at both levels attenuated but did not eliminate its interoception deactivation effect.
*Note on the TOST target (revised 2026-08-24):* the equivalence bounds in §12.4 were set on the H2 confidence effect and remain defined on it, because the ±PetCO2 and ±RVT subsamples (N = 24 and 30) do not support the same test on H3's 23 participants, of whom fewer still have capnography. The TOST is therefore reported as an equivalence statement about the **supportive** measure, and is not read as an equivalence statement about H3.

**H6 — MAIA moderation of the ACC (direct replication). Model A, N = 46.**
Higher MAIA interoceptive sensibility will predict **greater (spared) ACC activity** in the paced-related contrasts, replicating the IEAT's central individual-differences result (peak −4, +12, +36, k = 1205).
*Test:* MAIA total (scored as in the IEAT: mean of the six subscales omitting Not-Worrying and Not-Distracting; scores in `Repo/Data/maia_scores.csv`, scoring code `intero2023_maia.R`) as a group-level covariate on the H1 contrast from Model A; small-volume inference in a 10 mm sphere at the IEAT peak (−4, +12, +36), plus the atlas dACC ROI as a secondary mask. **One-sided** in the replication direction, as the direction is fully specified by the prior study.
*Coverage and power:* 46 of the 47 analysable participants have complete MAIA (sub-065 has none). At N = 46, one-sided α = .05: 80% power for r ≥ 0.362; at the IEAT's observed total-score correlation (r ≈ 0.32), power ≈ 70%. Adequately powered for the replication-sized effect; this is stated so a marginal result is read against it.
*Rationale:* the strongest form of replication available — same construct, same region, a priori coordinate from an independent sample. MAIA subscale analyses (Attention Regulation vs Noticing dissociation reported in the IEAT) are exploratory (§10).
*Administration note (registered):* this study administered the 37-item MAIA-2 with response boxes coded 1–6; scores are converted to the canonical 0–5 metric before use. The IEAT used the original MAIA; the six retained subscales are common to both.

---

## 4. Sample, exclusions, and analysable N

Exclusions are already applied and are not contingent on imaging results.

**Corrected 2026-08-20, before the registration was frozen and before any model was fitted.** The imaging inventory used to build the original run list predated the 2026-08-06 re-conversion of the dataset and was wrong in two ways, both verified against the fMRIPrep 24.1.1 derivatives of record (`VERIFICATION_REPORT.md`, `Preregistration/section4_correction_draft.md`). The criteria applied below are the ones already registered, applied to a corrected inventory.

**Excluded a priori:**

- Pilot/QA identifiers P00, P01, P50, P999
- **sub-051** — behavioural session only; no MRI acquired.
- **sub-056** — `task-breathing_run-1` is a truncated acquisition (229 volumes; the event record runs 28.3 s past the end of the scan), and `task-breathing_run-2` fails the registered motion gate (mean FD 0.471 mm, 26.7 % of volumes over 0.5 mm). No usable breathing run remains.
- Runs with no usable event timing or no trigger-locked physiology

**Reinstated by the correction:**

- **sub-082 and sub-104.** The earlier version excluded them because "run labels are shifted in the derivatives tree". That is true of the Windows `Data/Ready` snapshot, which mapped scanner runs positionally; the BIDS tree of record maps by scanner slot and is correct (Run1 → `task-localizer_run-1`, 196 vols; Run3 → `task-breathing_run-1`, 242 vols; Run5 → `task-authenticity`, 172 vols, for both). `Data/qc_ledger.csv` already graded both breathing runs as exact. Each contributes one gate-passing breathing run (mean FD 0.241 and 0.154).
- **sub-065 and sub-076.** Both have complete two-run breathing sessions (242 volumes each) that the stale inventory omitted; the ledger still carried sub-065 `task-localizer_run-1` at 167 volumes, the truncated file that re-conversion replaced with the full 196, and that single row is the only volume-count disagreement between the ledger and the dataset of record across 222 checked runs. sub-076 had no imaging rows at all. Both pass the registered motion gate on both runs. **sub-065 has no MAIA questionnaire and therefore does not enter H6.**

**Analysable sample (established, not projected):**

| | runs | participants |
|---|---|---|
| Breathing runs with events | 105 | 54 |
| + usable imaging | **103** | **53** |
| − sub-056 (truncated run-1, gated run-2) | 101 | 52 |
| **− motion gate (11 of the remaining 101 runs)** | **90** | **47** |
| + RVT regressor | 64 | 30 (all runs) |
| + PetCO2 regressor | 51 | 24 (all runs) |
| + both physiological regressors | 43 | 20 (all runs) |
| Verified pacer adherence (≥1 run) | — | 24 |

**Additional motion exclusion, applied before modelling:** a run is dropped if mean framewise displacement > 0.5 mm or if > 25% of its volumes exceed FD 0.5 mm. A participant is dropped if this leaves fewer than one usable breathing run. This criterion is fixed now and will not be adjusted after seeing results.

**The gate has been applied** (2026-08-14, using fMRIPrep's FD column only — no BOLD signal or contrast has been examined). Across the full 103-run set it fails 12 runs, one of which is sub-056 run-2, already removed on the line above. It removes exactly the same 5 participants entirely as in the registered version (sub-058, -062, -080, -085, -094); two more lose one of two runs. Gated runs are unambiguous (median 81 spike volumes vs 3 in retained runs; median FD 0.52 vs 0.18 mm). **Post-gate: 90 runs, 47 participants** (`Repo/FirstLevel/firstlevel_manifest.csv` plus `firstlevel_manifest_candidates.csv`). None of the four reinstated participants has breathing-run RVT or PetCO2, so the physiological rows above and the H5 TOST bounds are unchanged.

*Registered check on a paradigm-specific selection risk:* deep paced breathing could inflate FD and the gate could then preferentially remove adherent breathers. It does not: fd–adherence correlation across belt-usable runs r = −0.003 (n = 81), and gated runs show *lower* pacer adherence (median 0.230 vs 0.335).

**Run count is not an inclusion criterion.** A participant enters the analysable sample when the registered criteria — usable event timing, trigger-locked imaging, and the motion gate — are met by **at least one** breathing run, irrespective of how many runs they contribute. Where a participant contributes one run, the fixed-effects step of §8 is that run. This is stated as a general rule because it cannot otherwise be applied consistently: two participants of the original registered 43 — **sub-074 and sub-114** — already contribute exactly one breathing run, so excluding sub-082 and sub-104 for contributing one run each would apply a standard to those two participants alone rather than a registered criterion. Sensitivity analysis **S6** (§9) tests every hypothesis on two-run participants only, which removes all four together.

**Estimability rule.** A contrast counts for a participant when every cell it involves has **≥ 3 trials pooled across that participant's gate-passing runs**, and the contrast is estimable in the corresponding design (§7). This is the rule used throughout; it is stricter than matrix rank, and it is what produces H3's N.

**Final registered Ns (post-gate, post-correction, under the §7 models):**

| Contrast | model | N |
|---|---|---|
| C1 paced vs free (H1) | A | **47** |
| C2 aware vs unaware (H3, authoritative) | A | **23** |
| C3 confidence parametric (H2, supportive) | B | **44** |
| C5 salience high vs low (H4) | C | **43** |
| H6 MAIA moderation | A | **46** |
| C4 wrong-direction vs aware | — | 13 — **not analysed as a planned contrast** |
| S1 ±PetCO2 (basis on all retained runs) | — | **22** — TOST bound dz = 0.65 (80%) |
| S2 ±RVT (regressor on all retained runs) | — | **26** — TOST bound dz = 0.59 (80%) |

Three participants — **sub-066, sub-070 and sub-114** — give the same confidence rating on every retained paced trial (all 14 rated 6, all 14 rated 7, and all 7 rated 7 respectively). A constant modulator is not a regressor, so they are **excluded from H2 only** and retained for every other hypothesis. This is a property of the behaviour, not of the sample correction, and the rule is fixed now, before any model is fitted.

H6 power updates to: 80% for r ≥ 0.362 one-sided (N = 46); ≈ 70% at the IEAT's r ≈ 0.32. Under S6 (N = 42) it is r ≥ 0.378 and ≈ 67%. **The correction does not change the study's power situation materially and is not presented as if it does** — the case for reinstating these participants is that the inventory used to exclude them was wrong, which is a correctness argument, not a statistical rescue.

---

## 5. Design and event model

Per breathing run, the following events are modelled (already extracted and validated against a hand-checked participant, 80 assertions in `Repo/Tests/`):

- **Paced breathing epochs** (~20 s), one per change trial, split by condition
- **Free breathing epochs** (~20 s), modelled as a separate trial type
- **Rating window**, modelled explicitly rather than left in the implicit baseline
- **Valence/sub-rating events** as nuisance

Trial counts per participant (median): 7 aware-correct, 2 unaware, 1 wrong-direction, 2 No Change, 14 free. The design is short — 14 paced trials per participant — and that shortness is the reason the awareness contrast is estimable in only 23 participants (§3, §4).

**Wrong-direction trials** (median 1 per participant) are modelled as a **separate regressor of no interest** so they do not contaminate the aware or unaware regressors. They are not contrasted.

**No Change trials** (median 2) are modelled as their own regressor and are not treated as correct rejections.

**Trailing events are clipped to the scan end, counted, and never dropped** (registered 2026-08-24). Any modelled event extending past `n_vol × TR` has its duration truncated to the end of the acquisition; an event whose truncated duration is ≤ 0 — one beginning after the scanner stopped — is dropped, having no data to explain. These are runs in which the scanner stopped during a final rating, not runs with bad events: dropping the whole event would discard seconds of acquired signal and silently change the design, while extending the design matrix past the data is not possible. The count is reported per run. Across the 90 modelled runs, three carry a clipped event:

| run | volumes | events clipped | event over-run |
|---|---|---|---|
| sub-052 `task-breathing_run-1` | 244 | 1 | 11.0 s |
| sub-053 `task-breathing_run-1` | 250 | 2 | 25.6 s |
| sub-053 `task-breathing_run-2` | 260 | 1 | 5.6 s |

The larger over-runs elsewhere in the event files (sub-058, sub-094: 19–67 s) belong to runs the motion gate already removed, so they never reach a model.

**The paced side of C1 averages over the paced conditions present in the run** (registered 2026-08-24). C1 contrasts `free` against the unweighted mean of whichever of `paced_aware`, `paced_unaware` and `paced_nochange` have trials in that run. The earlier fixed 1/3 weighting is undefined in 18 of 91 runs — 17 of them missing `paced_unaware`, because the participant never answered "No Change" on a change trial — and the hypothesis is a claim about paced versus free breathing, not about the 1/3 weights. Precision is unaffected: median SE multiplier `sqrt(c'(X'X)⁻¹c)` is 0.0827 for the adaptive form and 0.0830 for the fixed form in runs where both are defined; the change moves H1 from N = 45 to N = 47, and the two recovered participants are recovered because a cell was empty, not because their data improved. *Stated as a limitation:* the paced side is a slightly different weighted average across participants. Since the weights carry no theoretical content, this is a smaller distortion than dropping 18 runs. The per-run composition is logged in `analysis/design_qc.csv`. `paced_wrongdir` and `paced_noresponse` remain regressors of no interest and stay out of the contrast.

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

## 7. First-level models

**Three models per run**, all in MNI152NLin2009cAsym space (fMRIPrep 24.1.1 derivatives, already generated), all sharing the same nuisance set, the same HRF (canonical + temporal derivative), the same boxcar-over-epoch-duration convention, and `drift_model = None` (the cosine columns in the confound sheet *are* the high-pass). They differ only in how the paced epochs are parcelled.

**Awareness and the confidence modulator are never in the same design** (registered 2026-08-24; rationale and evidence in `Preregistration/methods_decisions.md` decision 2). In the single model this document previously specified, the confidence modulator lives on the same epochs the awareness regressors partition, and the two are collinear by construction: median absolute r between the modulator and the aware-minus-unaware direction is **0.437**, exceeds 0.5 in **26 of 70** runs, reaches **0.942** at worst, and makes **9 of 91 designs rank-deficient** — every one of them recovering full rank when the modulator alone is removed. Under the split, rank deficiency is **0 of 90 in every model**. This is the removal of a defect that made nine designs unfittable as written, not a tuning improvement.

**Model A — the awareness model (authoritative). Yields H1, H3, H6.**

1. Paced breathing — aware (correct direction)
2. Paced breathing — unaware (responded "No Change")
3. Paced breathing — wrong direction *(no interest)*
4. Paced breathing — No Change trials
5. Free breathing
6. Rating window
7. Valence/sub-rating events *(no interest)*

No parametric modulator.

**Model B — the confidence model (supportive). Yields H2 and nothing else.**
The paced epochs collapse to a single **`trial`** regressor carrying the mean-centred **confidence** modulator; `free`, `rating` and `subrating` are unchanged from Model A. `trial` spans **all** paced epochs, including No Change trials — that is, exactly the epochs the confidence modulator was originally defined on, so the split changes the design without changing what confidence modulates. Restricting `trial` to change trials would additionally cost one participant (sub-112, whose confidence varies only once No Change trials are counted), giving N = 43 rather than 44.

Model B contains no awareness regressors, so nothing remains for the modulator to alias against. **This separates the regressors, not the constructs** — see H2 in §3 for what may and may not be concluded from it.

**Model C — the salience model. Yields H4.**
The paced **change** epochs are split into high-salience (step) and low-salience (ramp) rather than by awareness; `paced_nochange`, `free`, `rating` and `subrating` are unchanged from Model A. No parametric modulator. Model C is added because §7 as previously written contained no salience-split regressor, so H4 had nothing to be tested on. Salience is not crossed with awareness (§3, H4).

**Parametric modulators.** Confidence (1–7) enters Model B only, mean-centred over the trials actually modelled in that run. **Arousal** (1–7) is not a registered test given its flat behavioural profile and is reported as exploratory (§10); where reported it is fitted in a separate model of the same form as Model B, never added to A or C.

**Nuisance regressors (base model, all participants, identical across A/B/C):**

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

**Design QC, computed before the freeze and before any model was fitted.** `Scripts/intero2023_designqc.py` and `Scripts/intero2023_twomodel_qc.py` build every design matrix and report rank, per-contrast SE multipliers, VIF among the task regressors, and the modulator/awareness correlation; outputs are committed as `analysis/design_qc.csv` and `analysis/design_qc_twomodel.csv`. No BOLD image is read by either.

---

## 8. Group model, ROIs, and inference

**Group model:** one-sample t-test on first-level contrast images, participant as the unit. Where a participant contributes two runs, run-level contrasts are combined at a fixed-effects intermediate level before entering the group model.

**Why fixed-effects combination (decided 2026-08-14, before any imaging analysis):**

- *Run-as-replication is not available for the awareness contrast.* Within a single run the cells collapse (median 3 aware / 1 unaware per run); requiring ≥3 trials per cell in both runs leaves exactly 1 participant, versus the 23 that pooling yields (§4). Pooling is the condition for H3's estimability, not a trade-off against replication.
- *The runs are exchangeable.* Direction accuracy 0.800 vs 0.754 (p = 0.22), confidence 5.51 vs 5.31, P(say change) 0.726 vs 0.713 — a slight non-significant run-2 softening, no heterogeneity that would make an average misleading. Each run keeps its own GLM and baseline, so only the contrast is averaged.
- *Pooling materially improves per-participant reliability.* Cross-run individual differences are stable (direction accuracy r = 0.688, confidence r = 0.736, n = 42/51); Spearman-Brown gives the pooled estimate reliability ≈ 0.82 vs ≈ 0.69 for one run.
- *Run as a random factor is degenerate with two levels* and is not used.

**Registered reliability diagnostic (descriptive, not a gate):** for C1 (paced vs free, Model A) and C3 (confidence parametric, Model B), run-1 vs run-2 agreement of first-level contrast estimates — spatial correlation of contrast maps within participants, and cross-participant ICC of ROI means. Reported alongside the behavioural benchmark (r ≈ 0.7). H3 is excluded from this diagnostic because it is estimable in both runs for only 1 participant; this is stated now so its absence is not read as an omission later.

**A priori ROIs** (Harvard–Oxford probabilistic atlas, thresholded at 50%, bilateral where applicable):

1. Anterior insula (primary for H3/H2, with dACC)
2. Posterior insula
3. Dorsal anterior cingulate / mid-cingulate (primary for H3/H2, with anterior insula)
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

| N | hypothesis | α = 0.05 (ROI) | α = 0.001 (cluster-forming) |
|---|---|---|---|
| 47 | H1 (Model A) | 0.42 | 0.64 |
| 44 | H2 (Model B) | 0.43 | 0.67 |
| 43 | H4 (Model C) | 0.44 | 0.67 |
| 24 | ±PetCO2 (S1) | 0.60 | 0.96 |
| 23 | **H3 (Model A)** | **0.61** | **0.98** |

**H3 at N = 23 requires a large effect (d ≈ 1.0) at whole-brain thresholds, which is why it is ROI-restricted.** It is nonetheless the *authoritative* awareness test (§3): being the better-identified contrast and being the better-powered one are different properties, and this document deliberately assigns primacy on identification rather than on power. The consequence is stated in advance in §3 and repeated here — **a null H3 is not to be reinterpreted afterwards, and a null H3 beside a positive H2 is not evidence that confidence rather than awareness drives the effect**, because that comparison (N = 23 vs 44) is confounded with power.

---

## 9. Registered sensitivity analyses

Each is fully specified now and will be reported whatever the outcome. Where a sensitivity analysis names a hypothesis, it is refitted in **that hypothesis's own model** (§7) — H1/H3/H6 in Model A, H2 in Model B, H4 in Model C — with the physiological or subsetting change applied identically to all of them. S6 is ordered last because it was added in v0.2; the numbering of S1–S5 is unchanged so that references to them in the v0.1 record still resolve.

**S1 — ±PetCO2 (the primary methodological question, N = 24).**
Fit H1–H4 with and without the PetCO2 basis in the 24 participants who have it on all runs (H3 in whichever of them are among its 23; that overlap is reported rather than assumed). Compare (a) peak z and cluster extent in each ROI, (b) whether each hypothesis's conclusion changes. Because this is a within-participant model comparison, it is better powered than the between-group comparison originally considered.

**S2 — ±RVT (N = 30).** As S1, with the RVT regressor.

**S3 — CO2-available vs CO2-unavailable participants.** Test whether the 24 participants with usable capnography differ from the remainder on the uncorrected H2 contrast (H2 rather than H3, because the subgroup comparison needs the larger sample to be worth running at all). **Registered caveat:** CO2 availability is near-collinear with cohort and protocol version (the capnograph outage was an unbroken date block from 2023-10-12 to 2024-01-18), so this comparison **cannot isolate CO2** and will be reported as a cohort check only, not as evidence about physiological correction.

**S4 — Adherence conditioning (N = 24).** Repeat **H3 and H2** restricted to participants with verified pacer adherence. Adherence is a run-level property (per-epoch classification is underpowered at AUC 0.66–0.69), and this is registered as robustness, not as the primary analysis.

**S6 — Two-run participants only (N = 43).**
Every hypothesis is additionally tested on the subset contributing **two** breathing runs. The subset is defined by run count alone, so it removes **sub-074, sub-114, sub-082 and sub-104** together — two participants from the original registered sample and two of the reinstated ones — and cannot be used to remove the reinstated participants selectively.

| | participants | H1 | H3 | H2 | H4 | H6 |
|---|---|---|---|---|---|---|
| **Primary** (corrected sample) | 47 | 47 | 23 | 44 | 43 | 46 |
| **S6** (two-run participants only) | 43 | 43 | 23 | 41 | 43 | 42 |

*Interpretation fixed in advance:* agreement between the primary and S6 results is reported as such; divergence is reported as a limitation of the single-run participants' precision and **not** used to select which analysis to present.

**S5 — Protocol version.** v1 (n = 221 trials) and v2 (n = 1218 trials) differ in rating structure. Include protocol version as a group-level covariate; report whether the **H3 and H2** conclusions change.

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
5. ~~**MAIA scoring check.**~~ **Resolved 2026-08-14, updated 2026-08-24:** MAIA-2 complete (all 37 items) for all 58 scored participants. Under the corrected sample, 46 of the 47 analysable participants are covered — sub-065 has none. Scores in `Repo/Data/maia_scores.csv`. **H6 finalised at N = 46.**
6. ~~**Imaging inventory.**~~ **Resolved 2026-08-20:** sub-065, sub-076, sub-082 and sub-104 reinstated; sub-051 and sub-056 named as excluded. Rationale in §4 and `Preregistration/section4_correction_draft.md`; evidence in `VERIFICATION_REPORT.md`.
7. ~~**One model or two.**~~ **Resolved 2026-08-24:** three models (§7). Awareness and confidence are never in the same design. Rationale in `Preregistration/methods_decisions.md` decision 2.
8. ~~**Which awareness test is authoritative.**~~ **Resolved 2026-08-24:** the categorical contrast H3, with H2 supportive, and with the power cost and the H3 subsample's selection both registered in advance (§3, §8).
9. ~~**C1's paced side, and trailing events.**~~ **Resolved 2026-08-24:** average over the paced conditions present in the run; clip trailing events to the scan end and log the count (§5).
