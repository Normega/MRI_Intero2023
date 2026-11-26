# Preregistration — Neural correlates of respiratory change awareness (MRI_Intero2023)

**Status:** final for registration (submitted to OSF 2026-11-26). Any change after the freeze is logged in `Preregistration/deviations.md`.
**Version:** 0.8
**Date:** 2026-11-25 (v0.7: 2026-11-25; v0.6: 2026-11-25; v0.5: 2026-11-25; v0.4: 2026-11-25; v0.3: 2026-11-24; v0.2: 2026-08-24; v0.1: 2026-08-13)
**Analysis code:** `Repo/R/`, `Repo/Scripts/` (all preprocessing and behavioural code already written and versioned)

**Changes in v0.2, all made before the registration was frozen and before any model was fitted.** Each is evidenced in `VERIFICATION_REPORT.md` and recorded with its rationale in `Preregistration/methods_decisions.md`:

1. **§4** — the imaging inventory is corrected: sub-065, sub-076, sub-082 and sub-104 are reinstated, sub-051 and sub-056 are named as excluded, and the sample cascade and N table are rebuilt. Post-gate sample: **90 runs, 47 participants** (was 84 / 43).
2. **§7** — the single first-level model is split into three (**A** awareness, **B** confidence, **C** salience). Awareness and the confidence modulator are never estimated in the same design. This removes all 9 rank-deficient designs the single model produced.
3. **§3** — the awareness primacy is reversed: the categorical contrast (**H3**) is authoritative and the parametric confidence effect (**H2**) is supportive, with the cost of that choice registered in advance.
4. **§5** — C1's paced side averages over the paced conditions present in each run; trailing events are clipped to the scan end, counted, and never dropped.
5. **§9** — sensitivity analysis **S6** (two-run participants only) is added.

**Changes in v0.3 (2026-11-24), all made before the registration was frozen and before any model was fitted.** They come from an independent re-check of every registered number against the event sheets, manifests and behaviour files, and a re-reading of the IEAT paper. No BOLD image was read. Rationale for each is in `Preregistration/methods_decisions.md` decisions 5–9.

6. **§3, §8: H6 is re-framed as a two-sided conceptual extension.** The IEAT's MAIA effect was a covariate on Active Interoception vs Active Exteroception. This task has no exteroceptive condition, and the paced − free contrast H6 uses is the analogue of the IEAT's Matching vs Interoception contrast, which the IEAT never related to MAIA. The prior study therefore does not fix the sign, and the one-sided test is withdrawn.
7. **§5, §7: the design has no implicit baseline, and this is now registered.** Modelled events cover a median 94% of each run (minimum 63%). Only between-condition contrasts are interpreted. "Deactivation" means relative to free breathing, and no condition-vs-rest contrast will be reported.
8. **§4, §8, §9, §12.4: the physiological-subset Ns are corrected to their post-gate values.** These are PetCO2 22, RVT 26, both 18, and verified adherence 26 (superseded by v0.4 item 12). The several places that carried pre-gate values (24 / 30 / 20 / 24) are fixed. The TOST bounds are restated at the true Ns: dz = ±0.65 and ±0.60 at 80% power.
9. **§3: H5's claim about H3 is descriptive.** Only 12 of H3's 23 participants had PetCO2 and 14 had RVT (13 and 15 after v0.4). The confirmatory H5 test is the TOST on H2.
10. **§4: run-level exceptions for Model B are registered.** These are runs whose confidence is constant within the run.
11. **Minor corrections (v0.3).** H2's falsification criterion includes the dACC. H4 gains an explicit inference set and a null outcome. The §5 trial medians are updated to the corrected sample. The direction-accuracy definition is stated. The IEAT's insula laterality is corrected (bilateral in its Results).

**Changes in v0.4 (2026-11-25), pre-freeze, before any model was fitted.** These follow a re-examination of the physiological recordings (`Analysis/pmu_rvt_validation.md`; code changes in commits `4aae7cd` and `bd88f7f`). No BOLD image was read. Rationale is in `methods_decisions.md` decisions 10–11.

12. **§4, §6, §9: physiological coverage recovered by fixing three pipeline bugs, with the registered pipeline and thresholds otherwise unchanged.** 7 RVT and 5 PetCO2 breathing runs are recovered:
    - The physio cache-key regex could not parse sub-070's part label, so sub-070's breathing runs had been built from the wrong recording.
    - A stale hard-coded exclusion had blocked sub-082 and sub-104.
    - sub-065 and sub-076 had been skipped because the old imaging snapshot did not contain them.

    Biopac RVT on every retained run: 26 → **30** participants. PetCO2: 22 → **25**. Both: 18 → **21**.
13. **§6, §9: the Siemens PMU respiration cushion is registered as a fallback RVT source, used in a new robustness analysis S2b only.** The confirmatory S2/H5 test stays on the Biopac belt. PMU RVT couples to the global signal at least as well as the Biopac belt, but the two belts' RVT time courses agree poorly. The saturation gate was revised after some fallback-run quality metrics had been seen. That sequence is disclosed in §6.
14. **§3, §4, §12.4: the TOST N is the H2 subset, and the bounds are restated at it.** v0.2/v0.3 quoted the whole-sample physio N (22 / 26) for a test defined on the H2 effect; the correct H2 ∩ RVT N was 24. At v0.4: PetCO2 N = 25, dz = ±0.61; Biopac RVT N = 27, dz = ±0.59 (80% power).

**Changes in v0.5 (2026-11-25), pre-freeze, before any model was fitted.** No BOLD intensity was read. Rationale is in `methods_decisions.md` decisions 12–13.

15. **§3, §8: the ROI masks are now constructible, and H1 has a fixed inference set.** v0.4 named Harvard–Oxford masks the atlas does not contain: it has a single insula label (no anterior/posterior division) and no dACC/mid-cingulate label (its anterior cingulate division runs from pregenual to mid-cingulate cortex). H1's masks and its Bonferroni count were not stated. The masks are now defined from the CAREN consensus network atlas (Doucet, Lee & Frangou, 2019) by network × AAL-3 region membership, following the procedure used in the lab's MindLock study. The salience network supplies the anterior insula and dACC, and the somatomotor network supplies the posterior insula, central operculum and sensorimotor cortex. H1 is tested in three masks, Bonferroni ×3. A Harvard–Oxford version is registered as an atlas-sensitivity check. The masks were built and checksummed before the freeze (`ROIs/roi_manifest.csv`).
16. **§3, §7, §8, §10: H7 is added, a confirmatory test of between-network coupling.** H7a is paced vs free (N = 47) and H7b is aware vs unaware (N = 23), each for salience–somatomotor and salience–executive coupling. The method follows the lab's MindLock connectivity analysis, with two changes for this design. Task-evoked responses are removed before correlating, because the epochs here are short. And a global-signal interpretation gate is fixed in advance, because paced breathing moves CO2 and H1 predicts lower regional activity. Feasibility was counted from event timing only. Rationale is in `methods_decisions.md` decision 13.

**Changes in v0.6 (2026-11-25), pre-freeze, before any model was fitted.** Decided on trial counts and confound sheets only; no BOLD intensity was read. Rationale is in `methods_decisions.md` decisions 14–16.

17. **§3, §4, §8, §9: H3 now includes every participant in whom the contrast exists, with a precision-weighted group model.** The v0.5 rule (≥ 3 trials per cell, pooled) excluded 24 of 47 participants, 23 of them on the unaware cell alone. Participants answer "No Change" on only 28% of change trials, and the design offers at most 12. The rule also selected toward less frequent detectors. H3 now enters every participant with at least one aware and one unaware trial in the same gate-passing run: **N = 45**. sub-070 and sub-114 never answered "No Change" on a change trial and cannot enter. The participant's contrast is the fixed-effects combination of the runs that contain both cells. The group model for H3 is FSL FLAME 1 on each participant's effect and variance, so a participant with few trials gets little weight instead of being dropped. The v0.5 rule and its 23 participants are kept as sensitivity analysis **S7** and as the **H7b** sample, which still needs ≥ 3 epochs per cell to estimate a correlation.
18. **§7: the CSF and white-matter regressors are regenerated from the fMRIPrep derivatives of record.** In the 95 sheets of the registered participants they came from an earlier anatomical pass (median r with the record 0.985 and 0.996; `VERIFICATION_REPORT.md` §2a). Every other column already reproduced bit-exactly. The nuisance set is unchanged, and no gate or N depends on these columns (`Scripts/intero2023_tissue_confounds.R`, `Analysis/tissue_confounds_regen.csv`).
19. **§10: an exploratory graded-awareness model is added**, in which a low-confidence aware trial counts as less aware than a high-confidence one.

**Changes in v0.7 (2026-11-25), pre-freeze, before any model was fitted.** Decided from design matrices and trial counts only; no BOLD intensity was read. Rationale is in `methods_decisions.md` decisions 17–18.

20. **§3, §4, §8: H3's first level is one joint model of both runs.** Under v0.6, a run with aware trials but no unaware trial could not form the contrast, so 56 of 290 aware trials (19%, in 14 participants) were discarded. The joint model shares the task regressors across runs and keeps the nuisance regressors and intercept run-specific. Shared free, rating and sub-rating conditions link the runs, so every trial enters. By design algebra this raises H3's total precision by 12.5%, and cuts the contrast variance by a median 22% in the 14 affected participants (`Scripts/intero2023_h3_joint_qc.py`, `Analysis/h3_joint_design_qc.csv`). The entry rule is unchanged in effect: the same 45 participants. The v0.6 per-run form is kept as a robustness check.
21. **§5 and the design QC: contrast standard errors corrected.** Both design-QC scripts scaled the contrast by the column SDs instead of dividing by them. They also built designs with the physiological columns, which belong to the sensitivity models only. Both are fixed and the QC re-run. Every rank verdict, estimability verdict and N is unchanged. The SE values quoted since v0.2 were wrong by about a factor of 3. The C1 statement in §5 is corrected: both forms have median SE 0.235. The registered comparison it supports, that the adaptive and fixed forms are equally precise, still holds.

**Changes in v0.8 (2026-11-25), pre-freeze, before any model was fitted.** Rationale is in `methods_decisions.md` decision 19.

22. **§3 H5, §4, §9, §12.4: the confirmatory H5 equivalence test moves to H3.** H3 is the authoritative awareness test, so the registered answer to "is the awareness effect respiratory or CO2-driven?" has to be about H3. That confound is also sharpest for H3: aware and unaware trials may differ in how closely breathing followed the pacer. At v0.5's N = 23 the H3 subsets (13 and 15) could not carry the test. At N = 45 they are 25 and 28, as large as H2's. The H3 TOST is confirmatory **when H3 is supported**, with bounds dz = ±0.61 (PetCO2, N = 25) and ±0.58 (Biopac RVT, N = 28). If they disagree, H3 governs and H2's TOST is supportive. If H3 is not supported, H5 is read from the H2 TOST and stated as a conclusion about the supportive measure only. The design-QC outputs are now written inside the repository (`Analysis/`); the copies committed before were from before the v0.7 SE fix.

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

**Computed before the freeze, and used to write v0.2 and v0.3:** the design matrices themselves. These are rank, per-contrast SE multipliers, collinearity among task regressors, event coverage of the scan, trial counts, and estimability (`Scripts/intero2023_designqc.py`, `Scripts/intero2023_twomodel_qc.py`, outputs in `analysis/`). These are properties of the timing and the confound sheets. No BOLD voxel is read by either script, and no quantity below is a function of any contrast estimate. Run-level image quality (framewise displacement, tSNR) is used for exclusion decisions, as §4's registered motion gate already required.

The ROI masks (§8) were built before the freeze from atlas files and fMRIPrep's run-level `desc-brain_mask` images. The brain masks fix the target grid and the group coverage mask. They are automated preprocessing outputs, and no BOLD intensity is read to build them.

The global signal *has* been used, for one narrow purpose: validating that the physiological regressors have the physiologically correct sign (§6). No task regressor has been related to any imaging measure.

Because the behavioural results are known, hypotheses H2–H4 below are **informed by observed behaviour**. They are predictions about *brain* data, which remain untested, but they are not blind to the behavioural pattern. Anything genuinely exploratory is marked as such in §10 and will be reported as exploratory regardless of outcome.

---

## 1. Study aims

The task manipulates the rate of paced breathing (faster / slower / no change, delivered as either an abrupt step or a gradual ramp) and asks participants to report what their breathing did. It also includes unpaced free-breathing periods. The aims:

1. Identify the neural response to **voluntary paced breathing** relative to free breathing.
2. Identify the correlates of **awareness of a respiratory change** — the central question.
3. Determine whether those correlates survive **correction for physiological (vascular) confounds**, which in this paradigm are partly manipulated by the task itself.
4. **Replicate key results of the direct predecessor study**, the IEAT (Farb, Zuo & Price, 2023, eNeuro 10(6):ENEURO.0088-23.2023), in a larger sample with an added awareness manipulation.
5. Test whether the salience network couples more strongly with the somatomotor and executive networks when breathing is paced, and when a change is noticed (H7, added v0.5).

### 1b. Relation to the predecessor study (IEAT)

The IEAT (N = 22, 44 scans) compared attention to the breath against a difficulty-matched visual task, and included an **Active Matching** condition — breathing paced to an external stimulus — that is the direct ancestor of this study's paced epochs. Its results constrain several predictions here:

- **Interoceptive attention produced widespread cortical deactivation**, not activation of interoceptive cortex, with preserved performance.
- **Exogenous pacing (Active Matching vs Active Interoception) drove additional deactivation *bilaterally* along the insula/operculum pathway, somatomotor cortex, ventral occipital cortex, and cerebellum** (IEAT Results, Fig. 5, Table 10). The IEAT abstract singles out the right insula as implicated only by this contrast, not by interoception vs exteroception; its Results describe the insula/operculum effect as bilateral, and H1 is specified bilaterally accordingly.
- **Self-reported interoceptive sensibility (MAIA) predicted sparing from deactivation in the anterior cingulate** and left-lateralised language regions, **on the Active Interoception vs Active Exteroception contrast** ("greater MAIA scores were associated with reduced deactivation during Active Interoception relative to Active Exteroception"). The ACC cluster retained as the IEAT's PPI seed is reported in its Results at k = 1205, (−4, +12, +36). Its Methods give a different description of the dACC cluster (k = 794, (−4, +58, +14)), and the paper does not reconcile the two. This document uses the Results coordinate, because it is the cluster the IEAT carried forward. ROI-level MAIA correlations were |r| ≈ 0.27–0.32 (IEAT Table 9). ACC–DAN connectivity increased during interoception.
- Its discussion draws a distinction directly relevant here: **anterior insula tracks liminal detection accuracy** (heartbeat-style paradigms), while **attention to salient respiratory signals implicates the ACC**. The present task's respiratory changes are salient and supraliminal (direction accuracy 0.78), which strengthens the a priori case for the ACC alongside the insula in the awareness hypotheses.

Correspondence between designs: this study's **paced epochs ≈ Active Matching** (exogenous respiratory control); **free epochs ≈ Active Interoception** (endogenous rhythm under interoceptive monitoring). The present design adds what the IEAT could not ask: trial-by-trial *awareness* of an imposed respiratory change.

MAIA-2 items are present in this study's screener battery, so the sensibility moderation can be carried over (H6). It can be carried over only as a **conceptual extension**, though, because this task has no exteroceptive condition to reproduce the IEAT's moderated contrast.

---

## 2. Key behavioural findings that constrain the imaging hypotheses

These are established, not predicted:

| Finding | Value |
|---|---|
| Three-category identification accuracy (chance = 1/3) | **0.525**, κ = 0.271, χ²(4) = 141.7, p = 1e-29 |
| Ordinal SDT sensitivity | **0.698** probit units per stimulus step, z = 13.2 |
| Direction accuracy on true-change trials, among trials where a direction was reported (i.e. excluding "No Change" responses) | **0.778**; per-participant median 0.833, 85% above chance. Over *all* true-change trials, including "No Change" responses, the correct-direction rate is 0.562 |
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

Regions are specified a priori as ROIs in §8. Directional predictions are one-sided in the sense stated; whole-brain tests remain two-sided. H1–H6 are each assigned to one of the three first-level models defined in §7, and H7 uses the connectivity model at the end of §7; **awareness and confidence are never estimated in the same model** (§7, and `Preregistration/methods_decisions.md` decision 2).

**H1 — Exogenous respiratory control (replication of IEAT Active Matching). Model A, N = 47.**
Two competing directions exist for the paced vs free contrast, and this document commits to the replication direction. The naive volitional-control account predicts paced > free activation in SMA, sensorimotor cortex, and insula. The IEAT found the opposite: exogenous pacing produced greater **deactivation** along the insula/operculum pathway, somatomotor cortex, ventral occipital cortex, and cerebellum relative to endogenous breathing.
*Registered prediction:* **paced < free** (relative deactivation) in the insula/operculum and somatomotor ROIs, replicating the IEAT.
*Inference set (registered v0.5):* three bilateral masks, all drawn from the somatomotor network (§8): **posterior insula, central operculum and sensorimotor cortex**. Small-volume FWE, one-sided paced < free, Bonferroni ×3, plus the two-sided whole-brain test. *Supported if* paced < free survives in at least one of the three masks. *Falsified if* it survives in none. The opposite tail (paced > free) is also tested in the same three masks at the same correction, and a surviving effect there is the failed replication described next. "Deactivation" here means lower activity during paced than during free epochs. The design has no rest baseline (§5), so nothing is claimed about either condition relative to rest.
*The test is diagnostic either way:* a paced > free effect in these regions would contradict the predecessor study in a larger sample and will be reported as a failed replication, not reframed. This hypothesis was revised from the volitional-control direction on 2026-08-14, on reading the IEAT results and **before any imaging analysis**.

**H3 — Awareness, categorical (AUTHORITATIVE). Model A, N = 45 (v0.6; was 23).**
Change trials answered with the correct direction ("aware") > change trials answered "No Change" ("unaware") in bilateral anterior insula and dACC. ROI-restricted inference only (three masks, Bonferroni ×3, §8). The first level is one joint model of the participant's gate-passing runs (v0.7, §8). The group model is precision-weighted (FSL FLAME 1, §8), because participants contribute between 1 and 8 unaware trials.

*Why this is the authoritative awareness test (revised 2026-08-24, before any imaging analysis).* H3 is a within-participant comparison of two classes of trial. H2 is a parametric effect that is **confounded with that same comparison**: confidence and awareness covary in the behaviour, with a within-participant point-biserial r(confidence, aware) of median 0.316 (IQR 0.126–0.570, positive in 37 of 43 participants; aware trials rated a median 0.50 points higher). A measure that carries the confound cannot arbitrate it, so the categorical contrast is designated the confirmatory test and the parametric one is designated supportive. The earlier version of this document had the primacy the other way round, on power grounds alone.

*Who enters, and why the weighting (registered v0.6).* A participant enters H3 with at least one aware and one unaware trial across their gate-passing runs. That is 45 of the 47 participants. The two who cannot enter, **sub-070 and sub-114**, never answered "No Change" on a change trial. Within those 45, unaware counts run from 1 to 8 (median 3). Every trial enters, because the first level is a joint model of both runs (§8). A per-run model would discard the aware trials in runs with no unaware trial: 56 of 290 (19%), in 14 participants. A participant with one unaware trial has a far noisier contrast than one with eight. The FLAME group model estimates each participant's within-participant variance and the between-participant variance and weights each participant by the inverse of their sum. Weak participants are thus down-weighted rather than excluded.

*What this changes about the sample.* The v0.5 rule (≥ 3 trials per cell, pooled) kept 23 participants. Those 23 said "No Change" more often (median 4 unaware trials against 2; unaware share of their aware + unaware trials 0.36 against 0.19 for the 24 excluded), so the v0.5 sample over-represented less frequent detectors. The v0.6 sample excludes only two participants, both by necessity. The residual selection is that the most frequent detectors contribute the least precise contrasts and so carry the least weight. The weight distribution is reported with the result. The v0.5 set is retained as sensitivity analysis S7 (§9), with its list fixed here:

> sub-057, sub-060, sub-061, sub-064, sub-065, sub-066, sub-067, sub-068, sub-071, sub-073, sub-076, sub-077, sub-078, sub-084, sub-100, sub-102, sub-106, sub-108, sub-112, sub-120, sub-130, sub-132, sub-134

Both sets, with per-participant trial counts and the runs that contribute, are frozen in `Data/h3_sets.csv`.

*Falsified if* no aware > unaware effect survives small-volume correction in either anterior insula ROI or the dACC.

**H2 — Awareness, parametric (SUPPORTIVE). Model B, N = 44.**
Trial-by-trial confidence on paced trials will correlate positively with activity in **bilateral anterior insula** and **dorsal anterior cingulate**.
*Rationale:* confidence is metacognitively valid here (type-2 AUC 0.661), and these regions are the canonical substrate for interoceptive salience and confidence.

*Registered statement of what H2 can and cannot show.* Model B contains no awareness regressors, so the confidence modulator is estimable in every run — but that is a property of the design, not of the constructs. **Model B's confidence effect absorbs the awareness effect rather than being adjusted for it**, and no claim of independence is made from it. A confidence effect in Model B is consistent with an awareness effect, a graded-confidence effect, or both, and Model B cannot distinguish them. Interpretation of H2 is conditioned on H3: where the two overlap spatially, the parsimonious reading is one effect, not two.

*Two inferences are ruled out in advance*, so the N asymmetry between the authoritative and supportive tests cannot be read as a result:

1. **A null H3 alongside a positive H2 is not evidence that confidence rather than awareness drives the effect.** The Ns are now similar (45 and 44), but they are not comparable in precision. H3 contrasts two cells of about 3 and 4 trials, whereas H2's modulator uses every paced trial. The comparison therefore remains confounded with power.
2. **The S7 estimate (the 23-participant v0.5 subsample) is not to be generalised beyond that subsample**, per the selection stated above, whichever way it comes out.

*Falsified if* no positive confidence effect survives small-volume correction in any of the three masks: either anterior insula ROI or the dACC (Bonferroni ×3, the same inference set as H3). *(v0.3: the dACC was previously omitted from this criterion although it is tested.)*

**H4 — Salience acts on detection. Model C, N = 43.**
High-salience (step) > low-salience (ramp) change trials will elicit greater activity in the salience network (anterior insula, dACC).
*Rationale:* behaviourally, salience affected *noticing* (p = 0.043) but not *direction discrimination* (p = 0.26), so its neural effect is predicted to fall on detection-related rather than discrimination-related regions.
*Registered note (added 2026-08-24):* salience is **not crossed with awareness**. The median participant contributes 4 high-salience and 8 low-salience change trials, so a 2 × 2 split would leave cells of about one trial. Model C splits the change epochs by salience and by nothing else; the awareness × salience interaction is not tested, and is not exploratory either — it is not estimable.
*Inference set:* the same three masks as H3 (left and right anterior insula, dACC; small-volume FWE, Bonferroni ×3), plus the whole-brain test.
*Supported if* high > low survives small-volume correction in at least one of the three masks. *Falsified if* it survives in none. An effect found only in sensory or discrimination regions at whole-brain level, with no mask surviving, is also falsifying: it is reported as a salience effect outside the predicted network.

**H5 — Physiological correction attenuates but does not abolish awareness effects.**
Adding physiological regressors (§6) will reduce cluster extent in H1–H4 but the anterior insula and dACC effects will survive.
*Confirmatory test (revised v0.8, 2026-11-25; supersedes the 2026-11-24 H2-only rule):* the TOST of §12.4 on the **H3** aware > unaware effect, in the H3 participants who have the regressor on every run in their joint model. That is N = 25 for ±PetCO2 and N = 28 for ±RVT, using the Biopac belt only. The decision rule is fixed in advance:
- *H3 supported* (uncorrected model, §3 H3): the H5 verdict is the §12.4 grid cell of the H3 TOST, read in the ROI(s) where H3 was supported. The H2 TOST (N = 25 / 27) is computed the same way and is **supportive**. If the two disagree, the disagreement is reported and H3 governs.
- *H3 not supported, H2 supported:* there is no categorical awareness effect to attenuate. The H5 verdict is the H2 TOST cell, stated as a conclusion about the confidence measure, not about awareness. The H3 TOST is reported descriptively.
- *Neither supported:* H5 is not tested. Both comparisons are reported descriptively.
PetCO2 and RVT are separate comparisons, each with its own verdict, and no correction is applied between them (unchanged from v0.2).
This is directional and consequential: PetCO2 is a **mediator** in this design, not a pure nuisance — the pacing manipulation changes end-tidal CO2 by construction. Regressing it out therefore removes part of the effect of interest for H1 in particular. Both corrected and uncorrected models will be reported side by side, and neither is designated "the" result for H1. Precedent: the IEAT found RVT associated with widespread deactivation, and that controlling respiration rate and RVT at both levels attenuated but did not eliminate its interoception deactivation effect.
*Note on the TOST target (revised v0.8):* from v0.2 to v0.7 the confirmatory target was H2, because v0.5's 23 H3 participants included only 13 with PetCO2 and 15 with RVT, too few for a usable bound (80% power needed dz ≈ 0.84–0.92). With H3 at N = 45 the subsets are 25 and 28, and the bounds match H2's. The target moves to H3 for two reasons. H3 is the authoritative awareness test, and H2's confidence effect bundles awareness with confidence, so equivalence on H2 says little about H3. The respiratory confound is also specific to the aware-vs-unaware split: detection may track trials on which breathing followed the pacer more closely. *Registered cost:* H3's per-participant estimates are noisier than H2's, so the third grid cell (*indeterminate*) is more likely for H3 than for H2. The bound is in dz units, so its power is the same; what changes is how wide it is in raw units. The corrected–uncorrected correlation is reported for both so that width can be read.

**H6 — MAIA and the ACC (conceptual extension of the IEAT; two-sided). Model A, N = 46.**
*What is carried over, and what is not.* In the IEAT, higher MAIA predicted reduced ACC deactivation on the **Active Interoception vs Active Exteroception** contrast. This task has no exteroceptive condition, so that contrast cannot be reproduced. The contrast available here is paced vs free, the analogue of the IEAT's Active Matching vs Active Interoception, and the IEAT did not relate MAIA to it. H6 therefore tests whether the IEAT's individual-differences region, defined in advance, is MAIA-sensitive on this study's control contrast. It is **not a direct replication**, and was labelled as one in v0.2 in error (corrected 2026-11-24, before any model was fitted).
*Why two-sided.* The prior study does not fix the sign on this contrast, and the two plausible mappings predict opposite signs:
- If MAIA spares the ACC during interoceptive monitoring as such, both paced and free epochs carry that monitoring, and no difference is predicted.
- If free epochs are the closer analogue of Active Interoception, sparing during free epochs predicts a **negative** MAIA effect on paced − free.
- If the sparing reflects resistance to the deactivation that exogenous pacing drives (H1), it predicts a **positive** effect.

No sign is registered, and the test is two-sided.
*Test:* MAIA total is entered as a demeaned group-level covariate on the C1 contrast (paced − free) from Model A. It is scored as in the IEAT: the mean of the six subscales, omitting Not-Worrying and Not-Distracting. Scores are in `Repo/Data/maia_scores.csv`; the scoring code is `intero2023_maia.R`. Inference is small-volume, two-sided α = .05, in a 10 mm sphere at the IEAT Results coordinate (−4, +12, +36), with the CAREN dACC mask (§8) reported as a secondary, descriptive mask. The support criterion uses the sphere only.
*Supported if* the MAIA effect survives small-volume correction in the sphere, in either direction. The sign is then reported and interpreted against the three mappings above. *Not supported if* it does not survive.
*Coverage and power:* 46 of the 47 analysable participants have complete MAIA (sub-065 has none). At N = 46, two-sided α = .05, 80% power requires |r| ≥ 0.40. At the IEAT's ROI-level |r| ≈ 0.27–0.32, power is ≈ 44–59%. **H6 is underpowered for an IEAT-sized effect, and a null is not evidence against the IEAT result**, which concerned a different contrast.
*Rationale:* a fixed a priori region and construct from an independent sample, applied to the contrast this design affords. MAIA subscale analyses (the Attention Regulation vs Noticing dissociation reported in the IEAT) are exploratory (§10).
*Administration note (registered):* this study administered the 37-item MAIA-2 with response boxes coded 1–6; scores are converted to the canonical 0–5 metric before use. The IEAT used the original MAIA; the six retained subscales are common to both.

**H7 — Between-network coupling in the attention and awareness contexts (added v0.5; both parts confirmatory). Connectivity model (§7, final subsection).**
*Rationale.* Monitoring an imposed respiratory change, and becoming aware of one, should require the salience network to exchange signal with two other networks: the one that carries the bodily input and the one that uses it. Two pairs are named.
- **SAL–SMN:** somatomotor input to salience. The IEAT's pacing effect fell on the somatomotor-network regions that H1 tests.
- **SAL–ECN:** salience-to-executive broadcast, the standard network account of access (Menon & Uddin, 2010; Sridharan, Levitin & Menon, 2008).

Networks are the whole CAREN salience (SAL), somatomotor (SMN) and executive (ECN) networks (§8). *Prior evidence is weak and is stated now:* in MindLock, the analogous SMN–SAL test was null (dz = +0.17), and coupling changes appeared only in pairs involving ECN.

**H7a — attention context (paced vs free). N = 47.** *Prediction:* coupling is greater during paced than during free epochs, for SAL–SMN and SAL–ECN.

**H7b — awareness context (aware vs unaware). N = 23, the participants with ≥ 3 epochs in each cell (the S7 set, §3 H3).** A correlation needs several epochs per cell, so H7b keeps the v0.5 rule when H3 relaxes it (v0.6). *Prediction:* coupling is greater on aware than on unaware change epochs, for the same two pairs.

*Test:* a two-sided paired t on the per-participant difference in Fisher-z coupling, Holm-corrected across the two pairs within each part. *Supported if* at least one pair shows the predicted direction and survives Holm. *Falsified if* neither does. An effect that survives Holm in the opposite direction is reported as contrary to the prediction, not reframed.

*Interpretation gate, fixed now.* A supported H7a or H7b is reported as a network-coupling effect only if it stays significant (two-sided p < .05, predicted direction) once the fMRIPrep global signal, processed identically, is partialled out of the pair correlation (§7). If it does not, it is reported as **a coupling change not separable from the global signal**. Two features of this design make that the live alternative:
- Paced breathing changes ventilation and CO2 by construction (§6), which moves the whole-brain signal.
- H1 predicts lower regional activity under pacing. On its own, that pushes measured coupling up for non-neural reasons, because shared noise becomes a larger share of the remaining variance (MindLock, Tests of Non-Neural Accounts).

The other non-neural checks (§7) are reported alongside, descriptively: co-activation, temporal amplitude, the five-network mean, and the physiological corrections. For H7a the physiological corrections use +RVT on 30 participants and +PetCO2 on 25. For H7b they use 15 and 13.

*Power.* At 80% power, the detectable dz is:
- **H7a:** 0.42 at α = .05, or 0.46 at Holm's first-step α = .025.
- **H7b:** 0.61 and 0.68.

MindLock's registered coupling effect was dz = 0.49, from about 3–4 times as many volumes per condition. Here each participant contributes about 140 volumes per condition for H7a, and medians of 67 aware and 44 unaware volumes for H7b. The within-participant estimates are therefore noisier. **H7b is powered only for a large effect, and a null H7b is not evidence against awareness-related coupling.**

*Not an individual-differences measure.* MindLock's per-participant coupling differences had no between-participant reliability (split-half −0.13 to +0.20) although each coupling estimate was reliable. No H7 difference score will be correlated with MAIA or behaviour in a confirmatory test.

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
| of which: Biopac RVT on every retained run (S2) | 61 | 30 |
| of which: RVT from Biopac or PMU fallback on every retained run (S2b) | 79 | 41 |
| of which: PetCO2 on every retained run (S1) | 51 | 25 |
| of which: both Biopac RVT and PetCO2, on every retained run | 43 | 21 |
| of which: verified pacer adherence (≥1 retained run) | — | 26 |

*(v0.4: physiological rows updated after the pipeline-bug fixes and the PMU fallback, §6. v0.3: the physiological rows are counted on the 90 post-gate runs. v0.2 carried pre-gate values here (64 / 51 / 43 runs; 30 / 24 / 20 participants; adherence 24), although the N table below already had the post-gate 22 and 26.)*

**Additional motion exclusion, applied before modelling:** a run is dropped if mean framewise displacement > 0.5 mm or if > 25% of its volumes exceed FD 0.5 mm. A participant is dropped if this leaves fewer than one usable breathing run. This criterion is fixed now and will not be adjusted after seeing results.

**The gate has been applied** (2026-08-14, using fMRIPrep's FD column only — no BOLD signal or contrast has been examined). Across the full 103-run set it fails 12 runs, one of which is sub-056 run-2, already removed on the line above. It removes exactly the same 5 participants entirely as in the registered version (sub-058, -062, -080, -085, -094); two more lose one of two runs. Gated runs are unambiguous (median 81 spike volumes vs 3 in retained runs; median FD 0.52 vs 0.18 mm). **Post-gate: 90 runs, 47 participants** (`Repo/FirstLevel/firstlevel_manifest.csv` plus `firstlevel_manifest_candidates.csv`). All four reinstated participants have belt-verified pacer adherence on their retained runs, and so enter S4. Their physiological regressors had not been built because of pipeline bugs, now fixed (§6, v0.4). sub-065, sub-082 and sub-104 now have Biopac RVT and PetCO2 on every retained run. sub-076 has Biopac RVT on run-1 and PMU RVT on run-2, so it enters S2b only, and it has no valid capnogram.

*Registered check on a paradigm-specific selection risk:* deep paced breathing could inflate FD and the gate could then preferentially remove adherent breathers. It does not: fd–adherence correlation across belt-usable runs r = −0.003 (n = 81), and gated runs show *lower* pacer adherence (median 0.230 vs 0.335).

**Run count is not an inclusion criterion.** A participant enters the analysable sample when the registered criteria — usable event timing, trigger-locked imaging, and the motion gate — are met by **at least one** breathing run, irrespective of how many runs they contribute. Where a participant contributes one run, the fixed-effects step of §8 is that run. This is stated as a general rule because it cannot otherwise be applied consistently: two participants of the original registered 43 — **sub-074 and sub-114** — already contribute exactly one breathing run, so excluding sub-082 and sub-104 for contributing one run each would apply a standard to those two participants alone rather than a registered criterion. Sensitivity analysis **S6** (§9) tests every hypothesis on two-run participants only, which removes all four together.

**Estimability rule.** A contrast counts for a participant when every cell it involves has **≥ 3 trials pooled across that participant's gate-passing runs**, and the contrast is estimable in the corresponding design (§7). This rule is stricter than matrix rank. It applies to every contrast except H3.

**H3 exception (registered v0.6; first level revised v0.7).** For C2, a participant enters with **≥ 1 aware and ≥ 1 unaware trial across their gate-passing runs**. The contrast comes from one joint first-level model of all those runs (§8), so every trial contributes. The H3 group model is precision-weighted (§8), so the low-count participants this admits are down-weighted rather than excluded. This gives N = 45. In these data, every one of the 45 also has a run containing both cells. The v0.5 rule is kept as S7 (§9) and is the rule for H7b.

**Final registered Ns (post-gate, post-correction, under the §7 models):**

| Contrast | model | N |
|---|---|---|
| C1 paced vs free (H1) | A | **47** |
| C2 aware vs unaware (H3, authoritative; precision-weighted) | A | **45** |
| C2 under the v0.5 ≥ 3-trial rule (S7; also the H7b sample) | A | 23 |
| C3 confidence parametric (H2, supportive) | B | **44** |
| C5 salience high vs low (H4) | C | **43** |
| H6 MAIA moderation | A | **46** |
| C4 wrong-direction vs aware | — | 13 — **not analysed as a planned contrast** |
| S1 ±PetCO2 (basis on all retained runs) | — | **25** (H3 25, H2 25) — TOST on H3, dz = ±0.61 (80%); H2 supportive, ±0.61 |
| S2 ±RVT, Biopac belt (all retained runs) | — | **30** (H3 28, H2 27) — TOST on H3, dz = ±0.58 (80%); H2 supportive, ±0.59 |
| S2b ±RVT, Biopac or PMU fallback (robustness) | — | **41** (H2 38, H3 39) |
| S4 verified pacer adherence | — | **26** (H2 25, H3 26) |

Three participants — **sub-066, sub-070 and sub-114** — give the same confidence rating on every retained paced trial (all 14 rated 6, all 14 rated 7, and all 7 rated 7 respectively). A constant modulator is not a regressor, so they are **excluded from H2 only** and retained for every other hypothesis. This is a property of the behaviour, not of the sample correction, and the rule is fixed now, before any model is fitted.

**Run-level exceptions in Model B (registered 2026-11-24).** The same rule applies run by run. A run whose paced confidence ratings are constant has an all-zero mean-centred modulator, so C3 is not estimable in it. Seven runs are affected:
- Five belong to the three participants already excluded from H2: sub-066 runs 1–2, sub-070 runs 1–2, and sub-114 run-1 (its only retained run).
- Two belong to participants who remain in H2: **sub-106 `task-breathing_run-1`** and **sub-112 `task-breathing_run-2`**. For these two, C3 is taken from the other run alone. The run still enters Models A and C and contributes to every other hypothesis.

The fixed-effects combination for C3 uses only runs in which the modulator varies. This is the general rule; these two runs are simply the ones it currently applies to.

H6 power (two-sided, §3): 80% for |r| ≥ 0.40 at N = 46, and ≈ 44–59% at the IEAT's ROI-level |r| ≈ 0.27–0.32. Under S6 (N = 42) the threshold is |r| ≥ 0.42. **The correction does not change the study's power situation materially and is not presented as if it does** — the case for reinstating these participants is that the inventory used to exclude them was wrong, which is a correctness argument, not a statistical rescue.

---

## 5. Design and event model

Per breathing run, the following events are modelled (already extracted and validated against a hand-checked participant, 80 assertions in `Repo/Tests/`):

- **Paced breathing epochs** (~20 s), one per change trial, split by condition
- **Free breathing epochs** (~20 s), modelled as a separate trial type
- **Rating window**, modelled explicitly rather than left in the implicit baseline
- **Valence/sub-rating events** as nuisance

Trial counts per participant (median, corrected 47-participant sample): 7 aware-correct, 3 unaware, 1 wrong-direction, 2 No Change, 14 free. v0.2 quoted 2 unaware, a figure carried over from the pre-correction sample. The design is short — 14 paced trials per participant — and that shortness is why the awareness contrast rests on few trials per participant. Under the v0.5 ≥ 3-trial rule it was estimable in only 23 participants. v0.6 admits 45 and weights them by precision (§3, §4, §8).

**Wrong-direction trials** (median 1 per participant) are modelled as a **separate regressor of no interest** so they do not contaminate the aware or unaware regressors. They are not contrasted.

**No Change trials** (median 2) are modelled as their own regressor and are not treated as correct rejections.

**There is no implicit baseline (registered 2026-11-24).** Paced, free, rating and sub-rating events follow one another almost without gaps. Across the 90 modelled runs they cover a median **94%** of the acquisition (range 63–98%). The unmodelled remainder consists mostly of ~2 s transitions. Every single-condition parameter is therefore estimated against a few seconds of gaps rather than against rest, and is not interpretable on its own. Consequences, fixed now:
1. Only **between-condition contrasts** (C1–C5 and the parametric C3) are interpreted.
2. "Deactivation" in H1 and "spared activity" in H6 mean *relative to the comparison condition*. They do not mean relative to rest.
3. **No condition-vs-rest contrast** (for example free > baseline, or paced > baseline) will be reported, as a finding or as an exploratory analysis.
4. Single-condition beta maps will not be shown.

The modelled rating windows are what make the paced and free epochs separable from the response period. Leaving them in the implicit baseline instead would only change which unmodelled period the parameters are relative to. Design rank is unaffected: all 90 designs are full-rank in every model (§7).

**Trailing events are clipped to the scan end, counted, and never dropped** (registered 2026-08-24). Any modelled event extending past `n_vol × TR` has its duration truncated to the end of the acquisition; an event whose truncated duration is ≤ 0 — one beginning after the scanner stopped — is dropped, having no data to explain. These are runs in which the scanner stopped during a final rating, not runs with bad events: dropping the whole event would discard seconds of acquired signal and silently change the design, while extending the design matrix past the data is not possible. The count is reported per run. Across the 90 modelled runs, three carry a clipped event:

| run | volumes | events clipped | event over-run |
|---|---|---|---|
| sub-052 `task-breathing_run-1` | 244 | 1 | 11.0 s |
| sub-053 `task-breathing_run-1` | 250 | 2 | 25.6 s |
| sub-053 `task-breathing_run-2` | 260 | 1 | 5.6 s |

The larger over-runs elsewhere in the event files (sub-058, sub-094: 19–67 s) belong to runs the motion gate already removed, so they never reach a model.

**The paced side of C1 averages over the paced conditions present in the run** (registered 2026-08-24). C1 contrasts `free` against the unweighted mean of whichever of `paced_aware`, `paced_unaware` and `paced_nochange` have trials in that run. The earlier fixed 1/3 weighting is undefined in 18 of 91 runs — 17 of them missing `paced_unaware`, because the participant never answered "No Change" on a change trial — and the hypothesis is a claim about paced versus free breathing, not about the 1/3 weights. Precision is unaffected: median SE multiplier `sqrt(c'(X'X)⁻¹c)` is 0.235 for both forms in the 72 runs where both are defined (corrected v0.7; earlier versions quoted 0.0827 and 0.0830 from a scaling error in the QC script, change log item 21); the change moves H1 from N = 45 to N = 47, and the two recovered participants are recovered because a cell was empty, not because their data improved. *Stated as a limitation:* the paced side is a slightly different weighted average across participants. Since the weights carry no theoretical content, this is a smaller distortion than dropping 18 runs. The per-run composition is logged in `Analysis/design_qc.csv`. `paced_wrongdir` and `paced_noresponse` remain regressors of no interest and stay out of the contrast.

---

## 6. Physiological covariates

Two independent physiological measures are available, both already built and validated:

**RVT** (respiration volume per time, Birn et al. 2006), convolved with the respiration response function (Birn et al. 2008), z-scored within run.
Validated: RVT ⊗ RRF predicts the global signal at **r = +0.223, positive in 86% of runs**, at zero lag with no lag selection. Available for 153 runs / 41 participants overall, counted across all tasks before v0.4; breathing-run coverage after v0.4 is in §4.

**Pipeline fixes (v0.4, 2026-11-25).** Three bugs had withheld valid Biopac data from the breathing runs. Fixing them recovers 7 RVT and 5 PetCO2 runs under the registered pipeline and thresholds, unchanged:
1. The physio cache-key regex could not parse a part label containing an underscore. sub-070's breathing runs had therefore been built from a 1341-s fragment recording; its derived rows were rebuilt from the correct one.
2. A hard-coded exclusion of sub-082 and sub-104, left over from the withdrawn "shifted labels" rationale, was removed.
3. Volume counts are now taken from the derivatives of record, so sub-065 and sub-076 are no longer skipped.

One knock-on outside these participants: the cohort-level cut-points that classify each epoch's *measured* breathing change moved slightly (0.526 / 1.898 → 0.544 / 2.044). This reclassifies `actual_change_physio` in 20 epochs from 15 other participants. That variable enters no registered test. `adherence_verified` is unchanged for everyone else.

**PMU fallback RVT (v0.4; robustness analysis S2b only).** Where a run has no usable Biopac belt, RVT can be taken from the Siemens PMU respiration cushion. The rule (code: `Scripts/intero2023_pmu_rvt.R`; validation: `Analysis/pmu_rvt_validation.md`):
- *Source priority.* Biopac if `rvt_rrf_z` covers ≥ 50% of volumes, otherwise PMU if it passes every gate, otherwise none.
- *Alignment.* Scan start in the PMU log = `AcquisitionTime − LogStartMDHTime + 0.62 s`. This predicts the Biopac-trigger-locked offset at r = 0.993 across the 46 runs with a reliable reference (residuals −1.1 to +2.0 s).
- *Gates.* The log covers ≥ 95% of the scan; band ratio ≥ 0.90; < 10% of samples on a 12-bit rail in the scan window; ≥ 5 kept breaths; RVT on ≥ 50% of volumes; inverted polarity.
- *Processing.* The same RVT ⊗ RRF chain as the Biopac route.

This gives 18 breathing runs, raising RVT coverage from 30 to 41 participants (H3: 15 → 21 of 23).

*Why PMU is robustness only.* On 51 dual-belt runs, PMU RVT ⊗ RRF predicts the global signal at median r = +0.34, positive in 90% of runs, against +0.22 for the Biopac belt on the same runs. But the two belts' RVT time courses agree poorly (median r = 0.17 over all 51; 0.05 on the 25 not used to develop the method). The PMU signal carries respiratory information, but it is not interchangeable with the belt the plan registered. The confirmatory comparison (S2, H5) therefore stays on the Biopac belt.

*Disclosure of the gate sequence.* Before any fallback run was scored, the saturation gate was set at < 5% of the whole PMU file. The first scoring printed every fallback run's quality metrics, including its RVT ⊗ RRF correlation with the global signal, and that gate refused 17 of 28 candidates. The dual-belt validation runs themselves reach 27% saturation, so the gate was re-derived from the validation runs only, and set at < 10% within the scan window. The fallback runs' metrics were not used to choose the cut, but they had been seen when it was chosen. The alignment-reliability filter above (reference sharpness ≥ 0.05) was also set after the alignment residuals were seen; it changes neither the +0.62 s constant nor any gate. The 7 runs at 10–15% saturation are not used. These steps involved the global signal only (the exception disclosed in §0); no task regressor was related to any imaging measure.

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
- CSF and white matter mean signals, taken from the fMRIPrep derivatives of record (regenerated v0.6; see change log item 18)
- Cosine drift terms (fMRIPrep, 128 s high-pass)
- Spike regressors for volumes with FD > 0.5 mm

**Physiological models (added in sensitivity analyses only, §9):**

- `+RVT`: RVT ⊗ RRF (1 regressor)
- `+PetCO2`: 4-lag HRF basis (4 regressors)
- `+both`

The global signal is **not** included as a nuisance regressor in any model.

**Design QC, computed before the freeze and before any model was fitted.** `Scripts/intero2023_designqc.py` and `Scripts/intero2023_twomodel_qc.py` build every design matrix and report rank, per-contrast SE multipliers, VIF among the task regressors, and the modulator/awareness correlation; outputs are committed as `Analysis/design_qc.csv` and `Analysis/design_qc_twomodel.csv` (written there from v0.8; the copies committed before v0.8 predate the v0.7 SE fix). No BOLD image is read by either.

**Connectivity model (H7; registered v0.5).** It is computed from the fMRIPrep preprocessed BOLD (MNI152NLin2009cAsym res-2), unsmoothed, and is separate from the GLM fits.

1. **Networks.** The CAREN whole-network masks `ROIs/net_SAL`, `net_SMN` and `net_ECN` are checksummed in `ROIs/roi_manifest.csv` and share no voxels. Voxel time series are denoised (step 2), and the network signal is the **median** across voxels at each volume. The mean is used as a sensitivity specification.
2. **Denoising, per run.**
   - A zero-phase Butterworth band-pass at 0.01–0.1 Hz is applied identically to the voxel time series and to every regressor (Lindquist et al., 2019).
   - The data are then regressed on (a) the base nuisance set of this section without the cosine terms, which the band-pass supersedes, and without the spike columns, which step 3 replaces by excluding the volumes; and (b) **all seven Model A task regressors** (canonical HRF + temporal derivative).
   - Removing the condition-locked mean evoked responses leaves the fluctuations around them ("background" connectivity; Al-Aidroos, Said & Turk-Browne, 2012). Without this step, epochs this short, abutting rating windows, would inflate every between-network correlation (Cole et al., 2019).
   - No global signal regression.
3. **Windows.** An epoch contributes the volumes acquired in [onset + 4 s, onset + duration + 4 s). Volumes with FD > 0.5 mm are excluded. Each network's segment is demeaned within each epoch, and the segments are concatenated across epochs and runs within a condition.
4. **Coupling.** Pearson r for each pair, condition and participant, then Fisher z.
   - **H7a:** paced (all four paced conditions) vs free.
   - **H7b:** aware vs unaware, in participants with ≥ 3 epochs in each cell (the 23 participants of the S7 set).
5. **Group test.** A paired t on the z difference. With two observations per participant this is equivalent to MindLock's `z ~ condition + (1 | participant)`.
6. **Non-neural checks, reported with every H7 result.**
   - (i) **Global signal (the interpretation gate, §3):** the fMRIPrep `global_signal` is put through the identical filter and regression and partialled out of each pair correlation within the concatenated windows. The effect is then re-tested.
   - (ii) The mean of the five CAREN networks, partialled out in the same way.
   - (iii) **Amplitude:** each network's temporal SD by condition, tested paired.
   - (iv) **Co-activation:** the network-mean Model A estimate for the same condition contrast.
   - (v) **Physiology:** steps 2–5 are repeated with RVT ⊗ RRF from the Biopac belt (the S2 participants) and with the 4-lag PetCO2 basis (S1) added to the regression.
7. **Pre-declared sensitivity specifications, all reported.**
   - **H7a, 8 specifications:** signal summary (median / mean) × window shift (4 / 6 s) × task removal (Model A canonical, or a finite impulse response basis of 2 s bins over 0–32 s for paced, free, rating and sub-rating).
   - **H7b, 4 specifications:** summary × shift. A FIR basis by awareness cell would not be estimable.
8. **Specificity (exploratory, §10).** All 10 pairs among the five CAREN networks are tested for both contrasts, FDR-corrected.

Feasibility was counted from the event timing only, before the freeze. Every participant has 7–14 paced and 7–14 free epochs, a median of about 140 volumes per condition (minimum 70). In the H7b set, the medians are 67 aware volumes (range 34–89) and 44 unaware volumes (range 27–77).

---

## 8. Group model, ROIs, and inference

**Group model:** one-sample t-test on first-level contrast images, participant as the unit. Where a participant contributes two runs, run-level contrasts are combined at a fixed-effects intermediate level before entering the group model. H3 is the one exception, described next.

**H3 group model (registered v0.6; first level revised v0.7).** H3 uses a mixed-effects model that weights participants by precision. The pipeline runs as follows:
1. **Joint first level (v0.7).** Each of the participant's gate-passing runs is scaled to percent signal change, and the runs are concatenated. One nilearn GLM is fitted with the Model A design built as follows (`build_joint_design()` in `Scripts/intero2023_h3_joint_qc.py`):
   - The task regressors and their temporal derivatives are **shared** across runs.
   - The base nuisance regressors and the intercept are **run-specific** (block-diagonal).
   - The noise model is AR(1), as for every other first level. The single sample at the junction between runs is not modelled specially.
   - The model assumes that each condition's response, and the noise variance, are the same in both runs. The behaviour supports the first: the runs are behaviourally exchangeable (§8).
   - A one-run participant's joint model is that run's Model A.
2. The aware − unaware effect, its variance and the residual degrees of freedom are taken from that model.
3. These enter **FSL FLAME 1** (`flameo --runmode=flame1`) as a one-group mean within the group brain mask. FLAME estimates the between-participant variance and weights each participant by 1 / (own variance + between-participant variance).
4. **ROI inference:** voxel-level FWE p < .05 within each mask by Gaussian random field theory, one-sided, Bonferroni ×3. The smoothness is estimated from the FLAME residuals (`smoothest`), and the resel count is each mask's.
5. **Whole-brain:** cluster-forming z > 3.09, cluster FWE p < .05 by GRF, two-sided.

Permutation inference is not used for H3, because a precision-weighted model is not exchangeable under sign-flipping. The participant weights (as a share of the total) are reported with the result. *Robustness, reported but not inferential:* each participant's mean effect and variance within each H3 mask are entered into a random-effects meta-analysis (REML), so the ROI result can be checked without the GRF assumptions. S6 and S7 for H3 use the same pipeline; S7 is additionally reported under the OLS model of v0.5. *Per-run robustness check (the v0.6 form):* the contrast is computed in each run that contains both cells, the runs are combined by fixed effects, and the result enters the same FLAME model. Agreement is reported; divergence is reported as dependence on the cross-run assumption, and is not used to choose between the forms.

**Why fixed-effects combination (decided 2026-08-14, before any imaging analysis):**

- *Run-as-replication is not available for the awareness contrast.* Within a single run the cells collapse (median 3 aware / 1.5 unaware per run); requiring ≥3 trials per cell in both runs leaves exactly 1 participant, versus the 23 that pooling yields (§4). Pooling is the condition for H3's estimability, not a trade-off against replication.
- *The runs are exchangeable.* Direction accuracy 0.800 vs 0.754 (p = 0.22), confidence 5.51 vs 5.31, P(say change) 0.726 vs 0.713 — a slight non-significant run-2 softening, no heterogeneity that would make an average misleading. Each run keeps its own GLM and baseline, so only the contrast is averaged.
- *Pooling materially improves per-participant reliability.* Cross-run individual differences are stable (direction accuracy r = 0.688, confidence r = 0.736, n = 42/51); Spearman-Brown gives the pooled estimate reliability ≈ 0.82 vs ≈ 0.69 for one run.
- *Run as a random factor is degenerate with two levels* and is not used.

**Registered reliability diagnostic (descriptive, not a gate):** for C1 (paced vs free, Model A) and C3 (confidence parametric, Model B), run-1 vs run-2 agreement of first-level contrast estimates — spatial correlation of contrast maps within participants, and cross-participant ICC of ROI means. Reported alongside the behavioural benchmark (r ≈ 0.7). H3 is excluded from this diagnostic because it is estimable in both runs for only 1 participant; this is stated now so its absence is not read as an omission later.

**A priori ROIs (registered v0.5, replacing the v0.4 Harvard–Oxford list).** Masks are defined from the **CAREN** consensus atlas of resting-state networks (Doucet, Lee & Frangou, 2019, *Hum Brain Mapp* 40:4577). Each mask is a set of network × AAL-3 regions taken from CAREN's own region breakdown, so no second atlas is intersected and no coordinate cut is chosen. The salience and somatomotor networks divide the insula on their own. Salience-network insula is centred at y ≈ +6 (anterior insula), and somatomotor-network insula at y ≈ −15 (posterior/mid insula).

| # | Mask | Definition (CAREN network ∩ AAL-3 regions) | Voxels (2 mm) | Centroid (MNI) | Used for |
|---|---|---|---|---|---|
| 1 | L anterior insula | Salience ∩ Insula_L | 588 | (−37, +5, +2) | H3, H2, H4 |
| 2 | R anterior insula | Salience ∩ Insula_R | 527 | (+38, +7, +1) | H3, H2, H4 |
| 3 | dACC / anterior MCC | Salience ∩ (Cingulate_Mid_L/R, ACC_sup_L/R) | 661 | (0, +6, +36) | H3, H2, H4; H6 secondary |
| 4 | Posterior insula, bilateral | Somatomotor ∩ Insula_L/R | 285 | (±37, −15, +12) | H1 |
| 5 | Central operculum, bilateral | Somatomotor ∩ Rolandic_Oper_L/R | 531 | (±49, −15, +14) | H1 |
| 6 | Sensorimotor cortex, bilateral | Somatomotor ∩ (Precentral_L/R, Postcentral_L/R) | 4232 | bilateral, z ≈ +49 | H1 |
| 7 | SMA, bilateral | Somatomotor ∩ Supp_Motor_Area_L/R | 406 | (0, −13, +58) | descriptive |
| 8 | **IEAT ACC sphere** | 10 mm radius at MNI (−4, +12, +36), the IEAT Results coordinate of its MAIA-moderated ACC cluster (see §1b for the IEAT's alternative Methods coordinate) | 528 | (−4, +12, +36) | H6 |

*Inference sets, Bonferroni-corrected within each:* H3, H2 and H4 use masks 1–3 (×3). H1 uses masks 4–6 (×3). H6 uses mask 8 alone, with mask 3 reported descriptively. Mask 7 enters no test.

*Construction* (`Scripts/intero2023_build_rois.py`, run with FSL's Python):
- The atlas is resampled to the fMRIPrep MNI152NLin2009cAsym res-2 grid (97 × 115 × 97, 2 mm) by sform, nearest neighbour.
- Each mask is intersected with the group brain mask, the intersection of the 90 analysed runs' fMRIPrep brain masks. Coverage is 100% for every mask.
- Masks 1–7 share no voxels with each other. The IEAT sphere shares 209 of its 528 voxels with mask 3; H6 is a separate test.
- CAREN is distributed on the FSL MNI152 2 mm grid, not in MNI152NLin2009cAsym. The template difference is small, is not corrected, and is declared here, as it was in MindLock.
- Signal dropout cannot be assessed without reading BOLD intensities. Each mask's coverage within the group second-level mask is therefore reported after the freeze, but the masks are not altered by it.
- **The masks are frozen by checksum.** `ROIs/roi_manifest.csv` holds each mask's voxel count, centroid, and SHA-256 of file and data. `ROIs/roi_sources.csv` holds the SHA-256 of every atlas file used. A mask whose checksum does not match may not be used. **The mask files are not published, by decision (2026-11-25)**, neither in the repository nor on OSF, because they are derived from a third-party atlas. The generating script and the checksums are published. Anyone with the CAREN atlas can rebuild the masks and verify each one against its registered checksum.

*Atlas-sensitivity set (reported, never used for inference):* Harvard–Oxford cortical maxprob-thr50, built by the same script:
- Anterior insula: HO insula with y > 0, split by hemisphere.
- Posterior insula: HO insula with y ≤ 0.
- dACC: the HO anterior cingulate division.
- Central operculum: the HO central opercular cortex.
- Sensorimotor cortex: HO pre- and postcentral gyri.

The H1–H4 ROI tests are repeated in these masks. Disagreement with the CAREN result is reported, not resolved.

*Network masks for H7:* the whole CAREN SAL, SMN, ECN, DMN and VIS networks (`ROIs/net_*`; 4156, 9033, 7154, 9952 and 8966 voxels), built by the same script. They share no voxels. The H3/H2/H4 masks lie inside SAL, and the H1 masks inside SMN. H7 does not use the ROI masks.

**ROI restriction for H3** (resolved 2026-08-14): H3 is tested in three masks only — left anterior insula, right anterior insula, and dACC (Bonferroni ×3). The dACC is retained alongside the insula on two grounds: the IEAT found the ACC, not the insula, to be the region tied to breath-focused interoceptive awareness; and its liminal-vs-salient distinction predicts ACC involvement for exactly the kind of salient respiratory signal this task uses. Narrowing to anterior insula alone would have encoded the heartbeat-detection literature's prediction while discarding the direct predecessor's.

**Smoothing** (resolved 2026-08-14): 6 mm FWHM, matching the IEAT for replication comparability.

**Thresholds:**

- **Whole brain:** cluster-forming p < 0.001 uncorrected, cluster-level FWE p < 0.05
- **ROI:** small-volume FWE p < 0.05 within each ROI, Bonferroni-corrected across the number of ROIs tested for that hypothesis. For H3, voxel FWE comes from GRF on the FLAME result (see above).
- **Sidedness:** H1 is one-sided (paced < free) in its ROIs. H2 is one-sided (positive). H3 and H4 are one-sided in the stated direction. H6 is **two-sided**. H7a and H7b are two-sided paired tests on Fisher z, Holm-corrected across their two network pairs, with the predicted direction required for support (§3). All whole-brain tests are two-sided.
- No uncorrected results will be reported as findings; any uncorrected map shown will be labelled descriptive

**Power.** Smallest detectable Cohen's d at 80% power (one-sample t):

| N | hypothesis | α = 0.05 (ROI) | α = 0.001 (cluster-forming) |
|---|---|---|---|
| 47 | H1 (Model A) | 0.42 | 0.64 |
| 45 | **H3 (Model A; equal-weight approximation)** | **0.43** | **0.66** |
| 44 | H2 (Model B) | 0.43 | 0.67 |
| 43 | H4 (Model C) | 0.44 | 0.67 |
| 41 | ±RVT, PMU-extended (S2b) | 0.45 | 0.69 |
| 30 | ±RVT, Biopac (S2) | 0.53 | 0.83 |
| 26 | adherence (S4) | 0.57 | 0.91 |
| 25 | ±PetCO2 (S1) | 0.58 | 0.93 |
| 23 | S7 (v0.5 H3 set) and H7b | 0.61 | 0.98 |

**H3's row assumes equal per-participant precision. Its effective power is lower, because participants with one or two unaware trials contribute little weight.** No closed-form figure is registered for that loss. v0.5's H3 at N = 23 needed d ≈ 1.0 at whole-brain thresholds, and H3 remains ROI-restricted. It is nonetheless the *authoritative* awareness test (§3): being the better-identified contrast and being the better-powered one are different properties, and this document deliberately assigns primacy on identification rather than on power. The consequence is stated in advance in §3 and repeated here — **a null H3 is not to be reinterpreted afterwards, and a null H3 beside a positive H2 is not evidence that confidence rather than awareness drives the effect**, because that comparison is confounded with per-participant precision (§3).

---

## 9. Registered sensitivity analyses

Each is fully specified now and will be reported whatever the outcome. Where a sensitivity analysis names a hypothesis, it is refitted in **that hypothesis's own model** (§7) — H1/H3/H6 in Model A, H2 in Model B, H4 in Model C — with the physiological or subsetting change applied identically to all of them. S6 is ordered last because it was added in v0.2; the numbering of S1–S5 is unchanged so that references to them in the v0.1 record still resolve.

**S1 — ±PetCO2 (the primary methodological question, N = 25).**
Fit H1–H4 with and without the PetCO2 basis in the 25 participants who have it on all retained runs (H1 25, H2 25, H3 25, H4 22, H6 24; the H5 TOST is on H3 when H3 is supported, §3 H5). Compare (a) peak z and cluster extent in each ROI, (b) whether each hypothesis's conclusion changes. Because this is a within-participant model comparison, it is better powered than the between-group comparison originally considered.

**S2 — ±RVT, Biopac belt (N = 30).** As S1, with the Biopac RVT regressor (H1 30, H2 27, H3 28, H4 26, H6 29). This is the confirmatory RVT comparison, and the H5 TOST uses it (H3 N = 28, dz = ±0.58; H2 N = 27, ±0.59).

**S2b — ±RVT with the PMU fallback (robustness; N = 41).** As S2, but runs without a usable Biopac belt take RVT from the Siemens PMU cushion under the rule in §6 (H1 41, H2 38, H3 39, H4 37, H6 40). The TOST is also computed here (H3 N = 39, dz = ±0.49; H2 N = 38, ±0.49) and reported as robustness. If S2 and S2b disagree, the disagreement is reported. S2 remains the registered result.

**S3 — CO2-available vs CO2-unavailable participants.** Test whether the 25 participants with usable capnography differ from the remainder on the uncorrected H2 contrast (H2 rather than H3, because the subgroup comparison needs the larger sample to be worth running at all). **Registered caveat:** CO2 availability is near-collinear with cohort and protocol version (the capnograph outage was an unbroken date block from 2023-10-12 to 2024-01-18), so this comparison **cannot isolate CO2** and will be reported as a cohort check only, not as evidence about physiological correction.

**S4 — Adherence conditioning (N = 26; H2 25, H3 26).** Repeat **H3 and H2** restricted to participants with verified pacer adherence on at least one retained run. Adherence is a run-level property (per-epoch classification is underpowered at AUC 0.66–0.69), and this is registered as robustness, not as the primary analysis.

**S6 — Two-run participants only (N = 43).**
Every hypothesis is additionally tested on the subset contributing **two** breathing runs. The subset is defined by run count alone, so it removes **sub-074, sub-114, sub-082 and sub-104** together — two participants from the original registered sample and two of the reinstated ones — and cannot be used to remove the reinstated participants selectively.

| | participants | H1 | H3 | H2 | H4 | H6 |
|---|---|---|---|---|---|---|
| **Primary** (corrected sample) | 47 | 47 | 45 | 44 | 43 | 46 |
| **S6** (two-run participants only) | 43 | 43 | 42 | 41 | 43 | 42 |

*Interpretation fixed in advance:* agreement between the primary and S6 results is reported as such; divergence is reported as a limitation of the single-run participants' precision and **not** used to select which analysis to present.

**S7 — The v0.5 H3 rule (N = 23; registered v0.6).** H3 is repeated in the 23 participants with ≥ 3 aware and ≥ 3 unaware trials pooled across runs (list in §3). It is fitted twice: with the FLAME model of §8 and with the OLS one-sample model of v0.5. *Interpretation fixed in advance:* agreement with the primary H3 is reported as such. Divergence is reported as a dependence on the inclusion rule and is not used to select which result to present. S7 carries the v0.5 selection toward less frequent detectors (§3), so it is not a cleaner test, only a stricter one.

**S5 — Protocol version.** v1 (n = 221 trials) and v2 (n = 1218 trials) differ in rating structure. Include protocol version as a group-level covariate; report whether the **H3 and H2** conclusions change.

---

## 10. Explicitly exploratory (will be labelled as such)

- Any brainstem or thalamic effect (coverage and physiological noise are unfavourable at 2 mm with this protocol)
- Free-breathing epochs related to measured respiratory change — the per-epoch respiration estimate has a noise floor (sd of log period ratio = 0.399) comparable to the imposed signal, so this is underpowered by construction
- Arousal as a parametric modulator
- Any relationship between neural measures and questionnaire/screener variables beyond the registered MAIA test (H6); MAIA subscale dissociations (Attention Regulation vs Noticing, per the IEAT) are exploratory
- Direction-specific effects (Faster vs Slower), given the behavioural asymmetry (0.614 vs 0.511)
- ACC-seeded PPI connectivity with the dorsal attention network during paced vs free breathing (conceptual extension of the IEAT's H4)
- **Graded awareness (Model A-g; added v0.6).** A variant of Model A in which aware and unaware change epochs form one regressor. Its parametric modulator is g = +confidence on aware trials and −confidence on unaware trials (confidence 1–7), mean-centred within run. A low-confidence aware trial therefore sits close to the low-confidence unaware trials, and a confident aware trial sits furthest from them. Wrong-direction, No Change, free, rating and sub-rating are modelled as in Model A. It is fitted in the joint two-run form and tested in the H3 masks with the H3 group model (N = 45). The shared modulator column is mean-centred within each run. Unlike Model B's confidence modulator, g is anchored to the categorical outcome, so it asks whether confidence grades the awareness contrast rather than standing in for it. It cannot rescue or overturn H3. Two variants are also reported: confidence on aware trials only, and rank-transformed g.
- Specificity of H7: all 10 CAREN network pairs for paced vs free and aware vs unaware, FDR-corrected (§7, connectivity step 8)
- Any connectivity analysis other than H7, including dynamic or seed-based methods

---

## 11. Deviations protocol

Any departure from this document will be recorded in `Repo/Preregistration/deviations.md` with the date, the reason, and whether it was made before or after the affected result was seen. Analyses added after seeing imaging results will be reported as exploratory without exception.

---

## 12. Open decisions for co-author input

1. ~~**Combining runs.**~~ **Resolved 2026-08-14:** fixed-effects within participant, with a registered descriptive run-reliability diagnostic. Rationale in §8.
2. ~~**Smoothing kernel.**~~ **Resolved 2026-08-14:** 6 mm FWHM, matching the IEAT.
3. ~~**H3 ROI restriction.**~~ **Resolved 2026-08-14:** bilateral anterior insula + dACC (three masks, Bonferroni ×3). The dACC is retained on the IEAT's evidence; see §8.
4. ~~**Formal equivalence test for H5.**~~ **Resolved 2026-08-14.** Registered criterion:
   - *Test:* TOST (two one-sided tests, α = .05 each) on the **paired within-participant difference** between the corrected and uncorrected ROI contrast estimates, extracted per participant from each of the three primary ROIs. *Target (v0.8):* the H3 aware > unaware effect from the joint first level, with the H2 confidence effect as the supportive test; the decision rule is in §3 H5. *Weighting:* the paired difference is equal-weighted for both, although H3's group model is precision-weighted. Weighting it would need the covariance of two estimates from nearly identical models, which FLAME does not supply. The equal-weighted test is the more conservative choice. *Significance of the corrected effect (the grid's second cell):* the same one-sided ROI test as the hypothesis itself, fitted on the corrected model in the same subset (FLAME 1 for H3).
   - *Bounds, set by what the samples can support at 80% power:* for **H3** (confirmatory, v0.8), **dz = ±0.61** for ±PetCO2 (N = 25) and **dz = ±0.58** for ±RVT on the Biopac belt (N = 28); S2b (N = 39) uses ±0.49. For **H2** (supportive), ±0.61 (N = 25) and ±0.59 (N = 27); S2b (N = 38) ±0.49. Simulated power at the H3 bounds is 0.81, 0.82 and 0.83 (200,000 draws, true difference 0). ±0.57 at N = 28 gives 0.807; the bound is rounded outward for margin, as v0.3 did for H2. Tighter bounds are not honestly supportable: dz = 0.3 would require N ≈ 90. *(v0.4, 2026-11-25: Ns are now the H2 participants in each subset, after the physio pipeline fixes. v0.3 stated dz = ±0.65 / ±0.60 at the whole-sample N = 22 / 26, although the H2 ∩ RVT N was 24. Ns corrected 2026-11-24. v0.2 stated these bounds at the pre-gate N = 24 / 30, with 82% / 84% power; at the true post-gate Ns the same bounds give 80%, by simulation.)*
   - *Interpretation grid, fixed in advance:* TOST passes → correction leaves the effect **statistically equivalent**; TOST fails but the corrected effect remains significant in the ROI → **attenuated but present**; TOST fails and the corrected effect is not significant → **indeterminate** (not "abolished" — the design cannot distinguish abolition from power loss at these bounds).
   - *Context for the bound's meaning:* the paired SD is the SD of the difference between two estimates from nearly identical models. If the two estimates correlate at r = 0.95, dz = 0.61 on the difference corresponds to ≈ 0.19 between-participant SD of the estimate itself — a materially tight bound despite the modest N. The observed corrected–uncorrected correlation will be reported alongside the TOST so the bound can be read in raw units.
5. ~~**MAIA scoring check.**~~ **Resolved 2026-08-14, updated 2026-08-24:** MAIA-2 complete (all 37 items) for all 58 scored participants. Under the corrected sample, 46 of the 47 analysable participants are covered — sub-065 has none. Scores in `Repo/Data/maia_scores.csv`. **H6 finalised at N = 46.**
6. ~~**Imaging inventory.**~~ **Resolved 2026-08-20:** sub-065, sub-076, sub-082 and sub-104 reinstated; sub-051 and sub-056 named as excluded. Rationale in §4 and `Preregistration/section4_correction_draft.md`; evidence in `VERIFICATION_REPORT.md`.
7. ~~**One model or two.**~~ **Resolved 2026-08-24:** three models (§7). Awareness and confidence are never in the same design. Rationale in `Preregistration/methods_decisions.md` decision 2.
8. ~~**Which awareness test is authoritative.**~~ **Resolved 2026-08-24:** the categorical contrast H3, with H2 supportive, and with the power cost and the H3 subsample's selection both registered in advance (§3, §8).
9. ~~**C1's paced side, and trailing events.**~~ **Resolved 2026-08-24:** average over the paced conditions present in the run; clip trailing events to the scan end and log the count (§5).
10. ~~**H6 framing.**~~ **Resolved 2026-11-24:** H6 is a two-sided conceptual extension, not a direct replication (§3). See `methods_decisions.md` decision 5.
11. ~~**Baseline.**~~ **Resolved 2026-11-24:** there is no implicit baseline. Only between-condition contrasts are interpreted, and no condition-vs-rest contrast will be reported (§5).
12. ~~**H5 and H3.**~~ **Resolved 2026-11-24; revised v0.8 (2026-11-25):** the confirmatory H5 test is the TOST on H3 when H3 is supported, with H2's TOST supportive. If H3 is not supported, H5 is read from H2 and stated about the confidence measure only (§3 H5). See `methods_decisions.md` decision 19.
13. ~~**Model B run exceptions.**~~ **Resolved 2026-11-24:** C3 is combined only over runs whose confidence modulator varies. This affects sub-106 run-1 and sub-112 run-2 (§4).
14. ~~**Physiological coverage.**~~ **Resolved 2026-11-25:** three pipeline bugs are fixed (§6). PMU fallback RVT is registered for robustness analysis S2b only. The confirmatory RVT comparison uses the Biopac belt. See `methods_decisions.md` decisions 10–11.
15. ~~**ROI masks and H1's inference set.**~~ **Resolved 2026-11-25:** CAREN network × region masks (§8), built and checksummed before the freeze. H1 uses posterior insula, central operculum and sensorimotor cortex (×3). Harvard–Oxford is the atlas-sensitivity set. See `methods_decisions.md` decision 12.
16. ~~**Between-network coupling.**~~ **Resolved 2026-11-25:** H7a (paced vs free) and H7b (aware vs unaware) are confirmatory, for SAL–SMN and SAL–ECN, Holm ×2 within each, with a global-signal interpretation gate (§3, §7). See `methods_decisions.md` decision 13.
17. ~~**H3 inclusion.**~~ **Resolved 2026-11-25:** every participant with both cells in at least one run (N = 45), precision-weighted with FSL FLAME 1. The v0.5 ≥ 3-trial rule is sensitivity analysis S7 and the H7b rule. See `methods_decisions.md` decision 14.
18. ~~**Tissue confounds.**~~ **Resolved 2026-11-25:** `csf` and `white_matter` are regenerated from the derivatives of record. See `methods_decisions.md` decision 15.
19. ~~**H3 first level.**~~ **Resolved 2026-11-25:** one joint model of both runs, with shared task and run-specific nuisance regressors. The v0.6 per-run form is a robustness check. See `methods_decisions.md` decision 17.
20. ~~**Design-QC standard errors.**~~ **Resolved 2026-11-25:** scaling error and physio columns fixed; QC re-run; no verdict or N changed. See `methods_decisions.md` decision 18.
21. ~~**H5 target.**~~ **Resolved 2026-11-25:** the TOST moves to H3, conditional on H3 being supported (item 12; §3 H5). See `methods_decisions.md` decision 19.
