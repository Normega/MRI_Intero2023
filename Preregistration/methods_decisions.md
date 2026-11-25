# Registered methods decisions

**Status:** decided 2026-08-24, pre-freeze. Nothing here has been estimated on
BOLD data. Every number below comes from design matrices, event sheets and
behavioural counts only — `Scripts/intero2023_designqc.py` and
`Scripts/intero2023_twomodel_qc.py`, outputs in `analysis/design_qc.csv` and
`analysis/design_qc_twomodel.csv`.

This document records the choices that `VERIFICATION_REPORT.md` §7 identified
as blocking the freeze, the evidence each rests on, and the exact prereg text
each one changes. It is the amendment log: the prereg body is edited to match,
and this file preserves *why*.

Sample throughout: **90 runs / 47 participants** — the 84 gate-passing
registered runs plus the 6 runs of the four reinstated participants (decision
4). sub-056 remains excluded.

---

## Decision 1 — H1's paced side averages over the paced conditions present in the run

**Decided:** `paced_gt_free` is the contrast of **free** against the unweighted
mean of the paced conditions **actually present in that run** — drawn from
`paced_aware`, `paced_unaware`, `paced_nochange`. `paced_wrongdir` and
`paced_noresponse` remain regressors of no interest and stay out of the
contrast, as registered.

**Why the registered form fails.** The prereg fixes the paced side at a 1/3
weight on each of the three conditions. That contrast is undefined in any run
where one of them has no trials, which is **18 of 91 runs** — 17 of them
missing `paced_unaware`, because the participant never answered "No Change" on
a change trial. Those runs are not bad data; they are runs in which one cell of
a fixed weighting happens to be empty. The registered contrast is not a claim
about the 1/3 weights, it is a claim about paced breathing versus free
breathing, and the fixed form implements that claim only by accident of trial
counts.

**What it costs:** nothing in precision. Median SE multiplier
`sqrt(c'(X'X)^-1 c)` is **0.0827** for the adaptive form and **0.0830** for the
fixed form in runs where both are defined.

**What it buys:** H1 goes from **N = 45** to **N = 47**, and the two
participants recovered (sub-070, sub-114) are recovered because a cell was
empty, not because their data improved.

**Stated as a limitation:** the paced side is a slightly different weighted
average in different participants. Given that the weights carry no theoretical
content, this is a smaller distortion than dropping 18 runs. Runs where only
one paced condition is present contribute a single-condition contrast; the
per-run composition is logged in `analysis/design_qc.csv`
(`se_n_paced_conditions`).

**Prereg edit:** §5/§7 contrast table — replace the fixed 1/3 specification of
`paced_gt_free` with the adaptive definition, and add the per-run composition
to the reported outputs.

---

## Decision 2 — Two first-level models: awareness and confidence are never in the same design

**Decided.** Two first-level models are fitted per run on the same
preprocessed data. A **third**, Model C, was added while applying this decision
to the prereg body: §7 as registered contained no salience-split regressor, so
H4 had nothing to be tested on at all. Model C splits the paced *change* epochs
by salience (high/low) and by nothing else — salience is not crossed with
awareness, because the median participant has 4 high-salience and 8
low-salience change trials and a 2 × 2 split would leave cells of about one
trial. Model C is rank-deficient in 0 of 90 runs and gives H4 **N = 43**.

**Model A — authoritative.** The registered condition boxcars —
`paced_aware`, `paced_unaware`, `paced_wrongdir`, `paced_nochange`, `free`,
`rating`, `subrating` — SPM canonical HRF + temporal derivative, no parametric
modulator, all confound-sheet columns, `drift_model=None`. Model A yields
**H1** (`paced_gt_free`, decision 1), **H3** (`aware_gt_unaware`), **H4**
(salience), and **H6** (MAIA moderation of H1).

**Model B — supportive.** The paced epochs collapse to a single `trial`
regressor carrying the mean-centred confidence parametric modulator; `free`,
`rating`, `subrating` and the confounds are unchanged. Model B yields **H2**
(`confidence_mod`) and nothing else.

**Scope of `trial`:** all paced epochs, *including* `paced_nochange` — i.e.
exactly the epochs the registered parametric modulator was defined on. Two
reasons. It keeps Model B's modulator identical in scope to the registered one,
so the split is a change of design only, not a change of what confidence
modulates; and restricting `trial` to change trials would drop a further
participant (H2 N = 43 rather than 44, losing sub-112, whose confidence varies
only once no-change trials are counted).

### Why the split

The registered single model puts the awareness boxcars and the confidence
modulator in one design, where they are collinear by construction — the
modulator lives on the same epochs the awareness regressors partition.
Measured, in the registered design:

| | value |
|---|---|
| median absolute r(`conf_pm`, aware − unaware direction) | **0.437** |
| runs with absolute r > 0.5 | **26 / 70** |
| worst run | **0.942** |
| **rank-deficient designs** | **9 / 91** |

All 9 rank deficiencies are caused by `conf_pm` alone: removing it restores
full rank in every one. Under the split, rank deficiency is **0 / 90 in Model
A and 0 / 90 in Model B**. This is not a tuning improvement; it is the removal
of a defect that made nine designs unfittable as written.

### What the split does *not* fix — stated for the record

Separating the regressors does not separate the constructs. Confidence and
awareness covary in the behaviour: the within-participant point-biserial
`r(confidence, aware)` has median **0.316** (IQR 0.126–0.570), positive in
**37 of 43** participants, with aware trials rated a median **0.50** rating
points higher than unaware ones. Model B's confidence effect therefore
**absorbs** the awareness effect rather than being adjusted for it. The
registration states this plainly: **Model B's confidence estimate is not
independent of awareness, and no claim of independence is made from it.** A
confidence effect surviving in Model B is consistent with an awareness effect,
a graded-confidence effect, or both, and the two cannot be told apart within
Model B.

Concretely: any interpretation of H2 must be conditioned on H3. If Model A
shows an awareness effect in the same voxels where Model B shows a confidence
effect, the parsimonious reading is one effect, not two.

### The hierarchy, and its cost

The registration designates **aware vs. unaware (Model A, H3) authoritative**
and **confidence parametric (Model B, H2) supportive**. This inverts the
original prereg, which made H2 primary and H3 secondary. The rationale is that
the categorical contrast is a within-participant comparison of two trial
classes, whereas the parametric effect is confounded with that same
comparison — a supportive measure cannot arbitrate its own confound.

**The cost is recorded here rather than discovered later.** Under the registered
≥3-trials-per-cell rule, the authoritative test runs on **N = 23** and the
supportive test on **N = 44**. The 23 are not a random subset: their pooled
median unaware count is **4**, versus **2** for the 24 participants who do not
reach the threshold (median unaware *share* of change trials 0.36 vs 0.19).
The H3 sample therefore **over-represents participants who detected changes
less often**. Two consequences are registered in advance:

1. **H3's estimate does not generalise to the full sample.** It is an effect in
   participants with enough unaware trials to estimate it, which is by
   definition a poorer-detecting subgroup. This is stated wherever H3 is
   reported, and the H3 subsample's detection rate is reported next to it.
2. **A null H3 alongside a positive H2 is not evidence that confidence, not
   awareness, drives the effect.** At N = 23 versus N = 44 the comparison is
   confounded with power. This specific inference is ruled out in advance.

The 23 participants are fixed now and listed in
`analysis/design_qc_twomodel.csv` so the subsample cannot drift:

> sub-057, sub-060, sub-061, sub-064, sub-065, sub-066, sub-067, sub-068,
> sub-071, sub-073, sub-076, sub-077, sub-078, sub-084, sub-100, sub-102,
> sub-106, sub-108, sub-112, sub-120, sub-130, sub-132, sub-134

### Participants with no confidence variance

`sub-066` (all 14 paced trials rated 6), `sub-070` (all 14 rated 7) and
`sub-114` (all 7 rated 7) give a constant modulator, which is not a regressor.
They are **excluded from H2 only** and retained for every other hypothesis.
This rule is fixed now, before any model is fitted.

### Resulting Ns

| hypothesis | model | N |
|---|---|---|
| H1 `paced_gt_free` (adaptive) | A | **47** |
| H3 `aware_gt_unaware` | A | **23** |
| H6 MAIA × H1 | A | **46** |
| H2 `confidence_mod` | B | **44** |
| H4 salience high > low | C | **43** |

**Prereg edit:** §5 gains the two-model specification; §7 assigns each contrast
to its model; §3 swaps the H2/H3 primacy and carries the two registered
inference constraints above; the N table is replaced.

---

## Decision 3 — Trailing events are clipped to the scan end, counted, never dropped

**Decided:** any modelled event extending past `n_vol × TR` has its duration
truncated to the acquisition end. Events whose truncated duration is ≤ 0 —
i.e. those beginning after the scanner stopped — are dropped, since they have
no data to explain. The count is reported per run.

**Why:** these are runs where the scanner stopped during a final rating, not
runs with bad events. Dropping the whole event discards seconds of real,
acquired signal and silently changes the design; extending the design matrix
past the data is not possible. Truncation keeps every acquired volume
explained by the condition that actually occupied it.

**Extent — small and now on the record.** Across the 90 modelled runs, **3
runs** carry a clipped event:

| run | volumes | events clipped | events over-run by |
|---|---|---|---|
| sub-052 `task-breathing_run-1` | 244 | 1 | 11.0 s |
| sub-053 `task-breathing_run-1` | 250 | 2 | 25.6 s |
| sub-053 `task-breathing_run-2` | 260 | 1 | 5.6 s |

The larger over-runs in the events files (sub-058, sub-094: 19–67 s) belong to
runs already removed by the registered motion gate, so they never reach a
model.

**Prereg edit:** §5 event model gains the truncation rule; the per-run clipped
count joins the reported design QC.

---

## Decision 4 — sub-065, sub-076, sub-082 and sub-104 are reinstated

**Decided:** all four are admitted, on the grounds set out in
`Preregistration/section4_correction_draft.md` and `VERIFICATION_REPORT.md` §6,
together with the single-run decision rule and the two-run-only sensitivity
analysis already drafted there. **That analysis is registered as S6, not S3**:
the draft named it S3, which collides with the existing S3 (CO2-available vs
CO2-unavailable participants). S1–S5 keep their numbers so that references to
them in the v0.1 record still resolve. sub-056 stays excluded — its run-2 fails the registered
motion gate (FD 0.471, 26.7 % of volumes spiking) and its run-1 over-runs the
acquisition by 28 s.

**Why the four were missing:** the prereg's imaging inventory predates the
2026-08-06 re-conversion. `Data/qc_ledger.csv` still records sub-065
`task-localizer_run-1` at 167 volumes; it is 196 on disk, and that single row
is the only volume-count disagreement between the ledger and this dataset
across 222 checked runs. sub-076 had no imaging rows at all. sub-082 and
sub-104 were excluded on a stale slot-mapping rationale that this dataset
contradicts: Run1 → `localizer_run-1` (196), Run3 → `breathing_run-1` (242),
Run5 → `authenticity` (172), verified directly.

**This is a sample change, not a consistency fix**, and it is the one item that
had to be settled before the freeze rather than after. Their sheets are in
`FirstLevel/` with a separate manifest,
`FirstLevel/firstlevel_manifest_candidates.csv`; the registered
`firstlevel_manifest.csv` is untouched, so the reinstatement stays an explicit
amendment rather than a silent absorption.

**Blindness:** the decision rule in the §4 draft uses only design/timing
evidence and contrast-free image quality — never a contrast estimate. sub-082's
H1 SE multiplier is 0.0836 and sub-104's 0.0717, both inside the range of the
registered sample (8 and 10 registered participants respectively are less
precise); tSNR 30.2 and 29.1 (68th and 46th percentile); both have PEPOLAR
distortion correction.

**Prereg edit:** §4 exclusion list, sample cascade and N table are replaced by
`section4_correction_draft.md`, including sensitivity analysis S6. Under the
two-model split the corrected Ns are one to two higher than that draft's table,
which predates it: H1 47 (draft 46), H2 44 (42), H3 23 (23), H4 43 (43), H6 46
(45). S6 itself gives 43 / 23 / 41 / 43 / 42 on 43 participants.

---

# v0.3 decisions (2026-11-24, pre-freeze)

These five came out of an independent re-check of every registered number against
`FirstLevel/` event sheets, the two manifests, `Data/breathing_trials.csv` and
`Data/maia_scores.csv`, and a re-reading of the IEAT paper (Farb, Zuo & Price,
2023). No BOLD image was read. Every registered N in §4 and S6 was
reproduced exactly: 90 runs, 47 participants, the motion gate failing 12 of 103
runs, H1 47, H3 23 (identical ID list), H2 44, H4 43, H6 46, and S6
43 / 23 / 41 / 43 / 42. The decisions below concern what those numbers support.

## Decision 5 — H6 is a two-sided conceptual extension

**Evidence.** The IEAT's MAIA covariate was entered on the Active Interoception
vs Active Exteroception contrast (Methods, "Hypothesis testing"). Its Results
read: "greater MAIA scores were associated with reduced deactivation during
Active Interoception relative to Active Exteroception". This task has no
exteroceptive condition. H6 as written in v0.2 tested MAIA on paced − free,
which is the analogue of the IEAT's Matching vs Interoception contrast. The IEAT
never related MAIA to that contrast. So v0.2's justification for a one-sided
test ("the direction is fully specified by the prior study") does not hold.
Under the mapping §1b itself uses (free ≈ Active Interoception), the IEAT effect
would, if anything, predict a *negative* MAIA effect on paced − free. The
handoff had specified a positive one.

**Decided.** H6 is relabelled a conceptual extension and tested two-sided in the
same a priori sphere. The power cost is stated: |r| ≥ 0.40 for 80% at N = 46,
≈ 44–59% at the IEAT's ROI-level |r| ≈ 0.27–0.32 (IEAT Table 9). The coordinate
is kept at the IEAT Results value (−4, +12, +36, k = 1205, the cluster it carried
forward as a PPI seed). The paper's Methods give a different dACC description
(k = 794, (−4, +58, +14)); that discrepancy is recorded in §1b rather than
resolved by choice.

## Decision 6 — No implicit baseline: only between-condition contrasts

**Evidence.** Modelled events (paced, free, rating, sub-rating) cover a median
**94.3%** of each run's acquisition across the 90 runs (min 63.2%, max 98.0%).
The remainder is mostly ~2 s transitions. Every single-condition parameter is
therefore relative to a few seconds of gaps.

**Decided.** Only between-condition contrasts are interpreted. "Deactivation"
(H1) and "spared" (H6) are relative to the comparison condition. No
condition-vs-rest contrast is reported, even as exploratory, and single-condition
beta maps are not shown. This rules out, before any result exists, the most
natural post-hoc workaround for H6 (a free-vs-rest contrast standing in for the
IEAT's interoception contrast).

## Decision 7 — Physiological-subset Ns corrected to post-gate values

**Evidence.** Counted on the 90 post-gate runs: RVT on every retained run, 26
participants (54 runs); PetCO2, 22 (46 runs); both, 18 (38 runs); verified pacer
adherence on ≥1 retained run, 26. v0.2's §4 cascade, §8 power table, §9 (S1, S2,
S3, S4) and §12.4 carried the pre-gate values (30 / 24 / 20 / 24), while the §4
N table and the handoff had 22 / 26. For adherence, `Data/prereg_sample.csv`
gives 24, but that file predates the correction. Against the retained runs, the
four reinstated participants all carry `adherence_verified` and the gated
sub-062 and sub-094 drop out, which gives 26. TOST power by simulation (paired,
α = .05 per side, true difference 0): N = 22 needs dz ≥ 0.65 and N = 26 needs
dz ≥ 0.60 for 80%.

**Decided.** The bounds are dz = ±0.65 (PetCO2) and ±0.60 (RVT) at 80% power. The
0.59 in the v0.2 N table and OSF drafts was at the 80% boundary and is rounded to
the defensible value.

## Decision 8 — H3 under physiological correction is descriptive

**Evidence.** Of H3's 23 participants, 12 have PetCO2 on all retained runs and 14
have RVT. 80% TOST power at those Ns needs dz ≈ 0.92 and 0.84.

**Decided.** The confirmatory H5 test is the TOST on H2 (already the case in
v0.2). The H5 hypothesis text no longer predicts that H3's effects "will
survive". The corrected-vs-uncorrected H3 comparison is reported descriptively,
with no survival or abolition conclusion.

## Decision 9 — Model B run exceptions and minor corrections

**Model B.** Seven runs have constant paced confidence within the run, so the
mean-centred modulator is all zeros there. Five belong to the three H2-excluded
participants. Two do not: **sub-106 run-1** and **sub-112 run-2**. For C3, these
two participants are combined over their other run only. The general rule is
that the C3 fixed-effects step uses only runs in which the modulator varies.

**Minor:**
- H2's falsification criterion now includes the dACC, matching its inference set.
- H4 gains an explicit inference set (the H3 masks), a support criterion and a
  null outcome.
- §5 trial medians use the corrected sample (unaware 3, not 2).
- The §2 direction-accuracy figure (0.778) is labelled as conditional on
  reporting a direction (0.562 over all change trials).
- §1b's description of the IEAT insula effect follows the IEAT Results
  (bilateral) rather than its abstract (right).

---

# v0.4 decisions (2026-11-25, pre-freeze)

## Decision 10 — Physio pipeline bugs fixed; recovered runs adopted

**Evidence.** A run-by-run audit of all 90 analysable breathing runs (no BOLD
read; fMRIPrep `global_signal` used only as the disclosed validity check) found
three bugs that had withheld valid Biopac data. Code fixes and rebuilt tables
are in commit `4aae7cd`:

1. **The cache-key regex** `_part-([^_]*)_physio` could not parse sub-070's
   part label `_A`. The scripts silently read a 1341-s fragment recording
   instead, which gave run-1 its "40% coverage" and left run-2 empty. sub-083
   (`_PartA`; not in the analysable set) had the same problem.
2. **A stale hard-coded exclusion** blocked sub-082 and sub-104 in the RVT, RRF
   and covariate scripts.
3. **Missing volume counts:** sub-065 and sub-076 had none, because the old
   imaging snapshot lacked them. Counts now come from the derivatives of record.

Recovered under the registered pipeline and thresholds, otherwise unchanged:
RVT for sub-065 r1/r2, sub-076 r1, sub-082 r1, sub-104 r1 and sub-070 r1/r2.
PetCO2 for sub-065 r1/r2, sub-082 r1, sub-104 r1 and sub-070 r2. sub-070 r1
PetCO2 fails on the correct recording (swing 0.56 < 0.8) and stays out.
Biopac-RVT participants rise from 26 to 30, PetCO2 from 22 to 25, and both
from 18 to 21.

**Verification.** Re-running the unfixed code in the same environment reproduces
HEAD for every table that doesn't depend on the global signal, and 152/153
Biopac regressors byte-for-byte. The exception is sub-065 localizer run-1,
whose re-conversion took it from 167 to 196 volumes. The 25 changed confound
sheets gained physio columns only; all other columns are textually identical
and no existing physio column changed. Outside the affected participants
exactly one thing moves: the cohort cut-points for `actual_change_physio`
shift, reclassifying 20 epochs in 15 participants. That variable enters no
registered test.

**Not adopted:** the 2 PetCO2 runs that miss the swing floor by < 0.01, and the
2 that would need the capnogram rate floor lowered. Registered thresholds stand.

**Also corrected:** the TOST N was quoted as the whole-sample physio N, but the
test is defined on H2. It is now the H2 ∩ subset N: PetCO2 25 (dz ±0.61) and
Biopac RVT 27 (dz ±0.59).

## Decision 11 — Siemens PMU RVT as a fallback source, robustness only (S2b)

**Evidence.** `Analysis/pmu_rvt_validation.md` (commit `bd88f7f`).
- *Alignment:* `AcquisitionTime − LogStartMDHTime + 0.62 s` predicts the
  Biopac-locked offset at r = 0.993 on 46 runs with a reliable reference.
- *Global-signal coupling:* on 51 dual-belt runs, PMU RVT⊗RRF couples to the
  global signal at median +0.34, positive in 90%, against Biopac's +0.22 on
  the same runs.
- *Inter-belt agreement is poor:* the two belts' RVT time courses correlate at
  median 0.17 over all 51 runs, and 0.05 on the 25 not used to develop the
  method. An earlier figure of 0.34 came from 27 runs pre-selected for good
  raw-trace alignment, and does not hold.

**Decided (by the PI, 2026-11-25).** The PMU fallback (rule and gates in prereg
§6) supplies RVT for 18 further runs, giving 41 participants and H3 coverage of
21/23. It is used only in robustness analysis **S2b**. The confirmatory RVT
comparison (S2) and the H5 TOST stay on the Biopac belt, because PMU RVT is not
interchangeable with the registered measure. The 7 PMU runs at 10–15%
saturation are not used.

**Disclosure (also in prereg §6).**
- *Saturation gate.* The gate pre-set at < 5% of the whole file refused 17 of
  28 candidates. It was replaced by < 10% in the scan window, derived from the
  dual-belt validation runs only. The fallback runs' quality metrics, including
  their global-signal r, had been printed before the replacement.
- *Alignment-reliability filter.* This (sharpness ≥ 0.05) was set after the
  residuals were seen. It changes no constant or gate.

---

# v0.5 decisions (2026-11-25, pre-freeze)

## Decision 12 — ROI masks from the CAREN atlas; H1's inference set fixed

**Problem.** §8 of v0.4 named six Harvard–Oxford (HO) masks at a 50% threshold.
Two of them cannot be built from that atlas:
- HO has a single "Insular Cortex" label. The registered anterior insula (the
  primary mask for H2/H3/H4) and posterior insula have no HO definition, so an
  analyst would have had to choose a cut-off after the freeze.
- HO has no dACC/mid-cingulate label. Its "Cingulate Gyrus, anterior division"
  runs from pregenual to mid-cingulate cortex (1531 voxels, centroid y = +19).
H1's text said "the insula/operculum and somatomotor ROIs" without naming the
masks or the Bonferroni count.

**Decision.** Masks are defined from the CAREN consensus atlas of resting-state
networks (Doucet, Lee & Frangou, 2019) by network × AAL-3 region membership.
This is the procedure locked for the lab's MindLock study (MindLock
`Docs/Preregistration_Deviations.md`, C1 and C6):
- Salience network ∩ Insula_L and ∩ Insula_R give the left and right anterior
  insula (centroids y ≈ +5 / +7).
- Salience ∩ (Cingulate_Mid, ACC_sup) gives dACC/aMCC (centroid (0, +6, +36)).
  The IEAT coordinate (−4, +12, +36) falls inside it.
- Somatomotor network ∩ Insula gives the posterior insula (centroid y ≈ −15).
  Somatomotor ∩ Rolandic_Oper gives the central operculum, and ∩
  Precentral/Postcentral gives sensorimotor cortex.

**Why CAREN.** CAREN is a consensus of published network atlases rather than
a single parcellation. Its networks split the insula along the
anterior/posterior functional division without any coordinate chosen by us.
It is the same source and procedure as MindLock, so the lab's two studies
share ROI definitions.

**H1 inference set.** The three bilateral somatomotor-network masks are
posterior insula, central operculum and sensorimotor cortex, Bonferroni ×3.
They are the three regions named in the IEAT's pacing result ("insula/operculum
pathway, somatomotor cortex"). The opposite tail is tested in the same masks
at the same correction, because H1 is registered as diagnostic in both
directions.

**Atlas sensitivity.** An HO version is registered for checking (insula split
at y = 0; HO anterior cingulate; central opercular cortex; pre/postcentral)
and is never used for inference.

**Freeze mechanism.** `Scripts/intero2023_build_rois.py` builds every mask on
the fMRIPrep res-2 grid, and the masks are checksummed:
- `ROIs/roi_manifest.csv` holds voxels, centroid, and SHA-256 of each mask
  file and its data.
- `ROIs/roi_sources.csv` holds the SHA-256 of every atlas file used.
- `ROIs/roi_overlap.csv` holds the pairwise shared voxels.
MindLock lost a set of primary masks when one variant overwrote another at a
fixed path, so every mask here has its own file name. Checksums are the
record.

The NIfTI files are gitignored because CAREN's redistribution terms are not
established and the repository is public. Building the masks read the
fMRIPrep brain masks only, never BOLD intensities.

---

## Applied

`Preregistration/intero2023_imaging_prereg.md` was rewritten to **v0.2** on
2026-08-24 to carry all four decisions: §3 (hypotheses and primacy), §4 (sample
and Ns), §5 (contrast definition and trailing events), §7 (the three models),
§8 (power table and the H3 selection statement), §9 (S6), §12 (resolved items
6–9), plus a change log in the header. This document remains the record of
*why*; the prereg body is now the record of *what*.

Decisions 5–9 were applied as **v0.3** on 2026-11-24: §0, §1b, §2, §3 (H1, H2, H4,
H5, H6), §4, §5, §8, §9, §12 (items 10–13), and the header change log. The OSF
draft, `osf_template_answers.md` and `ANALYSIS_HANDOFF.md` were brought in line
in the same pass.

Decisions 10–11 were applied as **v0.4** on 2026-11-25 (prereg header, §3 H5, §4, §6, §8, §9 S1/S2/S2b/S3, §12.4, §12 item 14), with the code and data changes in commits `4aae7cd` and `bd88f7f`.

Decision 12 was applied as **v0.5** on 2026-11-25 (prereg header, §0, §3 H1 and H6, §8, §12 item 15). The OSF draft, `osf_template_answers.md` and `ANALYSIS_HANDOFF.md` were updated in the same pass.

## What is still open

- Optional and manifest-neutral: regenerating the 95 confound sheets so `csf`
  and `white_matter` come from the derivatives of record rather than the earlier
  anatomical pass (`VERIFICATION_REPORT.md` §2a). Staged at
  `analysis/verify_regen/`, not adopted.

**The freeze gate remains in force.** No first-level model is fitted until the
OSF registration is frozen (`ANALYSIS_HANDOFF.md` §0).
