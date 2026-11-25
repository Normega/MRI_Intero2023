# Registered methods decisions

**Status:** decided 2026-08-24, pre-freeze. Nothing here has been estimated on
BOLD data. Every number below comes from design matrices, event sheets and
behavioural counts only — `Scripts/intero2023_designqc.py` and
`Scripts/intero2023_twomodel_qc.py`, outputs in `Analysis/design_qc.csv` and
`Analysis/design_qc_twomodel.csv`.

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
`sqrt(c'(X'X)^-1 c)` is **0.235** for both the adaptive and the fixed form, in
the 72 runs where both are defined. *(Corrected 2026-11-25, decision 18. This
entry originally quoted 0.0827 and 0.0830, from a scaling error in the QC
script; the conclusion is unchanged.)*

**What it buys:** H1 goes from **N = 45** to **N = 47**, and the two
participants recovered (sub-070, sub-114) are recovered because a cell was
empty, not because their data improved.

**Stated as a limitation:** the paced side is a slightly different weighted
average in different participants. Given that the weights carry no theoretical
content, this is a smaller distortion than dropping 18 runs. Runs where only
one paced condition is present contribute a single-condition contrast; the
per-run composition is logged in `Analysis/design_qc.csv`
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
`Analysis/design_qc_twomodel.csv` so the subsample cannot drift:

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
participant-level H1 SE multiplier is 0.242 and sub-104's 0.250. Both are inside
the range of the registered sample (median 0.175, worst 0.263), but near its
imprecise end: 1 of the 43 registered participants is less precise than each.
That is expected for one-run participants, and S6 removes them together with
the other two. *(Corrected 2026-11-25, decision 18. This entry originally quoted
0.0836 and 0.0717, with 8 and 10 registered participants less precise, from
the QC scaling error. The reinstatement rests on the inventory being wrong, not
on these figures.)* tSNR 30.2 and 29.1 (68th and 46th percentile); both have PEPOLAR
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

*Superseded v0.8 (decision 19):* at H3 N = 45 the subsets are 25 and 28, and
the confirmatory TOST moves to H3.

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

The NIfTI files are not published, in the repository or on OSF. This was
decided by the PI on 2026-11-25: the masks derive from a third-party atlas.
The generating script and the checksums are published. Building the masks read the
fMRIPrep brain masks only, never BOLD intensities.

## Decision 13 — H7: between-network coupling, paced vs free and aware vs unaware

**Request.** Test for greater inter-network coupling in the breathing-attention
and breathing-awareness contexts, using the MindLock connectivity analysis as
the template. The PI chose:
- **Contexts:** paced vs free as the attention contrast (H7a), and aware vs
  unaware as the awareness contrast (H7b).
- **Pairs:** SAL–SMN and SAL–ECN.
- **Status:** both parts confirmatory.

**What is carried over from MindLock:**
- whole CAREN network masks;
- median network signal, with the mean as a sensitivity specification;
- a 0.01–0.1 Hz band-pass applied identically to data and regressors;
- per-epoch demeaning before concatenation;
- Fisher z, with a paired test equivalent to `z ~ condition + (1|participant)`;
- no global signal regression in the primary pipeline;
- the non-neural-account checks, reported with the result;
- an exploratory all-pairs specificity analysis with FDR;
- the rule that the per-participant coupling difference is not an
  individual-differences measure.

**What changes, and why:**
1. **Task-evoked responses are regressed out before correlating** (Model A
   regressors). MindLock's 50 s blocks made this optional. Here the epochs are
   about 20 s and abut rating windows, so the evoked responses at every switch
   are shared across networks and would inflate every correlation (Cole et
   al., 2019).
2. **The global-signal check is a pre-committed interpretation gate rather
   than a post-hoc test.** Paced breathing changes ventilation and CO2 by
   construction. H1 also predicts lower regional activity under pacing, which
   is exactly the condition MindLock identified as pushing measured coupling
   up for non-neural reasons. A supported H7 is called a network-coupling
   effect only if it survives partialling the identically processed global
   signal.
3. **Spike volumes are excluded from the windows** instead of being regressed,
   because spike regressors do not survive band-pass filtering.

**Feasibility, counted from event timing before the freeze** (no BOLD read):
- Epochs are about 20 s at TR = 2 s, with a 4 s window shift.
- Every participant has 7–14 paced and 7–14 free epochs, a median of 139 and
  141 volumes (minimum 70).
- In the 23 H7b participants (the v0.5 H3 set), the medians are 67 aware and 44 unaware volumes
  (ranges 34–89 and 27–77).
- MindLock had about 500 volumes per mode, at a shorter TR.

**Power at 80%:**
- **H7a (N = 47):** dz 0.42 at α .05 and 0.46 at Holm's first step.
- **H7b (N = 23):** dz 0.61 and 0.68.
- **Caveat:** MindLock's registered coupling effect was dz 0.49, measured on
  more data. H7b is registered as confirmatory knowing that it detects only
  large effects, and its null is pre-declared uninformative.

**Prior evidence, stated so that it is not overread.** MindLock's SMN–SAL
test was null (dz 0.17). Its coupling changes were confined to ECN pairs.

---

# v0.6 decisions (2026-11-25, pre-freeze)

## Decision 14 — H3 enters every participant in whom the contrast exists, precision-weighted

**Request.** The PI asked why so few participants qualify for the aware vs
unaware contrast, and chose to include everyone with a precision weighting.

**Why the v0.5 rule bit.**
- The design offers at most 12 change trials per participant (6 per run).
- Participants answer "No Change" on only 28% of them. That is the liberal
  bias in the behaviour: P(say change) is about 0.70.
- Of the 24 participants the ≥ 3-trial rule excluded, 23 failed on the unaware
  cell alone, and 46 of 47 had ≥ 3 aware trials.
- The rule selected toward less frequent detectors. The included 23 had a
  median unaware share of aware + unaware trials of 0.36, against 0.19 for
  those excluded.

**Decision.**
- **Inclusion:** a participant enters H3 when at least one gate-passing run
  contains ≥ 1 aware and ≥ 1 unaware trial. This gives N = 45. sub-070 and
  sub-114 have no unaware trial and cannot enter.
- **Participant contrast:** *(superseded by decision 17: a joint two-run
  model)* the fixed-effects (inverse-variance) combination over the runs that
  contain both cells. In 28 of 45 participants both runs
  contribute. Aware trials in a run with no unaware trial are lost to the
  contrast: 56 of 290 (19%), in 14 participants. Every unaware trial is used.
- **Group model:** FSL FLAME 1 on each participant's effect, variance and
  degrees of freedom. ROI inference is voxel FWE by GRF within each mask,
  Bonferroni ×3, and whole-brain inference is GRF cluster FWE. Permutation is
  not used for H3, because a precision-weighted model is not exchangeable
  under sign-flipping. As a robustness check, the ROI means go into a REML
  random-effects meta-analysis.
- **S7:** the v0.5 rule and its 23 participants, fitted with FLAME and with
  the v0.5 OLS model.
- **H7b keeps the ≥ 3-epoch rule (N = 23).** A correlation over one or two
  epochs of about 10 volumes is not estimable in a useful sense.
- The sets are frozen in `Data/h3_sets.csv`
  (`Scripts/intero2023_h3_sets.R`). The S7 set reproduces the registered
  v0.5 list exactly.

**Why FLAME rather than a lower threshold.** At ≥ 2 trials per cell, N would
be 36, and the selection would remain (unaware share 0.33 vs 0.08). Precision
weighting is the standard treatment of unequal trial counts. It keeps
low-count participants at low weight instead of discarding them, and removes
most of the selection. The cost is a second group-level engine, used for H3
only. The handoff's no-mixing rule is amended for this one model.

**Consequences.**
- H3 Ns in the physiological subsets: S1 25, S2 28, S2b 39, S4 26, S6 42.
- The H5 comparison on H3 stays descriptive, although its subsets now match
  H2's. The confirmatory target was fixed on H2 before this change, and H3's
  per-participant estimates are less precise.
- H2's ruled-out inference ("a null H3 beside a positive H2 …") still
  holds. The Ns now match (45 vs 44), but the comparison is still confounded
  with per-participant precision.

**What was seen.** Trial counts from the event sheets only. No imaging
data.

## Decision 15 — Tissue confounds regenerated from the derivatives of record

**Request.** The PI confirmed that CSF should be controlled and that the
open regeneration item should be adopted.

**Decision.** `csf` and `white_matter` in all first-level confound sheets now
come from the fMRIPrep derivatives of record
(`Scripts/intero2023_tissue_confounds.R`). The script was written fresh,
rather than copying `analysis/verify_regen/`, because those staged sheets
predate the v0.4 physio fixes.
- **What changed:** 95 sheets (the registered participants). The 8 sheets of
  the reinstated participants were already built from the record.
- **Agreement with the old columns:** median r 0.985 / 0.996, with minima
  0.887 / 0.936 (`Analysis/tissue_confounds_regen.csv`).
- **Checks:**
  - Every rewritten value equals the staged regeneration.
  - No other column changed in any sheet.
  - Each sheet's motion and FD columns were matched to the fMRIPrep file
    before rewriting.
- **Unchanged:** the manifest, gates and Ns. The nuisance set itself is
  unchanged: both columns were always in it.

## Decision 16 — Exploratory graded-awareness model

**Request.** The PI asked to explore treating a low-confidence aware response
as less aware.

**Decision.** Model A-g (exploratory, prereg §10):
- **Regressor:** one regressor for aware and unaware change epochs together.
- **Modulator:** g = +confidence on aware trials and −confidence on unaware
  trials (confidence on a 1–7 scale), mean-centred within run.
- **Other conditions:** as in Model A.
- **Test:** the H3 masks, with the H3 group model, N = 45.
- **Variants:** aware-only confidence, and rank-transformed g.

**Behavioural background.** Confidence is high in both cells. Aware trials
are mostly rated 5–7, unaware trials mostly 4–6, and few trials are rated
below 4. The low-confidence end of g is therefore sparse.

**Why it is exploratory.** g is a monotone function of the H3 category plus
within-category confidence, so it overlaps H3 by construction. It also
overlaps H2's modulator, which already absorbs the awareness effect. It
can describe whether confidence grades the awareness contrast. It cannot
adjudicate H3.

# v0.7 decisions (2026-11-25, pre-freeze)

## Decision 17 — H3's first level is one joint model of both runs

**Request.** The PI chose the combined-run approach, so that every trial
contributes to H3.

**The alternative considered and rejected.** One option was to estimate aware
and unaware activity separately, each against the implicit baseline, and
subtract them at the group level. The events tile about 94% of every run, so
a single-condition estimate is measured against a few seconds of gaps, and
baseline differences between runs do not cancel. By design algebra its
variance is a median 4.4 times that of the within-run contrast (78% less total
precision). The contrast is therefore kept at the first level.

**Decision.**
- **Design:** one nilearn GLM per participant over all gate-passing runs:
  - each run is scaled to percent signal change, then the runs are
    concatenated;
  - the Model A task regressors and their derivatives are shared across runs;
  - the base nuisance regressors and the intercept are run-specific;
  - the noise model is AR(1).
- **Group level:** the aware − unaware effect, variance and dof enter FLAME 1
  as before.
- **Reference implementation:** `build_joint_design()` in
  `Scripts/intero2023_h3_joint_qc.py`.
- **What it buys:** total H3 precision +12.5%. In the 14 participants whose
  aware trials a per-run model discards (56 of 290), the contrast variance
  falls by a median 22%. In the other 31 it falls by a median 7%, because
  free, rating and sub-rating are estimated from both runs.
- **What it assumes:** each condition's response and the noise variance are
  the same in both runs. The runs are behaviourally exchangeable (§8).
  Concatenation alone would not make them comparable; the
  percent-signal-change scaling per run does.
- **Robustness:** the v0.6 per-run form (fixed effects over runs with both
  cells) is reported alongside.
- **Model A-g** is fitted in the same joint form.
- **N is unchanged at 45.** The entry rule now reads "≥ 1 of each cell across
  runs", and in these data that is the same set (`Data/h3_sets.csv`, script
  v1.1, which asserts it).

**What was seen.** Design matrices and trial counts only.

## Decision 18 — Design-QC standard errors corrected

**The bug.** `Scripts/intero2023_designqc.py` and
`Scripts/intero2023_twomodel_qc.py` standardise the design columns for
numerical stability. Because beta_std = sd × beta, a contrast c on the raw
betas is c / sd on the standardised ones. Both scripts used c × sd. SEs were
therefore wrong whenever the contrasted columns differed in scale, by about a
factor of 3 in practice. The same vector was used in the estimability check
of rank-deficient designs.

**Second fix.** Both scripts built designs from every sheet column, which
includes the physiological regressors. Those belong to the sensitivity models
only, so the QC now drops them.

**Re-run, 2026-11-25:**
- **Unchanged:**
  - The single-model QC still finds 9 / 91 rank-deficient designs, with the
    confidence-modulator / awareness median |r| at 0.437.
  - The two-model QC still finds 0 / 90 rank-deficient designs in every
    model.
  - Every registered N is reproduced.
- **Corrected figures:**
  - C1 median SE 0.235 for both forms (was quoted as 0.0827 / 0.0830).
  - sub-082 0.242 and sub-104 0.250, against a registered median of 0.175;
    1 of 43 registered participants is less precise than each (was quoted as
    8 and 10). The entries of decisions 1 and 4 are corrected in place.
- **No registered decision rested on the absolute values.** The C1 claim was
  comparative and holds. The reinstatement rested on the inventory.
- The per-participant estimates in `Analysis/h3_joint_design_qc.csv` were
  computed with the corrected formula from the start.

## Decision 19 — The confirmatory H5 equivalence test moves to H3

**Evidence.** Decision 8 kept H3's physiological comparison descriptive because
only 12–15 of v0.5's 23 H3 participants had the regressors. From v0.6 H3 has
45 participants, and the subsets are 25 (PetCO2) and 28 (Biopac RVT), as
large as H2's 25 and 27. Paired TOST power by simulation (α = .05 per side,
true difference 0, 200,000 draws): N = 25 at ±0.61 gives 0.81; N = 28 at
±0.58 gives 0.82 (±0.57 gives 0.807); N = 39 at ±0.49 gives 0.83.

**Why H3 and not H2.**
- H3 is the authoritative awareness test. Without this change, any claim that
  the awareness effect is not respiratory or CO2-driven would be unregistered.
- H2's confidence effect bundles awareness with confidence. Equivalence on H2
  is indirect evidence about H3.
- The confound is sharpest for H3. Detection may track trials on which
  breathing followed the pacer more closely, so aware and unaware trials may
  differ in physiology systematically. H2's confidence weighting is spread
  over all paced trials.

**Costs, registered.**
- H3's per-participant estimates are noisier, so "indeterminate" is more
  likely for H3. Power in dz units is unchanged; the bound is wider in raw
  units, and the corrected–uncorrected correlation is reported so it can be
  read.
- The H3 test only means something if H3 is supported, so the conditional
  rule is stated in advance.
- Two equivalence tests could disagree. H3 governs and H2 is supportive, so
  the verdict comes from one test, and no false-positive inflation arises.

**Decided.**
- *H3 supported:* the H5 verdict is the H3 TOST grid cell, read in the ROI(s)
  where H3 was supported. Bounds: ±0.61 (PetCO2, N = 25) and ±0.58 (Biopac
  RVT, N = 28); S2b ±0.49 (N = 39). The H2 TOST is supportive, and H3 governs
  if they disagree.
- *H3 not supported, H2 supported:* the verdict is the H2 TOST cell, stated as
  a conclusion about the confidence measure only.
- *Neither supported:* both are descriptive.
- *Weighting:* the paired difference is equal-weighted, although H3's group
  model is FLAME-weighted. Weighting it would need the covariance of two
  estimates from nearly identical models, which FLAME does not supply. Equal
  weighting is the conservative choice.
- *Corrected-effect significance (grid cell 2):* the same one-sided ROI test
  as the hypothesis, on the corrected model in the same subset (FLAME 1 for
  H3).

**Also changed.** The design-QC scripts now write into the repository
(`Analysis/design_qc.csv`, `Analysis/design_qc_twomodel.csv`). The copies
committed before were from August, before the decision-18 SE fix. A re-run
reproduces the verified outputs byte for byte.

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

Decisions 12–13 were applied as **v0.5** on 2026-11-25 (prereg header, §0, §1, §3 H1, H6 and H7, §7 connectivity model, §8, §10, §12 items 15–16). The OSF draft, `osf_template_answers.md` and `ANALYSIS_HANDOFF.md` were updated in the same pass.

Decisions 17–18 were applied as **v0.7** on 2026-11-25 (prereg header, §3 H3, §4, §5, §8, §9 S7, §10 Model A-g, §12 items 19–20). `Data/h3_sets.csv` was rebuilt (script v1.1), and the QC outputs were re-run (written to `Analysis/` from v0.8).

Decisions 14–16 were applied as **v0.6** on 2026-11-25 (prereg header, §3 H2, H3, H5 and H7b, §4, §5, §7, §8, §9 S1–S7, §10, §12 items 17–18). The data changes are in `FirstLevel/*_desc-model_confounds.tsv` (tissue columns) and `Data/h3_sets.csv`. The OSF draft, `osf_template_answers.md` and `ANALYSIS_HANDOFF.md` were updated in the same pass.

Decision 19 was applied as **v0.8** on 2026-11-25 (prereg header, §3 H5, §4 N table, §8 design QC, §9 S1/S2/S2b, §12.4, §12 items 12 and 21). The OSF draft, `osf_template_answers.md` and `ANALYSIS_HANDOFF.md` were updated in the same pass.

## What is still open

- Nothing in the analysis plan. Three legacy tables in `Data/` hold stale
  rows and are read by no registered analysis (`Data/STALE_TABLES.md`).

**The freeze gate remains in force.** No first-level model is fitted until the
OSF registration is frozen (`ANALYSIS_HANDOFF.md` §0).
