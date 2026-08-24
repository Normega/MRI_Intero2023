# Registered methods decisions

**Status:** decided 2026-08-24, pre-freeze. Nothing here has been estimated on
BOLD data. Every number below comes from design matrices, event sheets and
behavioural counts only — `Scripts/intero2023_designqc.py` and
`Scripts/intero2023_twomodel_qc.py`, outputs in `analysis/design_qc.csv` and
`analysis/design_qc_twomodel.csv`.

This document records four choices that `VERIFICATION_REPORT.md` §7 identified
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
preprocessed data:

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
together with the single-run decision rule and sensitivity analysis **S3**
already drafted there. sub-056 stays excluded — its run-2 fails the registered
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
`section4_correction_draft.md`, including sensitivity analysis S3.

---

## What is still open

- Applying these edits into `Preregistration/intero2023_imaging_prereg.md`
  itself. This document and `section4_correction_draft.md` hold the replacement
  text; the prereg body has not yet been rewritten.
- Optional and manifest-neutral: regenerating the 95 confound sheets so `csf`
  and `white_matter` come from the derivatives of record rather than the earlier
  anatomical pass (`VERIFICATION_REPORT.md` §2a). Staged at
  `analysis/verify_regen/`, not adopted.

**The freeze gate remains in force.** No first-level model is fitted until the
OSF registration is frozen (`ANALYSIS_HANDOFF.md` §0).
