# Draft correction to prereg §4 — sample, exclusions, analysable N

**Status: APPLIED to `intero2023_imaging_prereg.md` v0.2 on 2026-08-24**, while
the registration was still unfrozen. Retained as the evidence record for that
amendment. Two things changed on application: the sensitivity analysis proposed
below as "S3" is registered as **S6** (S3 was already taken by the
CO2-availability cohort check), and the corrected Ns are one to two higher than
the table below because the two-model split (`methods_decisions.md` decision 2)
was adopted at the same time — H1 47, H2 44, H3 23, H4 43, H6 46. The
tables below are left at their pre-split values as the record of what was
proposed; §4 of the prereg carries the registered figures. Evidence:
`VERIFICATION_REPORT.md`.
This is a **pre-freeze correction**, not a deviation: no model has been fitted
and no BOLD image has been read. If the registration is already frozen when
this is read, everything here becomes a `deviations.md` entry instead.

---

## What is wrong with §4 as written

**(a) The a-priori exclusion of sub-082 and sub-104 describes data that no
longer exists.** §4 reads:

> **sub-082 and sub-104** — run labels are shifted in the derivatives tree (the
> file named `task-breath1` holds 172 volumes where the breathing slot expects
> 242; 172 is the authenticity length).

That is true of the Windows `Data/Ready` snapshot, which mapped scanner runs
positionally and so labelled the 3rd acquired run of a 3-run session as
`breath1`. The BIDS tree on the analysis machine maps by scanner slot number
and is correct:

```
sub-082:  Run1 -> task-localizer_run-1  (196 vols)
          Run3 -> task-breathing_run-1  (242 vols)
          Run5 -> task-authenticity     (172 vols)
sub-104:  identical structure
```

`Data/qc_ledger.csv` already agrees — it grades both breathing runs
`242 vs v2 mode 242 (0.0%)`, *exact*. The registration therefore excludes two
participants on the basis of a labelling defect its own QC ledger records as
absent. Both runs pass the registered motion gate (sub-082 fd_mean 0.241,
sub-104 0.154) and both events files fit inside the acquisition.

**(b) Two participants were never screened at all.** §4's "+ usable imaging =
95 runs / 48 participants" line silently drops sub-065 and sub-076, whose
breathing runs were absent from the imaging inventory used to build the
manifest. That inventory predates the 2026-08-06 re-conversion of the dataset;
it still carries sub-065 `task-localizer_run-1` at 167 volumes, which is the
truncated file that re-conversion replaced with the full 196. Both participants
have complete 5-run sessions and pass the registered gate on both breathing runs.

**(c) sub-056 is excluded correctly but silently.** It is inside the "runs with
no usable event timing" clause without being named or evidenced.

**(d) The registered Ns of 43 for H1 and H2 are not attainable** under the
model as §7 defines it — see the separate contrast-definition item below and
`VERIFICATION_REPORT.md` §4c–5. This is not a sampling problem and is not
fixed by (a)–(c).

---

## Proposed replacement text for the exclusion list

> **Excluded a priori:**
>
> - Pilot/QA identifiers P00, P01, P50, P999
> - **sub-051** — behavioural session only; no MRI acquired.
> - **sub-056** — `task-breathing_run-1` is a truncated acquisition (229
>   volumes; the event record runs 28.3 s past the end of the scan), and
>   `task-breathing_run-2` fails the registered motion gate (mean FD 0.471 mm,
>   26.7 % of volumes over 0.5 mm). No usable breathing run remains.
> - Runs with no usable event timing or no trigger-locked physiology
>
> **Corrected 2026-08-20, before the registration was frozen.** The imaging
> inventory used to build the original run list predated the 2026-08-06
> re-conversion of the dataset and was wrong in two ways, both verified against
> the fMRIPrep 24.1.1 derivatives of record:
>
> - **sub-082 and sub-104 are reinstated.** Their run labels are not shifted in
>   the dataset of record; each contributes one gate-passing breathing run.
> - **sub-065 and sub-076 are added.** Both have complete two-run breathing
>   sessions (242 volumes each) that the stale inventory omitted; both pass the
>   registered motion gate on both runs. sub-065 has no MAIA questionnaire and
>   therefore does not enter H6.
>
> No imaging result had been examined when this correction was made; the
> criteria applied are the ones already registered, applied to a corrected
> inventory.

## Proposed replacement for the analysable-sample table

| | runs | participants |
|---|---|---|
| Breathing runs with events | 105 | 54 |
| + usable imaging | **103** | **53** |
| − sub-056 (truncated run-1, gated run-2) | 101 | 52 |
| **− motion gate (11 of the remaining 101 runs)** | **90** | **47** |
| + RVT regressor | 64 | 30 (all runs) |
| + PetCO2 regressor | 51 | 24 (all runs) |
| + both physiological regressors | 43 | 20 (all runs) |

Across the full 103-run set the gate fails 12 runs; one of those is sub-056
run-2, already removed on the line above. It still removes exactly the same 5
participants entirely (sub-058, -062, -080, -085, -094 — unchanged from the
registered version). **Post-gate: 90 runs, 47 participants.** None of the four
reinstated participants has breathing-run RVT
or PetCO2, so the physiological rows and the H5 TOST bounds are unchanged.

## Proposed replacement for the final registered Ns

Counting a participant when the contrast is estimable in at least one
gate-passing run. Both columns assume the §7 contrast correction is adopted;
if it is not, use the "as written" figures from `VERIFICATION_REPORT.md` §5.

| Contrast | registered | **corrected** |
|---|---|---|
| C1 paced vs free (H1) | 43 | **46** |
| C3 confidence parametric (H2) | 43 | **42** |
| C2 aware vs unaware (H3) | 21 | **23** |
| C5 salience (H4) | 41 | **43** |
| H6 MAIA moderation | 43 | **45** |
| S1 ±PetCO2 | 22 | 22 (unchanged) |
| S2 ±RVT | 26 | 26 (unchanged) |

H2 lands at 42 rather than 46 because three participants (sub-066, sub-070,
sub-114) report a single confidence value across every retained run, leaving no
modulator to estimate, and one further participant's runs are rank-deficient.
That is a property of the behaviour, not of the sample correction, and it
should be stated in §7 rather than buried in an N.

**H6 power** updates from 80 % at r ≥ 0.374 (N = 43) to 80 % at **r ≥ 0.366**
(N = 45), one-sided. Power at the IEAT's r ≈ 0.32 rises from 67 % to 69 %.
The correction does not change the study's power situation materially and
should not be described as if it does.

---

## Companion corrections needed elsewhere in the prereg

These are separate amendments. They are listed here because the N table above
depends on the first of them.

1. **§7, `paced_gt_free`.** Defined as a fixed 3-way average
   `(aware + unaware + nochange)/3 − free`, it is undefined in 18 of 91 runs —
   17 lack `paced_unaware` because the participant never answered "No Change"
   on a change trial. Proposed: average over the paced conditions **present in
   that run**. Precision is unchanged (median SE multiplier 0.0830 vs 0.0840; corrected: 0.247 for both, decision 18).
2. **§7, confidence modulator.** Median |r| with the aware-minus-unaware
   direction is 0.437 across runs, exceeds 0.5 in 26 of 70, and is exactly 1 in
   9 runs, where the design becomes rank-deficient. H2 and H3 are not fully
   separable in a substantial minority of runs; the registration should say so,
   and should state what happens to those 9 runs and to the 3 participants with
   no confidence variance.
3. **§5, trailing events.** Two retained runs carry `free`/`rating` events that
   over-run the acquisition (sub-053 run-1 by 25.6 s, sub-052 run-1 by 11.0 s)
   because the scanner stopped mid-rating. Proposed: clip to the scan end and
   log; do not let the modelling library discard them silently.
4. **Confound provenance.** The committed `csf` and `white_matter` columns do
   not reproduce from the fMRIPrep derivatives of record (median r 0.985 and
   0.996). Every other column reproduces bit-exactly. Regenerate the sheets
   before fitting so the analysis matches its stated inputs.

---

## Decision rule for the reinstated participants, and its sensitivity analysis

**Added 2026-08-21, before the freeze and before any model was fitted.**

### Why this decision can be made now

Inclusion here is settled on evidence that is **blind to the effect being
tested**. Three tiers are distinguishable, and only the third is off-limits
before a registration is frozen:

1. **Design and timing** — trial counts, contrast estimability, design
   efficiency. No BOLD voxel is read.
2. **Image quality that is not a function of the tested contrast** — tSNR,
   framewise displacement, susceptibility-distortion correction, ghosting.
   §4 already relies on this tier: the registered motion gate is computed from
   fMRIPrep's `framewise_displacement`, which is BOLD-derived. Admitting FD
   while refusing tSNR would be arbitrary.
3. **Quantities computed from the contrast estimates** — whether including a
   participant moves a result. This is the only prohibited tier, and nothing
   below requires it.

### The rule

> **Single-run participants are retained, and the decision is not specific to
> any participant.** A participant enters the analysable sample when the
> registered criteria — usable event timing, trigger-locked imaging, and the
> motion gate — are satisfied by at least one breathing run, irrespective of
> how many runs they contribute. Where a participant contributes one run, the
> fixed-effects step is that run.

Stated this way the rule cannot be gerrymandered after the fact: it applies
identically to participants **already inside the registered sample**. Two of
the original 43 — **sub-074 and sub-114** — contribute exactly one breathing
run, and sub-114 additionally yields no estimable H1 or H2 contrast under the
§7 definition. Excluding sub-082 and sub-104 for contributing one run each
would therefore not apply a registered standard; it would apply a standard to
those two participants alone.

### The result-blind evidence, for the record

*(Corrected 2026-11-25, `methods_decisions.md` decision 18: the SE multipliers in this document came from a scaling error in the QC script.)* Corrected values for the adaptive C1, on the current sheets:
registered sample median 0.175, worst 0.263; sub-082 0.242 and sub-104 0.250, each with 1 of 43 registered participants less precise.

Subject-level SE multiplier for C1 as originally computed (design algebra only; lower is more precise;
registered sample median 0.0584, worst 0.1144):

| participant | SE multiplier | position in the registered sample |
|---|---|---|
| sub-082 | 0.0836 | 8 of 42 registered participants are less precise |
| sub-104 | 0.0717 | 10 of 42 registered participants are less precise |

Image quality (MRIQC, acquired June 2026, blind to any task contrast):
sub-082 tSNR 30.2 (68th percentile of breathing runs), sub-104 tSNR 29.1 (46th
percentile) — both at or above the cohort median. Both runs have a PEPOLAR
fieldmap and susceptibility-distortion correction applied. sub-082 carries more
low-amplitude drift (mean FD 0.252 mm) but only 0.4 % of its volumes exceed the
registered 0.5 mm threshold.

Neither participant is an outlier on any result-blind index.

### Registered sensitivity analysis S6 — two-run participants only

*(Proposed here as S3; registered as **S6** to avoid colliding with the existing S3.)*

Every hypothesis is additionally tested on the subset contributing **two**
breathing runs, which removes **sub-074, sub-114, sub-082 and sub-104**
together. The subset is defined by run count alone, so the two reinstated
participants cannot be removed without also removing two participants from the
original registered sample.

| | participants | H1 | H2 | H3 | H6 |
|---|---|---|---|---|---|
| **Primary** (corrected sample) | 47 | 46 | 42 | 23 | 45 |
| **S6** (two-run participants only) | 43 | 42 | 39 | 23 | 41 |

Interpretation is fixed in advance: agreement between the primary and S3
results is reported as such; divergence is reported as a limitation of the
single-run participants' precision and **not** used to select which analysis to
present.

### What this does and does not buy

Honestly stated: **the power gain is small.** H6 moves from 80 % power at
r ≥ 0.374 (registered N = 43) to r ≥ 0.362 (N = 46 under the adopted models); S6 sits at r ≥ 0.378.
The case for reinstating sub-082 and sub-104 is that the inventory used to
exclude them was wrong and the exclusion is inconsistent with participants
already in the sample. It is a correctness and transparency argument, not a
statistical rescue, and should not be presented as one.
