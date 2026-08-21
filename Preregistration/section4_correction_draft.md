# Draft correction to prereg §4 — sample, exclusions, analysable N

**Status: DRAFT for the registration holder. Not applied to
`intero2023_imaging_prereg.md`.** Evidence: `VERIFICATION_REPORT.md`.
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
   that run**. Precision is unchanged (median SE multiplier 0.0830 vs 0.0840).
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
