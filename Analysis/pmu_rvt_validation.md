# PMU respiration as a fallback RVT source — validation summary

Generated 2026-11-24 by `Scripts/intero2023_pmu_rvt.R`, before the registration
freeze and before any first-level model was fitted. Per-run data:
`Analysis/pmu_rvt_alignment.csv`, `Analysis/pmu_rvt_validation.csv`,
`Data/rvt_pmu_runs.csv`. No BOLD image was read. fMRIPrep's `global_signal`
column (the disclosed exception) was used only for the validity check below,
never for alignment or selection of any run.

## Rule

- **Source priority per run:** Biopac belt if `rvt_rrf_z` covers ≥ 50% of
  volumes; otherwise the Siemens PMU cushion if it passes every gate; otherwise none.
- **Alignment:** scan start in the PMU log = `AcquisitionTime − LogStartMDHTime + 0.62 s`.
  Both times are on the scanner's MDH clock. Only the derived offset is stored,
  never a time of day.
- **Gates:**
  - the log covers ≥ 95% of the scan
  - band ratio ≥ 0.90
  - < 10% of raw samples on a 12-bit rail within the scan window
  - ≥ 5 kept breaths
  - RVT defined on ≥ 50% of volumes
  - polarity inverted
- **Processing:** the same `belt_analyse → rvt_series / rvt_convolved` chain as the
  Biopac route.

## Alignment

Measured against the Biopac-trigger-locked offsets in `Data/pmu_wide_offsets.csv`:

| runs | prediction r | residual median | IQR | range |
|---|---|---|---|---|
| all 50 | 0.157 | −0.53 s | −0.82 to 1.59 | −259.9 to +18.8 s |
| reference reliable (sharpness ≥ 0.05), 46 | **0.993** | −0.60 s | −0.82 to 1.55 | −1.13 to +1.98 s |
| decision set, 26 (2026-11-24) | 0.997 | −0.62 s | −0.77 to −0.27 | −1.03 to +1.83 s |

The four runs excluded from the reliable set are ones whose *reference* is not a measurement:
- **sub-070 run-1:** its Biopac offset was computed from the wrong recording, because of the
  cache-key bug fixed in commit 4aae7cd. Sharpness 0.029, residual −260 s.
- **sub-089 runs 1–2:** the Biopac–PMU cross-correlation has no peak (sharpness
  0.041 and 0.003, |r| ≈ 0.2). Residuals +7 s and +19 s.
- **sub-085 run-1:** sharpness 0.013, residual +1.8 s.

The sharpness filter was set after these residuals were seen; it changes neither the
constant nor any gate. Residuals in the reliable set are bimodal, about −0.7 and
+1.6 s, roughly one TR apart. That is consistent with TR-level ambiguity in the Biopac
reference, which was searched on a 1 s grid. Every residual is within one TR. Over the
46 reliable runs the calibration implied is +0.60 s; the fixed constant is +0.62 s.

## Dual-belt agreement

Runs with both a Biopac RVT file and a PMU file:

| set | n | PMU–Biopac RVT⊗RRF r, median (IQR) | rate r | RVT⊗RRF ~ global signal, PMU | same, Biopac |
|---|---|---|---|---|---|
| decision set (judged on, 2026-11-24) | 27 | 0.34 (0.10–0.49) | 0.44 | +0.34 | +0.23 |
| further breathing runs | 25 | **0.05 (−0.15–0.30)** | 0.33 | +0.34 (92% > 0) | +0.22 (80% > 0) |
| all breathing runs | 51 | 0.17 (−0.04–0.38) | 0.41 | +0.34 (90% > 0) | +0.22 (82% > 0) |

**The decision set over-states RVT agreement.** It was chosen for sharp raw
Biopac–PMU cross-correlation, i.e. runs where the two traces already agreed. On the
25 runs added afterwards, the two RVT time courses barely agree (median r 0.05).
Breathing rate agrees better across the board (median 0.41).

What holds in every subset is the validity check the preregistration applies to RVT:
the PMU RVT⊗RRF predicts the global signal positively, as strongly as or more strongly
than the Biopac RVT does on the same runs (median +0.34 vs +0.22). The two transducers
sit on different parts of the trunk: the Biopac strain belt and the Siemens pneumatic
cushion. This analysis cannot say which RVT is closer to ventilation.

By window saturation, over all 51 dual-belt breathing runs:

| saturation | n | RVT agreement, median | GS r PMU | GS r Biopac |
|---|---|---|---|---|
| < 5% | 8 | 0.48 | 0.41 | 0.13 |
| 5–10% | 19 | 0.15 | 0.41 | 0.22 |
| 10–15% | 11 | 0.21 | 0.30 | 0.23 |
| ≥ 15% | 13 | 0.09 | 0.15 | 0.28 |

## Fallback result

29 of the 90 analysable runs lack a usable Biopac RVT. Of these:
- **18 runs (11 participants) pass and receive PMU RVT.** Window saturation is 0.4–8.6%;
  RVT⊗RRF vs global signal is +0.07 to +0.67, median +0.40, positive in 18/18.
- **10 runs are refused** for saturation ≥ 10%: sub-052 r1/r2, sub-053 r1/r2, sub-060 r1/r2,
  sub-075 r1/r2 and sub-120 r1/r2.
- **1 run (sub-054 r1) has no PMU file.**

The 18 runs:
- sub-055 r2
- sub-057 r1, r2
- sub-071 r1, r2
- sub-076 r2
- sub-081 r1, r2
- sub-087 r1
- sub-102 r1, r2
- sub-118 r1, r2
- sub-132 r1, r2
- sub-134 r1
- sub-136 r1, r2

## Provenance of the saturation gate (for the registration's disclosure)

1. **Pre-set gate.** Before any fallback run was scored, a gate of < 5% saturation over the
   whole PMU file was set.
2. **First scoring.** Scoring printed every fallback run's metrics, including its RVT⊗RRF vs
   global-signal r. The 5% gate refused 17 of 28 runs on saturation alone.
3. **Re-derivation.** The 26 dual-belt validation runs themselves reached 27% saturation. The
   gate was therefore re-derived from those validation runs alone, using PMU–Biopac RVT
   agreement by saturation bin (median r 0.63 / 0.48 / 0.22 / 0.34 for < 5 / 5–10 / 10–15 /
   ≥ 15%). It was set at < 10%, measured in the scan window rather than the whole file.
4. **What informed the cut.** The cut was not taken from the fallback runs' metrics, but those
   metrics had been seen when it was chosen.
5. **Later check.** Extending validation to all 51 dual-belt runs does not reproduce the
   5–10% agreement figure (0.15, against 0.48 in the decision set). The < 10% cut therefore
   rests on (a) the < 5% / ≥ 15% contrast in agreement and (b) the global-signal coupling,
   which stays ≥ +0.41 below 10% and falls to +0.15 at ≥ 15%. It does not rest on a clean
   agreement gradient.
6. **Runs at 10–15%** (7 runs; 10 refused in total with the ≥ 15% runs) are not used.
