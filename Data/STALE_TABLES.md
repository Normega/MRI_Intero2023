# Legacy tables with known stale rows

Written 2026-11-25, pre-freeze. Three tables in `Data/` were not rebuilt by the
v0.4 physiology fixes (commit `4aae7cd`). **No registered analysis reads any
of them.** They are kept as they were committed so that the history stays
inspectable. Do not use them as inputs. The current table for each is named
below.

| Table | Why it was not rebuilt | Known stale rows | Use instead |
|---|---|---|---|
| `physio_coverage.csv` | Written by `Scripts/intero2023_physio_qc.R` from the raw Biopac `.acq` files, which are not on this machine. | sub-065 `task-localizer_run-1` has `n_volumes` 167; the re-converted run has 196. sub-076 has `n_volumes` NA, because the old block inventory missed it. sub-070 and sub-083 rows were computed before the cache-key fix and have not been re-verified. | `qc_ledger.csv`, `resp_runs.csv`, and the manifests' `has_rvt` / `has_petco2` / `rvt_source` |
| `resp_run_adherence.csv` | No generating script exists in the repo. It is an early per-run adherence summary. | No rows for sub-070 or sub-083. Values predate the v0.4 cohort cut-point change. | `resp_runs.csv` (`adhere_*`) and `resp_epochs.csv` (`adherence_verified`, which defines S4) |
| `pmu_wide_offsets.csv` | No generating script exists in the repo. Its Biopac-derived offsets are the reference for the PMU alignment check in `Scripts/intero2023_pmu_rvt.R`. | sub-070 `task-breathing_run-1`: the offset came from the wrong recording (cache-key bug). The alignment check already excludes it with its `sharp >= 0.05` filter, and `Analysis/pmu_rvt_validation.md` documents that. | Keep as the validation reference. Nothing downstream changes. |

Rebuilding `physio_coverage.csv` would need the raw `.acq` files. Rebuilding
the other two would need new code written to reproduce tables whose original
code is lost. Neither would change any registered N or regressor, so neither
is planned.
