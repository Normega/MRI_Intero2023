# De-identification review ledger

Files listed in `.gitignore` under the review rule, why each is held back, and
what would clear it. Nothing here is lost — everything remains on the shared
drive — these files are simply not entering **git history** until cleared,
because scrubbing a committed file later requires a history rewrite.

| File | Why held | What would clear it |
|---|---|---|
| `Data/screener_clean.csv` | Raw Qualtrics export: MRI safety screening (implants, medications, mental-health history), weight, free-text answers | Column-level review; keep only analysis variables; drop free text |
| `Data/postscan_clean.csv` | Raw post-scan responses, unreviewed free text | Same |
| `Data/run_notes.csv` | RA session notes; free text may contain identifying anecdotes | Read-through and redaction |
| `Data/Data.zip` | Contents not inventoried | Inventory, then decide per file |
| `Data/log_sessions.csv`, `log_runs.csv`, `log_run_attempts.csv`, `log_qc_flags.csv` | Wallclock date-times per subject ID (scan date is an indirect identifier) | Shift dates to relative days from each participant's session, or drop time columns |
| `Data/physio_blocks.csv` | `abs_time` column (wallclock) | Same |
| `Data/physio_session_status.csv` | Session dates | Same |

**Cleared as de-identified** (committed): all code (`R/`, `Scripts/`, `Tests/`),
`Preregistration/`, `Figures/`, `FirstLevel/` (relative time only),
`Data/events/` (onsets relative to scan start), and derived tables keyed by
`sub-XXX` with no dates: `breathing_trials.csv`, `behaviour_*.csv`,
`maia_scores.csv`, `petco2_runs.csv`, `petco2_rrf_fit.csv`,
`petco2_fixedlag_fit.csv`, `resp_runs.csv`, `resp_epochs.csv`,
`resp_calibration.csv`, `resp_run_adherence.csv`, `rvt_runs.csv`,
`physio_coverage.csv`, `pmu_validation.csv`, `pmu_wide_offsets.csv`,
`prereg_sample.csv`, `events_coverage.csv`, `qc_ledger.csv`,
`screener_codebook.csv` (labels only, no responses), `Data/regressors/`.

The capnograph-outage calendar dates quoted in prose (preregistration,
comments) are equipment dates, not participant identifiers, and stay.

## Public-release sweep (2026-08-20)

Systematic scan of every committed file before the repository went public:

| Class | Method | Result |
|---|---|---|
| PID + scan-date strings | regex `P###_20######` over tracked files | Two code files used `P052_20231025` as a parsing example — a synthetic pairing (not P052's real date, but the date itself was another participant's). Replaced with impossible date `20991231`; tests re-verified. |
| Email addresses | regex over tracked files | Only the PI and REB office addresses inside consent-form text in `screener_codebook.csv` — public information, retained. |
| IP addresses | regex | none |
| Date-times in data files | regex over `*.csv`/`*.tsv` | none (all dated tables are excluded by `.gitignore`) |
| Identifier-like columns | header scan for name/dob/sex/age/phone/address/email | none (matches were substrings: cover*age*, st*age*, sali*ence*) |
| Remaining calendar dates | regex over docs | Collection window (2023-09 – 2024-06) and capnograph-outage bounds (2023-10-12 – 2024-01-18) in the preregistration — equipment/methods reporting that identifies no individual; retained deliberately. |

Standing rule for the public repo: nothing joins a `sub-XXX` key to a
calendar date, and the `.gitignore` boundary may only move after the
column-level reviews above.
