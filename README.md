# MRI_Intero2023 — reproducible analysis repository

Awareness of respiratory change during paced and free breathing, with
concurrent fMRI, respiration-belt, and expired-CO2 recording. This repository
holds the de-identified, shareable layer of the project: analysis code, derived
data products, the preregistration, and first-level modelling inputs. Raw data
(DICOM, .acq physiology, questionnaire exports) live outside the repository;
see `DEIDENT_REVIEW.md` for the boundary and its rationale.

## Layout

| Path | Contents |
|---|---|
| `R/` | Analysis modules: ID normalisation, log parsing, event extraction, physio (.acq, PMU), PetCO2, RVT, response-function convolution, QC ledger |
| `Scripts/` | Pipeline drivers, numbered by stage; each writes its outputs under `Data/` and appends to the QC ledger |
| `Tests/` | testthat suites (ground truth: hand-checked participant P130) |
| `Data/` | Derived, de-identified tables; `Data/events/` = BIDS-style event files; `Data/regressors/` = per-volume physiological regressors |
| `FirstLevel/` | GLM-ready covariate sheets, one events + one confounds file per run, row-aligned with fMRIPrep confounds; `firstlevel_manifest.csv` carries the registered gates |
| `Figures/` | QC and validation figures |
| `Preregistration/` | The imaging analysis plan, the OSF form text (`osf_template_answers.md`, filed as the update to https://osf.io/swc7p), and `deviations.md`, the deviations from that registration |

## Reproducing

Scripts assume R ≥ 4.6 with tidyverse, lme4, ordinal, signal, testthat, and the
`respkit` package (separate repository), plus Python `bioread` via reticulate
for .acq ingestion. Stage order and provenance are documented in each script
header. MRI preprocessing is fMRIPrep 24.1.1; its derivatives are consumed, not
regenerated, here.

## Preregistration

`Preregistration/intero2023_imaging_prereg.md` is the living document; the
registration of record is https://osf.io/swc7p (registered 2026-08-20 from
plan v0.1). Every later change is logged against the registered text in
`Preregistration/deviations.md`, with date, rationale, and whether any affected
result had been seen.
