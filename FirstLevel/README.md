# FirstLevel covariate sheets

One `_desc-model_events.tsv` + `_desc-model_confounds.tsv` pair per usable
breathing run. Confounds are row-aligned with the fMRIPrep confounds file
(one row per volume, TR = 2 s).

## events columns
- condition: paced_aware / paced_unaware / paced_wrongdir / paced_nochange /
  paced_noresponse / free / rating / subrating.
  wrongdir, noresponse and subrating are regressors of NO interest (prereg §5/§7).
- confidence, arousal: RAW 1-7, present on breathing epochs only.
  Mean-centre at design time, over the trials actually modelled.

## confounds columns
- motion (6) + derivative1 (6), framewise_displacement, csf, white_matter,
  cosine*, non_steady_state_outlier* : copied from fMRIPrep; row-1 NA -> 0.
- spike_NN : one column per volume with FD > 0.5 mm.
- rvt_rrf_z, rvt_z, resp_rate_c : belt physio, demeaned, NA -> 0 (mean
  imputation). Raw series with visible gaps: Repo/Data/regressors/*_resp.tsv.
- petco2_hrf_lag00/10/20/30 : the registered PetCO2 basis, demeaned, NA -> 0.
  The basis is present whole or absent whole (never a subset of lags).
- A physio column with <50% coverage in a run is omitted for that run;
  see firstlevel_manifest.csv (has_rvt / has_petco2).

## gates
- motion_gate_fail in the manifest applies the preregistered rule
  (mean FD > 0.5 mm or >25% of volumes above 0.5 mm). Files for gated runs
  are still written; exclusion is applied at model time, reviewably.

