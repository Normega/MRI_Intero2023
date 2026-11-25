# intero2023_covariates.R
# Builds the first-level covariate sheets specified in the preregistration §7:
# one events file and one confounds file per usable breathing run, aligned
# row-for-row with the fMRIPrep confounds.
#
# Produces, under Repo/FirstLevel/:
#   <sub>_<ent>_desc-model_events.tsv     onset/duration/condition + modulators
#   <sub>_<ent>_desc-model_confounds.tsv  one row per volume, GLM-ready
#   firstlevel_manifest.csv               one row per run: counts, gates, physio
#   README.md                             column conventions
#
# CONDITION MAPPING (prereg §5/§7). "No Change" responses on change trials are
# UNAWARE, not errors of commission; No Change *trials* are their own condition
# and are never treated as correct rejections (behavioural accuracy on them is
# at chance). Wrong-direction is kept separate so it contaminates neither.
#
#   paced_aware      change trial, response == instructed direction
#   paced_unaware    change trial, response == "No Change"
#   paced_wrongdir   change trial, directional response, wrong direction  [no interest]
#   paced_nochange   No Change trial, any response
#   free             free-breathing epoch
#   rating           post-epoch rating window (paced and free alike)
#   subrating        valence/sub-rating events                            [no interest]
#
# NA POLICY, deliberately explicit because it silently shapes the GLM:
#   * fMRIPrep motion derivatives and FD are NA on row 1 -> 0. Standard.
#   * Physio columns are demeaned within run and then NA -> 0, i.e. mean
#     imputation. The raw series with visible gaps stay in Repo/Data/regressors;
#     these sheets are the GLM-ready version. A physio column whose coverage is
#     below 50% of volumes is DROPPED for that run (mean-imputing most of a
#     regressor makes it a constant) and the manifest records it.
#   * Modulators (confidence, arousal) are written RAW; mean-centring is done
#     at design time per the prereg, so centring is over the trials actually
#     modelled, not over trials that may be excluded later.
#
# MOTION GATE (prereg §4): mean FD > 0.5 mm, or > 25% of volumes with
# FD > 0.5 mm. Flagged in the manifest and the QC ledger; files are still
# written so the gate remains reviewable rather than silently applied.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_ev   <- file.path(dir_out, "events")
dir_reg  <- file.path(dir_out, "regressors")
dir_fl   <- file.path(dir_root, "Repo", "FirstLevel")
dir_ready<- file.path(dir_root, "Data", "Ready")
dir_R    <- file.path(dir_root, "Repo", "R")
TR <- 2
FD_SPIKE <- 0.5           # mm, spike threshold AND gate threshold (prereg)
FD_GATE_FRAC <- 0.25
PHYSIO_MIN_COV <- 0.5

for (f in c("intero_ids.R", "intero_qc.R")) source(file.path(dir_R, f))
dir.create(dir_fl, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

# EXCLUDE of sub-082/sub-104 removed 2026-11-24: the label shift exists only in the
# old Data/Ready snapshot (prereg v0.2 §4). Kept empty so the filter below reads the same.
EXCLUDE <- character(0)
TASK_MAP <- c(breath1 = "task-breathing_run-1", breath2 = "task-breathing_run-2")
conf_idx <- tibble(path = list.files(dir_ready, pattern = "desc-confounds_timeseries[.]tsv$",
                                     recursive = TRUE, full.names = TRUE)) %>%
  mutate(bn = basename(path),
         subject = str_match(bn, "^(sub-[0-9]+)_")[, 2],
         rtask = str_match(bn, "_task-([a-z0-9]+)_")[, 2],
         bids_entity = unname(TASK_MAP[rtask])) %>%
  dplyr::filter(!is.na(bids_entity), !subject %in% EXCLUDE)

demean0 <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0, length(x)))
  x[ok] <- x[ok] - mean(x[ok])
  x[!ok] <- 0
  x
}
fit_len <- function(x, n) { length(x) <- n; x }   # pads with NA or trims

manifest <- list(); qc <- list()

for (i in seq_len(nrow(conf_idx))) {
  sb <- conf_idx$subject[i]; ent <- conf_idx$bids_entity[i]
  evf <- file.path(dir_ev, sprintf("%s_%s_events.tsv", sb, ent))
  if (!file.exists(evf)) next

  # ---- events ---------------------------------------------------------------
  ev <- suppressWarnings(read_tsv(evf, show_col_types = FALSE, progress = FALSE,
                                  col_types = cols(.default = col_character()))) %>%
    mutate(onset = as.numeric(onset), duration = as.numeric(duration),
           confidence = suppressWarnings(as.numeric(confidence)),
           arousal = suppressWarnings(as.numeric(arousal)))
  ev$condition <- with(ev, dplyr::case_when(
    trial_type == "paced_breathing" & change_direction == "No Change" ~ "paced_nochange",
    trial_type == "paced_breathing" & is.na(response)                 ~ "paced_noresponse",
    trial_type == "paced_breathing" & response == "No Change"         ~ "paced_unaware",
    trial_type == "paced_breathing" & response == change_direction    ~ "paced_aware",
    trial_type == "paced_breathing"                                   ~ "paced_wrongdir",
    trial_type == "free_breathing"                                    ~ "free",
    trial_type == "rating"                                            ~ "rating",
    trial_type == "valence_rating"                                    ~ "subrating",
    TRUE                                                              ~ NA_character_))
  ev_out <- ev %>%
    dplyr::filter(!is.na(condition), is.finite(onset), is.finite(duration)) %>%
    transmute(onset = round(onset, 3), duration = round(duration, 3), condition,
              epoch_index,
              confidence = ifelse(grepl("^paced|^free", condition), confidence, NA_real_),
              arousal    = ifelse(grepl("^paced|^free", condition), arousal, NA_real_),
              change_direction, salience, response)

  # ---- confounds ------------------------------------------------------------
  d <- suppressWarnings(read_tsv(conf_idx$path[i], show_col_types = FALSE,
                                 progress = FALSE,
                                 col_types = cols(.default = col_character())))
  n_vol <- nrow(d)
  motion <- c(paste0(rep(c("trans", "rot"), each = 3), "_", c("x", "y", "z")))
  keep <- c(motion, paste0(motion, "_derivative1"), "framewise_displacement",
            "csf", "white_matter", grep("^cosine", names(d), value = TRUE),
            grep("^non_steady_state_outlier", names(d), value = TRUE))
  keep <- intersect(keep, names(d))
  cf <- d[, keep] %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  # row-1 NAs on derivatives/FD -> 0
  cf <- cf %>% mutate(across(everything(), ~ ifelse(is.finite(.x), .x, 0)))

  fd <- cf$framewise_displacement
  spikes <- which(fd > FD_SPIKE)
  if (length(spikes)) for (k in seq_along(spikes)) {
    v <- rep(0L, n_vol); v[spikes[k]] <- 1L
    cf[[sprintf("spike_%02d", k)]] <- v
  }

  # ---- physio ---------------------------------------------------------------
  physio_note <- c()
  rf <- file.path(dir_reg, sprintf("%s_%s_resp.tsv", sb, ent))
  pf <- file.path(dir_reg, sprintf("%s_%s_petco2.tsv", sb, ent))
  has_rvt <- FALSE; has_co2 <- FALSE
  if (file.exists(rf)) {
    r <- suppressWarnings(read_tsv(rf, show_col_types = FALSE, progress = FALSE,
                                   na = "n/a"))
    if (abs(nrow(r) - n_vol) > 2) physio_note <- c(physio_note,
      sprintf("resp rows %d vs %d vols", nrow(r), n_vol))
    for (cc in c("rvt_rrf_z", "rvt_z", "resp_rate_brpm")) {
      v <- fit_len(r[[cc]], n_vol)
      if (mean(is.finite(v)) >= PHYSIO_MIN_COV) {
        cf[[sub("_brpm$", "_c", cc)]] <- demean0(v); has_rvt <- TRUE
      }
    }
  }
  if (file.exists(pf)) {
    p <- suppressWarnings(read_tsv(pf, show_col_types = FALSE, progress = FALSE,
                                   na = "n/a"))
    if (abs(nrow(p) - n_vol) > 2) physio_note <- c(physio_note,
      sprintf("petco2 rows %d vs %d vols", nrow(p), n_vol))
    co2_cols <- grep("^petco2_hrf_lag", names(p), value = TRUE)
    cov_ok <- vapply(co2_cols, function(cc)
      mean(is.finite(fit_len(p[[cc]], n_vol))) >= PHYSIO_MIN_COV, logical(1))
    if (all(cov_ok) && length(co2_cols)) {          # the basis enters whole or not at all
      for (cc in co2_cols) cf[[cc]] <- demean0(fit_len(p[[cc]], n_vol))
      has_co2 <- TRUE
    }
  }

  # ---- write ----------------------------------------------------------------
  write_tsv(ev_out, file.path(dir_fl, sprintf("%s_%s_desc-model_events.tsv", sb, ent)),
            na = "n/a")
  write_tsv(cf %>% mutate(across(everything(), ~ round(.x, 6))),
            file.path(dir_fl, sprintf("%s_%s_desc-model_confounds.tsv", sb, ent)),
            na = "n/a")

  fd_mean <- mean(fd); fd_frac <- mean(fd > FD_SPIKE)
  gate_fail <- fd_mean > FD_SPIKE | fd_frac > FD_GATE_FRAC
  cc_tab <- table(ev_out$condition)
  g <- function(x) if (x %in% names(cc_tab)) as.integer(cc_tab[[x]]) else 0L
  ev_end <- max(ev_out$onset + ev_out$duration)

  manifest[[length(manifest) + 1]] <- tibble(
    subject = sb, bids_entity = ent, n_vol = n_vol, scan_dur_s = n_vol * TR,
    events_end_s = round(ev_end, 1),
    n_aware = g("paced_aware"), n_unaware = g("paced_unaware"),
    n_wrongdir = g("paced_wrongdir"), n_nochange = g("paced_nochange"),
    n_noresponse = g("paced_noresponse"), n_free = g("free"),
    n_rating = g("rating") + g("subrating"),
    fd_mean = round(fd_mean, 3), fd_frac_gt_thresh = round(fd_frac, 3),
    n_spikes = length(spikes), motion_gate_fail = gate_fail,
    has_rvt = has_rvt, has_petco2 = has_co2,
    note = if (length(physio_note)) paste(physio_note, collapse = "; ") else NA_character_)

  qc[[length(qc) + 1]] <- qc_entry(sb, "imaging", "motion_gate",
    if (gate_fail) "fail" else if (fd_mean > 0.3) "degraded" else "good",
    entity = ent, value = fd_mean,
    note = sprintf("meanFD=%.3f fracFD>%.1f=%.2f spikes=%d", fd_mean, FD_SPIKE,
                   fd_frac, length(spikes)),
    source = "fMRIPrep framewise_displacement; prereg gate")
}

mf <- bind_rows(manifest)
write_csv(mf, file.path(dir_fl, "firstlevel_manifest.csv"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), bind_rows(qc)),
         file.path(dir_out, "qc_ledger.csv"))

writeLines(c(
  "# FirstLevel covariate sheets",
  "",
  "One `_desc-model_events.tsv` + `_desc-model_confounds.tsv` pair per usable",
  "breathing run. Confounds are row-aligned with the fMRIPrep confounds file",
  "(one row per volume, TR = 2 s).",
  "",
  "## events columns",
  "- condition: paced_aware / paced_unaware / paced_wrongdir / paced_nochange /",
  "  paced_noresponse / free / rating / subrating.",
  "  wrongdir, noresponse and subrating are regressors of NO interest (prereg §5/§7).",
  "- confidence, arousal: RAW 1-7, present on breathing epochs only.",
  "  Mean-centre at design time, over the trials actually modelled.",
  "",
  "## confounds columns",
  "- motion (6) + derivative1 (6), framewise_displacement, csf, white_matter,",
  "  cosine*, non_steady_state_outlier* : copied from fMRIPrep; row-1 NA -> 0.",
  "- spike_NN : one column per volume with FD > 0.5 mm.",
  "- rvt_rrf_z, rvt_z, resp_rate_c : belt physio, demeaned, NA -> 0 (mean",
  "  imputation). Raw series with visible gaps: Repo/Data/regressors/*_resp.tsv.",
  "- petco2_hrf_lag00/10/20/30 : the registered PetCO2 basis, demeaned, NA -> 0.",
  "  The basis is present whole or absent whole (never a subset of lags).",
  "- A physio column with <50% coverage in a run is omitted for that run;",
  "  see firstlevel_manifest.csv (has_rvt / has_petco2).",
  "",
  "## gates",
  "- motion_gate_fail in the manifest applies the preregistered rule",
  "  (mean FD > 0.5 mm or >25% of volumes above 0.5 mm). Files for gated runs",
  "  are still written; exclusion is applied at model time, reviewably.",
  ""), file.path(dir_fl, "README.md"))

cat(sprintf("\nruns written: %d  participants: %d\n", nrow(mf), n_distinct(mf$subject)))
cat("\n-- condition counts across all runs --\n")
print(as.data.frame(mf %>% summarise(across(starts_with("n_"), sum))))
cat("\n-- motion gate --\n")
print(as.data.frame(mf %>% count(motion_gate_fail)))
cat("\n-- physio availability --\n")
print(as.data.frame(mf %>% count(has_rvt, has_petco2)))
cat("\n-- runs with alignment notes --\n")
print(as.data.frame(mf %>% dplyr::filter(!is.na(note)) %>%
                      select(subject, bids_entity, note)))
cat(sprintf("\nRun ended: %s\n", Sys.time()))
