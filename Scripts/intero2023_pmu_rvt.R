# intero2023_pmu_rvt.R
# Siemens PMU respiration as a FALLBACK source of the RVT regressor, for
# breathing runs whose Biopac belt yields none. Added 2026-11-24, before the
# registration freeze and before any first-level model was fitted.
#
# SOURCE PRIORITY (per run, analysable 90-run set only)
#   1. Biopac belt, as built by intero2023_rvt.R (belt not graded "fail") and
#      with rvt_rrf_z defined on >= 50% of volumes -> rvt_source = "biopac".
#   2. Otherwise the scanner's PMU respiratory cushion, if it passes every gate
#      below -> rvt_source = "pmu".
#   3. Otherwise no RVT for the run.
#   A run never mixes sources; a participant may (run 1 Biopac, run 2 PMU).
#   Each series is z-scored within run either way, so amplitude units do not
#   cross runs.
#
# ALIGNMENT
#   scan start (s into the PMU log) = AcquisitionTime - LogStartMDHTime + 0.62
#   AcquisitionTime is the BOLD sidecar's (BIDS tree of record); both are on the
#   scanner's MDH clock. The +0.62 s (PMU_CAL_S, R/intero_pmu.R) is the negated
#   median residual against the Biopac-trigger-locked offsets of
#   Data/pmu_wide_offsets.csv over the 26 dual-belt breathing runs examined on
#   2026-11-24 (prediction r = 0.997; residual median -0.62 s, IQR -0.77 to
#   -0.27, range -1.03 to +1.83 s). The constant is FIXED; this script
#   re-estimates it over every available run only as a check and reports both.
#
# GATES (all must hold)
#   PMU log covers >= 95% of the scan window      PMU_MIN_COVERED
#   band_ratio >= 0.90                            PMU_MIN_BAND  (= registered
#                                                 "good" Biopac criterion)
#   < 10% of raw 12-bit samples in the scan       PMU_MAX_SAT
#     window on a rail (0 or 4095)
#   >= 5 kept breaths                             RVT_MIN_BREATHS (intero_rvt.R)
#   rvt defined on >= 50% of volumes              PMU_MIN_RVT_COV (= sheets' rule)
#   Polarity: inverted (PMU_SIGN = -1), chosen on agreement with Biopac RVT.
#   NO gate is placed on the global-signal correlation; it is reported only, as
#   for the Biopac route.
#
# PROVENANCE OF THE SATURATION GATE -- stated because it was not pre-set.
#   A gate of < 5% saturation over the WHOLE PMU file was set before the
#   fallback runs were scored. Scoring printed, for every fallback run, its
#   metrics including the RVT(x)RRF-vs-global-signal r, and the 5% gate failed
#   17 of 28 runs on saturation alone. Because the dual-belt validation runs
#   themselves reached 27% saturation, the gate was then re-derived from the
#   VALIDATION set only (PMU-Biopac RVT agreement by saturation bin: median r
#   0.63 / 0.48 / 0.22 / 0.34 for <5 / 5-10 / 10-15 / >=15 %), and re-applied
#   to saturation measured in the SCAN WINDOW rather than the whole file. The
#   cut (< 10%) was not chosen on the fallback runs' metrics, but those metrics
#   had been seen when it was chosen. Runs at 10-15% (7) are NOT used.
#   LATER CHECK (same day): extending the validation to all 51 dual-belt
#   breathing runs does not reproduce the decision set's RVT agreement (median
#   r 0.17 overall, 0.05 on the 25 runs added; 0.15 in the 5-10% bin vs 0.48).
#   The decision set was selected for sharp Biopac-PMU cross-correlation and
#   over-states agreement. Global-signal coupling of PMU RVT holds in every
#   subset (median +0.34 vs +0.22 for the Biopac RVT on the same runs). See
#   Analysis/pmu_rvt_validation.md.
#
# INPUTS   PMU .resp files (INTERO_PMU, searched recursively); BIDS sidecars
#          (INTERO_BIDS); fMRIPrep derivatives of record (INTERO_FMRIPREP:
#          volume count, and global_signal for the reported check only);
#          Data/regressors/*_resp.tsv (Biopac); FirstLevel manifests (run set);
#          Data/pmu_wide_offsets.csv (Biopac-derived offsets, validation).
# PRODUCES Data/regressors/<sub>_<entity>_desc-pmu_resp.tsv  rvt_z, rvt_rrf_z,
#            resp_rate_brpm -- same columns as the Biopac file
#          Data/rvt_pmu_runs.csv           every fallback candidate, metrics + use
#          Analysis/pmu_rvt_alignment.csv  offset prediction vs Biopac offsets
#          Analysis/pmu_rvt_validation.csv dual-belt agreement, per run
#          No time of day is written anywhere: offsets and residuals only.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(signal); library(readr); library(dplyr); library(tibble); library(purrr)
  library(stringr)
})

dir_root  <- Sys.getenv("INTERO_ROOT", "I:/Shared drives/MRI_Intero2023")
dir_repo  <- file.path(dir_root, "Repo")
dir_out   <- file.path(dir_repo, "Data")
dir_reg   <- file.path(dir_out, "regressors")
dir_an    <- file.path(dir_repo, "Analysis")
dir_R     <- file.path(dir_repo, "R")
dir_pmu   <- Sys.getenv("INTERO_PMU", file.path(dir_root, "Data", "MRI"))
dir_bids  <- Sys.getenv("INTERO_BIDS", "/media/normega/FB_Cluster_Data/intero2024/bids")
dir_ready <- Sys.getenv("INTERO_FMRIPREP", "/media/normega/FB_Cluster_Data/intero2024/derivatives/fmriprep")
TR <- 2
rk <- Sys.getenv("RESPKIT_R", "I:/Shared drives/TechnicalReference/respkit/R")
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R", "intero_resp.R", "intero_rrf.R",
            "intero_rvt.R", "intero_pmu.R")) source(file.path(dir_R, f))
dir.create(dir_an, showWarnings = FALSE)

SLOT <- c("task-localizer_run-1", "task-localizer_run-2", "task-breathing_run-1",
          "task-breathing_run-2", "task-authenticity")
pmu_idx <- pmu_file_index(dir_pmu) %>%
  mutate(subject = pid_bids(pid), bids_entity = SLOT[protocol_slot])

gs_of <- function(sb, ent) {   # global_signal ONLY: the disclosed exception
  f <- file.path(dir_ready, sb, "func", sprintf("%s_%s_desc-confounds_timeseries.tsv", sb, ent))
  if (!file.exists(f)) return(NULL)
  as.numeric(read_tsv(f, col_select = "global_signal", show_col_types = FALSE,
                      progress = FALSE)$global_signal)
}
acq_of <- function(sb, ent) bids_acq_seconds(file.path(dir_bids, sb, "func",
                                             sprintf("%s_%s_bold.json", sb, ent)))
read_pmu <- function(sb, ent) {
  pth <- pmu_idx$path[pmu_idx$subject == sb & pmu_idx$bids_entity == ent]
  if (!length(pth)) NULL else pmu_read(pth[1])
}
cov50 <- function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- read_tsv(f, na = "n/a", show_col_types = FALSE, progress = FALSE)
  mean(is.finite(r$rvt_rrf_z)) >= PHYSIO_MIN_COV
}
PHYSIO_MIN_COV <- 0.5

cat(sprintf("Run started: %s\nPMU files indexed: %d\n", Sys.time(), nrow(pmu_idx)))

# ---- 1. alignment check against the Biopac-derived offsets ---------------
wide <- read_csv(file.path(dir_out, "pmu_wide_offsets.csv"), show_col_types = FALSE)
al <- pmap_dfr(list(wide$subject, wide$bids_entity, wide$off_best_s), function(sb, ent, ob) {
  p <- read_pmu(sb, ent); acq <- acq_of(sb, ent)
  if (is.null(p) || !is.finite(acq)) return(NULL)
  tibble(subject = sb, bids_entity = ent, offset_raw_s = round(pmu_scan_offset_raw(p, acq), 3),
         offset_biopac_s = ob)
}) %>% mutate(residual_s = round(offset_raw_s - offset_biopac_s, 3))
# A Biopac-derived offset is a measurement only when the Biopac-PMU
# cross-correlation has a distinct peak. REFERENCE_OK: sharp >= 0.05 (best peak
# beats the next one > 5 s away by >= 0.05 in r). This filter was set AFTER the
# all-run residuals were seen (2026-11-24): the only large residuals (-260 s,
# +7 s, +19 s) were sub-070 run-1, whose Biopac offset was computed from the
# wrong recording (the cache-key bug), and sub-089 runs 1-2 (sharp 0.04 /
# 0.003, |r| ~0.2: no peak). It changes neither PMU_CAL_S nor any gate; both
# summaries are reported.
al <- al %>% left_join(wide %>% select(subject, bids_entity, r_best, sharp), by = c("subject", "bids_entity")) %>%
  mutate(reference_ok = sharp >= 0.05)
write_csv(al, file.path(dir_an, "pmu_rvt_alignment.csv"))
for (lab in c("all runs", "reference_ok")) {
  a <- if (lab == "all runs") al else al %>% dplyr::filter(reference_ok)
  cat(sprintf("\n== alignment, %s: %d runs  r = %.4f  residual median %.2f s, IQR [%.2f, %.2f], range [%.2f, %.2f]\n",
              lab, nrow(a), cor(a$offset_raw_s, a$offset_biopac_s), median(a$residual_s),
              quantile(a$residual_s, .25), quantile(a$residual_s, .75), min(a$residual_s), max(a$residual_s)))
  cat(sprintf("   calibration implied: %+.2f s (fixed PMU_CAL_S = %+.2f s)\n", -median(a$residual_s), PMU_CAL_S))
}

# ---- 2. dual-belt validation ---------------------------------------------
# DECISION_SET: the 27 runs the method was judged on (2026-11-24). FULL: every
# breathing run with both a Biopac RVT file and a PMU file, reported as a check.
DECISION_SET <- c(
  "sub-059 task-breathing_run-1", "sub-059 task-breathing_run-2", "sub-061 task-breathing_run-1",
  "sub-062 task-breathing_run-2", "sub-064 task-breathing_run-2", "sub-065 task-localizer_run-1",
  "sub-065 task-breathing_run-1", "sub-065 task-breathing_run-2", "sub-066 task-breathing_run-1",
  "sub-066 task-breathing_run-2", "sub-072 task-breathing_run-1", "sub-076 task-breathing_run-1",
  "sub-077 task-breathing_run-1", "sub-077 task-breathing_run-2", "sub-080 task-breathing_run-1",
  "sub-084 task-breathing_run-1", "sub-084 task-breathing_run-2", "sub-092 task-breathing_run-1",
  "sub-092 task-breathing_run-2", "sub-094 task-breathing_run-2", "sub-098 task-breathing_run-1",
  "sub-098 task-breathing_run-2", "sub-100 task-breathing_run-1", "sub-104 task-breathing_run-1",
  "sub-106 task-breathing_run-1", "sub-108 task-breathing_run-2", "sub-116 task-breathing_run-2")
bio_files <- tibble(path = list.files(dir_reg, pattern = "^sub-[0-9]+_task-[a-z_0-9-]+_resp[.]tsv$",
                                      full.names = TRUE)) %>%
  dplyr::filter(!grepl("_desc-pmu_", basename(path))) %>%   # Biopac files only
  mutate(subject = str_match(basename(path), "^(sub-[0-9]+)_")[, 2],
         bids_entity = str_match(basename(path), "^sub-[0-9]+_(.+)_resp[.]tsv$")[, 2]) %>%
  semi_join(pmu_idx, by = c("subject", "bids_entity")) %>%
  mutate(key = paste(subject, bids_entity),
         set = ifelse(key %in% DECISION_SET, "decision", "extra")) %>%
  dplyr::filter(set == "decision" | grepl("breathing", bids_entity))
val <- pmap_dfr(list(bio_files$path, bio_files$subject, bio_files$bids_entity, bio_files$set),
                function(bf, sb, ent, set) {
  bio <- read_tsv(bf, na = "n/a", show_col_types = FALSE, progress = FALSE)
  p <- read_pmu(sb, ent); acq <- acq_of(sb, ent)
  if (is.null(p) || !is.finite(acq)) return(NULL)
  r <- pmu_rvt_run(p, acq, nrow(bio), TR); if (!r$ok) return(NULL)
  g <- gs_of(sb, ent)
  tibble(subject = sb, bids_entity = ent, set = set, sat_win = round(r$sat_win, 4),
         band_ratio = round(r$a$quality$band_ratio, 4),
         r_rvt_rrf_pmu_vs_biopac = round(suppressWarnings(cor(r$rvtc, bio$rvt_rrf_z, use = "pairwise")), 4),
         r_rate_pmu_vs_biopac    = round(suppressWarnings(cor(r$rate, bio$resp_rate_brpm, use = "pairwise")), 4),
         r_gs_pmu    = if (is.null(g)) NA_real_ else round(cor_band(r$rvtc, g, 1 / TR), 4),
         r_gs_biopac = if (is.null(g)) NA_real_ else round(cor_band(bio$rvt_rrf_z, g, 1 / TR), 4))
})
write_csv(val, file.path(dir_an, "pmu_rvt_validation.csv"))
summ <- function(v) v %>% summarise(n = n(), agree_rvt_med = median(r_rvt_rrf_pmu_vs_biopac, na.rm = TRUE),
  agree_rate_med = median(r_rate_pmu_vs_biopac, na.rm = TRUE),
  gs_pmu_med = median(r_gs_pmu, na.rm = TRUE), gs_pmu_pos = mean(r_gs_pmu > 0, na.rm = TRUE),
  gs_bio_med = median(r_gs_biopac, na.rm = TRUE), gs_bio_pos = mean(r_gs_biopac > 0, na.rm = TRUE))
cat("\n== dual-belt validation ==\n-- decision set --\n"); print(as.data.frame(summ(val %>% dplyr::filter(set == "decision"))))
cat("-- all dual-belt breathing runs --\n"); print(as.data.frame(summ(val %>% dplyr::filter(grepl("breathing", bids_entity)))))
cat("-- by window saturation (all dual-belt breathing runs) --\n")
print(as.data.frame(val %>% dplyr::filter(grepl("breathing", bids_entity)) %>%
  mutate(bin = cut(sat_win, c(-Inf, .05, .10, .15, Inf), right = FALSE)) %>% group_by(bin) %>% summ()))

# ---- 3. fallback on the analysable set -----------------------------------
man <- bind_rows(read_csv(file.path(dir_repo, "FirstLevel", "firstlevel_manifest.csv"), show_col_types = FALSE),
                 read_csv(file.path(dir_repo, "FirstLevel", "firstlevel_manifest_candidates.csv"), show_col_types = FALSE)) %>%
  dplyr::filter(!motion_gate_fail, subject != "sub-056") %>% select(subject, bids_entity, n_vol)
cand <- man %>% dplyr::filter(!map2_lgl(subject, bids_entity, ~ cov50(file.path(dir_reg, sprintf("%s_%s_resp.tsv", .x, .y)))))
cat(sprintf("\n== fallback: %d of %d analysable runs lack a usable Biopac RVT ==\n", nrow(cand), nrow(man)))
fb <- pmap_dfr(list(cand$subject, cand$bids_entity, cand$n_vol), function(sb, ent, n_vol) {
  out <- tibble(subject = sb, bids_entity = ent, n_vol = n_vol)
  p <- read_pmu(sb, ent); acq <- acq_of(sb, ent)
  if (is.null(p) || !is.finite(acq)) return(out %>% mutate(used = FALSE, why = "no PMU file"))
  r <- pmu_rvt_run(p, acq, n_vol, TR)
  if (!r$ok) return(out %>% mutate(offset_s = r$offset_s, used = FALSE, why = "belt_analyse failed"))
  q <- r$a$quality; nb <- sum(r$a$bt$keep %in% TRUE); cv <- mean(is.finite(r$rvt)); g <- gs_of(sb, ent)
  why <- c(if (r$covered < PMU_MIN_COVERED) "log does not cover scan",
           if (!is.finite(q$band_ratio) || q$band_ratio < PMU_MIN_BAND) "band_ratio < 0.90",
           if (r$sat_win >= PMU_MAX_SAT) "saturation >= 10%",
           if (nb < RVT_MIN_BREATHS) "too few breaths", if (cv < PMU_MIN_RVT_COV) "rvt coverage < 50%")
  used <- !length(why)
  if (used) write_tsv(tibble(rvt_z = round(z_run(r$rvt), 6), rvt_rrf_z = round(z_run(r$rvtc), 6),
                             resp_rate_brpm = round(r$rate, 4)),
                      file.path(dir_reg, sprintf("%s_%s_desc-pmu_resp.tsv", sb, ent)), na = "n/a")
  out %>% mutate(offset_s = round(r$offset_s, 3), covered = round(r$covered, 4), sat_win = round(r$sat_win, 4),
                 band_ratio = round(q$band_ratio, 4), n_breaths = nb, rvt_cov = round(cv, 4),
                 rate_med_brpm = round(median(r$rate, na.rm = TRUE), 2),
                 r_gs = if (is.null(g)) NA_real_ else round(cor_band(r$rvtc, g, 1 / TR), 4),
                 used = used, why = if (used) NA_character_ else paste(why, collapse = "; "))
})
write_csv(fb, file.path(dir_out, "rvt_pmu_runs.csv"))
print(as.data.frame(fb %>% select(subject, bids_entity, sat_win, band_ratio, n_breaths, r_gs, used, why)))
cat(sprintf("\nPMU RVT written for %d runs (%d participants)\nRun ended: %s\n",
            sum(fb$used), n_distinct(fb$subject[fb$used]), Sys.time()))
