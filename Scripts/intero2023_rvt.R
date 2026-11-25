# intero2023_rvt.R
# Per-volume RVT and breathing rate from the Biopac belt, with the same
# global-signal sanity check applied to PetCO2.
#
# THE SANITY CHECK, AND WHY ITS SIGN IS THE INTERESTING PART
#
# RVT convolved with the Birn RRF is, by construction, a prediction of the BOLD
# signal: the RRF is the impulse response from RVT to BOLD. So the correlation
# with the global signal must be POSITIVE. PetCO2 through the same kernel came
# out negative, which is what exposed the kernel as wrong for CO2. If RVT comes
# out positive through that same kernel, the sign asymmetry is confirmed from
# both directions and neither result is a normalisation artefact.
#
# Belt runs graded "fail" are excluded outright -- RVT is an AMPLITUDE measure,
# so it is even less robust to a flat trace than rate is.
#
# Produces:
#   Repo/Data/rvt_runs.csv                        per run: coupling and quality
#   Repo/Data/regressors/<sub>_<entity>_resp.tsv  rvt_z, rvt_rrf_z, rate_brpm
#   Repo/Figures/rvt_check.png
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
  library(signal)
})

# Paths default to the Windows shared drive; INTERO_ROOT / RESPKIT_R / INTERO_FMRIPREP
# override them (2026-11-24: lets the same script run against the Linux tree).
dir_root  <- Sys.getenv("INTERO_ROOT", "I:/Shared drives/MRI_Intero2023")
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_fig   <- file.path(dir_root, "Repo", "Figures")
dir_reg   <- file.path(dir_out, "regressors")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")
dir_ready <- Sys.getenv("INTERO_FMRIPREP", file.path(dir_root, "Data", "Ready"))
TR <- 2

rk <- Sys.getenv("RESPKIT_R", "I:/Shared drives/TechnicalReference/respkit/R")
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R", "intero_resp.R", "intero_rrf.R",
            "intero_rvt.R")) source(file.path(dir_R, f))
dir.create(dir_reg, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

# EXCLUDE of sub-082/sub-104 removed 2026-11-24: their run labels are shifted only in
# the old Data/Ready snapshot; physio_blocks maps them by scanner slot and the BIDS
# tree of record agrees (prereg v0.2 §4). Kept as an empty vector so the filters
# below still read the same.
EXCLUDE <- character(0)
TASK_MAP <- c(loc1 = "task-localizer_run-1", loc2 = "task-localizer_run-2",
              breath1 = "task-breathing_run-1", breath2 = "task-breathing_run-2",
              auth1 = "task-authenticity")
conf <- tibble(path = list.files(dir_ready, pattern = "desc-confounds_timeseries[.]tsv$",
                                 recursive = TRUE, full.names = TRUE)) %>%
  mutate(bn = basename(path),
         subject = str_match(bn, "^(sub-[0-9]+)_")[, 2],
         bids_entity = confounds_entity(bn, TASK_MAP)) %>%   # old or BIDS naming
  dplyr::filter(!is.na(bids_entity), !subject %in% EXCLUDE)

blocks <- read_csv(file.path(dir_out, "physio_blocks.csv"), show_col_types = FALSE,
                   progress = FALSE) %>%
  dplyr::filter(kind == "bold", !is.na(bids_entity), is.na(dup_slot) | !dup_slot)

cfl <- list.files(dir_cache, pattern = "_physio[.]rds$", full.names = TRUE)
ckey <- physio_cache_index(cfl) %>% rename(cpath = path)   # see R/intero_ids.R

rows <- list(); qc <- list(); nreg <- 0L

for (j in seq_len(nrow(blocks))) {
  sb <- blocks$subject[j]; ent <- blocks$bids_entity[j]
  if (sb %in% EXCLUDE) next
  bp <- blocks$part[j]; if (is.na(bp)) bp <- "-"
  ck <- ckey %>% dplyr::filter(subject == sb, part == bp)
  # no exact part match -> skip; the old slice(1) fallback read a different recording
  if (!nrow(ck)) next
  rr <- try(readRDS(ck$cpath[1]), silent = TRUE)
  if (inherits(rr, "try-error") || is.null(rr$belt)) next

  fs <- rr$meta$fs; belt <- as.numeric(rr$belt$signal)
  s <- max(1L, round(blocks$onset_s[j] * fs)); e <- min(length(belt), round(blocks$end_s[j] * fs))
  if (e - s < fs * 60) next
  a <- belt_analyse(belt[s:e], fs, id = paste(sb, ent, sep = "_"))
  if (is.null(a)) next
  grade <- belt_grade(a$quality)
  if (grade == "fail") next            # amplitude measure; a flat trace is unusable

  cp <- conf %>% dplyr::filter(subject == sb, bids_entity == ent)
  n_vol <- if (nrow(cp)) {
    length(readLines(cp$path[1], warn = FALSE)) - 1L
  } else blocks$n_volumes[j]
  if (is.na(n_vol) || n_vol < 60) next

  rvt   <- rvt_series(a$bt, n_vol, TR)
  rate  <- rate_series(a$bt, n_vol, TR)
  rvtc  <- rvt_convolved(a$bt, n_vol, TR)

  g <- if (nrow(cp)) {
    d <- suppressWarnings(read_tsv(cp$path[1], show_col_types = FALSE, progress = FALSE,
                                   col_types = cols(.default = col_character())))
    suppressWarnings(as.numeric(d$global_signal))
  } else NULL

  r_raw <- if (!is.null(g)) cor_band(rvt,  g, 1 / TR) else NA_real_
  r_rrf <- if (!is.null(g)) cor_band(rvtc, g, 1 / TR) else NA_real_

  rows[[length(rows) + 1]] <- tibble(
    subject = sb, bids_entity = ent, task = str_match(ent, "task-([a-z]+)")[, 2],
    belt_grade = grade, n_vol = n_vol,
    n_breaths = sum(a$bt$keep %in% TRUE),
    cov_rvt = mean(is.finite(rvt)), cov_rate = mean(is.finite(rate)),
    rate_med = median(rate, na.rm = TRUE), rate_sd = sd(rate, na.rm = TRUE),
    r_rvt_raw = r_raw, r_rvt_rrf = r_rrf)

  out <- tibble(rvt_z = round(z_run(rvt), 6),
                rvt_rrf_z = round(z_run(rvtc), 6),
                resp_rate_brpm = round(rate, 4))
  write_tsv(out, file.path(dir_reg, sprintf("%s_%s_resp.tsv", sb, ent)), na = "n/a")
  nreg <- nreg + 1L

  qc[[length(qc) + 1]] <- qc_entry(sb, "physio", "rvt_coverage",
    if (mean(is.finite(rvt)) > 0.9) "good" else if (mean(is.finite(rvt)) > 0.6) "degraded" else "fail",
    entity = ent, value = mean(is.finite(rvt)),
    note = sprintf("%d kept breaths, belt %s", sum(a$bt$keep %in% TRUE), grade),
    source = "breath_table rvt interpolated to volumes")
}

rt <- bind_rows(rows)
write_csv(rt, file.path(dir_out, "rvt_runs.csv"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), bind_rows(qc)),
         file.path(dir_out, "qc_ledger.csv"))

cat(sprintf("\nruns with RVT: %d  participants: %d  regressor files: %d\n",
            nrow(rt), n_distinct(rt$subject), nreg))
cat("\n-- coverage --\n")
cat(sprintf("  RVT defined on median %.0f%% of volumes; rate on %.0f%%\n",
            100 * median(rt$cov_rvt), 100 * median(rt$cov_rate)))
cat("\n-- breathing rate per volume (brpm) --\n")
print(as.data.frame(rt %>% group_by(task) %>% summarise(n = n(),
  rate = round(median(rate_med, na.rm = TRUE), 1),
  within_run_sd = round(median(rate_sd, na.rm = TRUE), 2), .groups = "drop")))

cat("\n== SANITY CHECK: RVT through the RRF should predict global BOLD POSITIVELY ==\n")
cat(sprintf("  RVT raw       : median r = %+.3f   positive in %.0f%% of runs\n",
            median(rt$r_rvt_raw, na.rm = TRUE), 100 * mean(rt$r_rvt_raw > 0, na.rm = TRUE)))
cat(sprintf("  RVT (x) RRF   : median r = %+.3f   positive in %.0f%% of runs\n",
            median(rt$r_rvt_rrf, na.rm = TRUE), 100 * mean(rt$r_rvt_rrf > 0, na.rm = TRUE)))
cat("\n  compare PetCO2 through the SAME kernel, which must be and was negative.\n")
cat("\n-- by task --\n")
print(as.data.frame(rt %>% group_by(task) %>% summarise(n = n(),
  r_raw = round(median(r_rvt_raw, na.rm = TRUE), 3),
  r_rrf = round(median(r_rvt_rrf, na.rm = TRUE), 3), .groups = "drop")))

png(file.path(dir_fig, "rvt_check.png"), width = 1400, height = 460, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
hist(rt$r_rvt_rrf, breaks = 25, col = "grey80", border = "white",
     main = "RVT (x) RRF vs global signal", xlab = "r"); abline(v = 0, col = "red", lwd = 2)
hist(rt$rate_med, breaks = 25, col = "grey80", border = "white",
     main = "median breathing rate", xlab = "brpm")
boxplot(r_rvt_rrf ~ task, data = rt, col = "grey85", main = "by task", ylab = "r", xlab = "")
abline(h = 0, col = "red")
par(op); dev.off()
cat(sprintf("\nRun ended: %s\n", Sys.time()))
