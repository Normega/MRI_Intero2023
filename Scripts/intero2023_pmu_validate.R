# intero2023_pmu_validate.R
# Validates the Siemens PMU respiration stream against the Biopac belt before
# any analysis is allowed to depend on it.
#
# THREE THINGS ARE BEING TESTED AT ONCE, on runs where the Biopac belt is known
# good (belt_grade == "good" in resp_runs.csv):
#
#   1. THE PARSER. Two independent transducers on the same chest during the
#      same minutes must produce correlated traces. A parsing error -- a
#      desynchronised flag filter, a missed comment block -- destroys that
#      correlation, so a high r is strong evidence the decode is right.
#   2. THE ALIGNMENT ASSUMPTION. pmu_window() takes the first n_vol * TR
#      seconds on the assumption that logging starts with the scan. If that
#      holds the cross-correlation peaks at lag ~0. If the peak is displaced,
#      the assumption is wrong and the offset is measurable.
#   3. PMU QUALITY, on its own terms: saturation against the 12-bit rails and
#      detected breath rate.
#
# Only if this passes does intero2023_pmu.R ingest the stream for real.
#
# Produces: Repo/Data/pmu_validation.csv, Repo/Figures/pmu_validate.png
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root  <- "I:/Shared drives/MRI_Intero2023"
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_fig   <- file.path(dir_root, "Repo", "Figures")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")
dir_mri   <- file.path(dir_root, "Data", "MRI")
TR <- 2

rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R", "intero_resp.R", "intero_pmu.R"))
  source(file.path(dir_R, f))
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

idx <- pmu_file_index(dir_mri) %>%
  filter(!is_pilot(pid)) %>%
  left_join(INTERO_RUNS %>% distinct(protocol_slot, bids_task, bids_run),
            by = "protocol_slot") %>%
  mutate(subject = pid_bids(pid),
         bids_entity = ifelse(is.na(bids_run), paste0("task-", bids_task),
                              paste0("task-", bids_task, "_run-", bids_run)))
cat(sprintf("PMU files indexed: %d across %d participants\n",
            nrow(idx), n_distinct(idx$pid)))

rn <- read_csv(file.path(dir_out, "resp_runs.csv"), show_col_types = FALSE,
               progress = FALSE)
blocks <- read_csv(file.path(dir_out, "physio_blocks.csv"), show_col_types = FALSE,
                   progress = FALSE) %>%
  filter(kind == "bold", !is.na(bids_entity), is.na(dup_slot) | !dup_slot)

# validate on GOOD-belt runs only: a dead Biopac trace cannot referee anything
tgt <- rn %>% filter(belt_grade == "good") %>%
  inner_join(idx, by = c("subject", "bids_entity"))
cat(sprintf("runs with PMU and a good Biopac belt: %d (%d participants)\n",
            nrow(tgt), n_distinct(tgt$subject)))

cache_files <- list.files(dir_cache, pattern = "_physio[.]rds$", full.names = TRUE)
ckey <- tibble(cpath = cache_files,
               subject = sub("(_part-[^_]*)?_physio[.]rds$", "", basename(cache_files)),
               part = ifelse(grepl("_part-", basename(cache_files)),
                             sub("^.*_part-([^_]*)_physio[.]rds$", "\\1",
                                 basename(cache_files)), "-"))

bandpass_resp <- function(x, fs, lo = 0.1, hi = 0.5) {
  ok <- is.finite(x); if (sum(ok) < 100) return(rep(NA_real_, length(x)))
  xi <- stats::approx(which(ok), x[ok], xout = seq_along(x), rule = 2)$y
  b <- signal::butter(2, c(lo, min(hi, 0.45 * fs)) / (fs / 2), type = "pass")
  as.numeric(signal::filtfilt(b, xi - mean(xi)))
}

LAGS <- seq(-20, 20, by = 0.5)   # seconds
out <- list()

for (i in seq_len(nrow(tgt))) {
  sb <- tgt$subject[i]; ent <- tgt$bids_entity[i]
  p <- pmu_read(tgt$path[i])
  if (is.null(p)) next

  bi <- blocks %>% filter(subject == sb, bids_entity == ent) %>% slice(1)
  if (!nrow(bi)) next
  bp <- bi$part[1]; if (is.na(bp)) bp <- "-"
  cf <- ckey %>% filter(subject == sb, part == bp)
  if (!nrow(cf)) cf <- ckey %>% filter(subject == sb) %>% slice(1)
  if (!nrow(cf)) next
  rr <- try(readRDS(cf$cpath[1]), silent = TRUE)
  if (inherits(rr, "try-error") || is.null(rr$belt)) next

  fsb <- rr$meta$fs
  belt <- as.numeric(rr$belt$signal)
  s <- max(1L, round(bi$onset_s[1] * fsb)); e <- min(length(belt), round(bi$end_s[1] * fsb))
  if (e - s < fsb * 60) next
  bseg <- belt[s:e]
  n_vol <- bi$n_volumes[1]; if (is.na(n_vol)) n_vol <- floor((e - s) / fsb / TR)

  w <- pmu_window(p$signal, p$fs, n_vol, TR)
  if (w$covered < 0.5) next

  # common 20 Hz grid: well above the respiratory band, cheap to shift
  FSC <- 20
  tb <- seq(0, (length(bseg) - 1) / fsb, by = 1 / FSC)
  tp <- seq(0, (length(w$signal) - 1) / p$fs, by = 1 / FSC)
  xb <- stats::approx(seq(0, (length(bseg) - 1) / fsb, length.out = length(bseg)),
                      bseg, xout = tb, rule = 2)$y
  xp <- stats::approx(seq(0, (length(w$signal) - 1) / p$fs, length.out = length(w$signal)),
                      w$signal, xout = tp, rule = 2)$y
  xb <- bandpass_resp(xb, FSC); xp <- bandpass_resp(xp, FSC)

  rl <- vapply(LAGS, function(L) {
    sh <- round(L * FSC)
    n <- min(length(xb), length(xp))
    a <- xb[seq_len(n)]; b <- xp[seq_len(n)]
    if (sh > 0)      { a <- a[(sh + 1):n]; b <- b[1:(n - sh)] }
    else if (sh < 0) { a <- a[1:(n + sh)]; b <- b[(-sh + 1):n] }
    if (length(a) < FSC * 60) return(NA_real_)
    suppressWarnings(stats::cor(a, b, use = "complete.obs"))
  }, numeric(1))

  bidx <- which.max(abs(rl))
  out[[length(out) + 1]] <- tibble(
    subject = sb, bids_entity = ent, slot = tgt$protocol_slot[i],
    n_pmu = p$n, pmu_dur_s = p$n / p$fs, covered = w$covered,
    sat_lo = p$sat_lo, sat_hi = p$sat_hi,
    best_lag_s = LAGS[bidx], r_best = rl[bidx], r_lag0 = rl[LAGS == 0],
    mdh_window_s = (p$mdh_stop_ms - p$mdh_start_ms) / 1000,
    belt_iqr = tgt$filtered_iqr[i])
  if (length(out) %% 10 == 0) cat(sprintf("  ...%d runs\n", length(out)))
}

v <- bind_rows(out)
write_csv(v, file.path(dir_out, "pmu_validation.csv"))
cat(sprintf("\nvalidated runs: %d (%d participants)\n", nrow(v), n_distinct(v$subject)))

cat("\n== 1. PARSER: does PMU agree with the Biopac belt? ==\n")
cat(sprintf("  |r| at best lag: median %.3f  IQR [%.3f, %.3f]\n",
            median(abs(v$r_best), na.rm = TRUE),
            quantile(abs(v$r_best), .25, na.rm = TRUE),
            quantile(abs(v$r_best), .75, na.rm = TRUE)))
cat(sprintf("  runs with |r| > 0.5: %.0f%%   > 0.7: %.0f%%\n",
            100 * mean(abs(v$r_best) > 0.5, na.rm = TRUE),
            100 * mean(abs(v$r_best) > 0.7, na.rm = TRUE)))
cat(sprintf("  sign of r at best lag: %+.0f%% positive\n",
            100 * mean(v$r_best > 0, na.rm = TRUE)))

cat("\n== 2. ALIGNMENT: is the peak at lag 0? ==\n")
print(round(quantile(v$best_lag_s, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE), 2))
cat(sprintf("  within +/-1 s of zero: %.0f%%   within +/-2 s: %.0f%%\n",
            100 * mean(abs(v$best_lag_s) <= 1, na.rm = TRUE),
            100 * mean(abs(v$best_lag_s) <= 2, na.rm = TRUE)))
cat(sprintf("  median |r| lost by forcing lag 0: %.3f\n",
            median(abs(v$r_best) - abs(v$r_lag0), na.rm = TRUE)))

cat("\n== 3. PMU QUALITY ==\n")
cat(sprintf("  saturation at 0: median %.2f%%   at 4095: median %.2f%%\n",
            100 * median(v$sat_lo), 100 * median(v$sat_hi)))
cat(sprintf("  runs with >5%% saturation either rail: %.0f%%\n",
            100 * mean(v$sat_lo > 0.05 | v$sat_hi > 0.05)))
cat(sprintf("  log outlasts acquisition (MDH window > 2x run): %.0f%% of runs\n",
            100 * mean(v$mdh_window_s > 2 * v$pmu_dur_s, na.rm = TRUE)))

png(file.path(dir_fig, "pmu_validate.png"), width = 1400, height = 460, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
hist(v$r_best, breaks = 30, col = "grey80", border = "white",
     main = "PMU vs Biopac, r at best lag", xlab = "r"); abline(v = 0, col = "red")
hist(v$best_lag_s, breaks = 40, col = "grey80", border = "white",
     main = "best lag", xlab = "s"); abline(v = 0, col = "red", lwd = 2)
plot(v$belt_iqr, abs(v$r_best), log = "x", pch = 16, col = "#00000070",
     xlab = "Biopac filtered_iqr", ylab = "|r|", main = "agreement vs Biopac quality")
par(op); dev.off()

cat(sprintf("\nRun ended: %s\n", Sys.time()))
