# intero2023_resp_calibrate.R
# CALIBRATION PASS for the in-scanner respiration belt (Biopac DA100C).
#
# This script measures; it does not decide. It exists because respkit's
# BeforeYouRun.txt is explicit that the amplitude thresholds in resp_quality()
# are in the units of whatever transducer produced the signal, and are
# meaningless on a new one:
#
#   "Run resp_quality() across the whole sample with iqr_unusable/iqr_degraded
#    set to NA, plot filtered_iqr, pick a cut."
#
# Our belt is an uncalibrated strain gauge in Volts whose gain was set by hand
# each session, so its filtered_iqr has no reason to resemble the sample the
# defaults (0.05 / 0.20) were tuned on. Applying them blind would either pass
# every flatline or reject the whole cohort, and we would not be able to tell
# which.
#
# POLARITY IS MEASURED SEPARATELY BY RUN TYPE. respkit diagnoses polarity from
# the median duty cycle, on the assumption that inspiration is the shorter
# phase in quiet spontaneous breathing -- and warns that the test "is
# uninformative under a symmetric pacer". Our breathing runs are largely
# paced, so their duty cycle cannot settle the question. The localizer and
# authenticity runs are spontaneous, and those are what the polarity verdict
# is taken from; the breathing runs are reported alongside only to show the
# pacer's effect.
#
# Produces: Repo/Data/resp_calibration.csv   (one row per BOLD run)
#           Repo/Figures/resp_calibration.png
#           Repo/Figures/resp_traces_<sub>_<entity>.png  (visual check, step 4)
#
# Nothing here is written to the QC ledger: a calibration pass has no verdicts
# to record yet. The production pass (intero2023_resp.R) does that.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root  <- "I:/Shared drives/MRI_Intero2023"
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_fig   <- file.path(dir_root, "Repo", "Figures")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")

rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R")) source(file.path(dir_R, f))
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

# Trim is 15 s per end, not resp_quality's default 120. The default assumes a
# free-running recording whose ends may contain setup movement; our blocks are
# already bounded by the scanner trigger, and 120 s off each end of a 380 s
# localizer would leave 140 s to judge it on.
TRIM_S <- 15

blocks <- read_csv(file.path(dir_out, "physio_blocks.csv"),
                   show_col_types = FALSE, progress = FALSE) %>%
  filter(kind == "bold", !is.na(bids_entity), is.na(dup_slot) | !dup_slot)

cache_files <- list.files(dir_cache, pattern = "_physio[.]rds$", full.names = TRUE)
cache_key <- tibble(path = cache_files,
                    subject = sub("(_part-[^_]*)?_physio[.]rds$", "", basename(cache_files)),
                    part = ifelse(grepl("_part-", basename(cache_files)),
                                  sub("^.*_part-([^_]*)_physio[.]rds$", "\\1",
                                      basename(cache_files)), "-"))

rows <- list()

for (sb in sort(unique(blocks$subject))) {
  b_i <- blocks %>% filter(subject == sb)
  cf  <- cache_key %>% filter(subject == sb)
  if (!nrow(cf)) next

  for (j in seq_len(nrow(b_i))) {
    bp <- b_i$part[j]; if (is.na(bp)) bp <- "-"
    cfp <- cf %>% filter(part == bp); if (!nrow(cfp)) cfp <- cf[1, ]
    r <- try(readRDS(cfp$path[1]), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r$belt)) next

    fs   <- r$meta$fs
    belt <- as.numeric(r$belt$signal)
    s <- max(1L, round(b_i$onset_s[j] * fs))
    e <- min(length(belt), round(b_i$end_s[j] * fs))
    if (e - s < fs * 60) next
    seg <- belt[s:e]; dur <- (e - s) / fs
    ent <- b_i$bids_entity[j]

    rec <- resp_recording(seg, fs = fs, units = "V",
                          id = paste(sb, ent, sep = "_"))

    # Amplitude thresholds disabled. With iqr_* = NA the amplitude tests can no
    # longer return a verdict, so `quality` here reflects only the frequency,
    # flatline, saturation and band-ratio tests -- which BeforeYouRun says do
    # transfer as-is.
    q <- try(resp_quality(rec, trim_s = TRIM_S,
                          thresholds = list(iqr_unusable = NA, iqr_degraded = NA)),
             silent = TRUE)
    if (inherits(q, "try-error")) next

    # Detection, for the rate ceiling sanity check (step 3) and the duty cycle
    # that carries the polarity verdict (step 2). quality = FALSE because we
    # just ran it with the thresholds we actually want.
    a <- try(resp_analyse(rec, quality = FALSE), silent = TRUE)
    bt <- if (inherits(a, "try-error")) NULL else a$breaths

    dc <- if (!is.null(bt) && nrow(bt)) median(bt$duty_cycle, na.rm = TRUE) else NA_real_
    n_keep <- if (!is.null(bt) && "keep" %in% names(bt)) sum(bt$keep %in% TRUE) else NA_integer_

    rows[[length(rows) + 1]] <- tibble(
      subject = sb, bids_entity = ent, part = bp, dur_s = dur,
      task = str_match(ent, "task-([a-z]+)")[, 2],
      filtered_iqr      = q$filtered_iqr,
      band_ratio        = q$band_ratio,
      dominant_rate_bpm = q$dominant_rate_bpm,
      flatline_pct      = q$flatline_pct,
      saturation_pct    = q$saturation_pct,
      na_pct            = q$na_pct,
      quality_nonamp    = q$quality,
      n_breaths         = if (is.null(bt)) NA_integer_ else nrow(bt),
      n_breaths_keep    = n_keep,
      rate_bt_brpm      = if (is.null(bt) || !nrow(bt)) NA_real_ else 60 * nrow(bt) / dur,
      duty_cycle_med    = dc,
      flags             = paste(q$flags, collapse = "; ")
    )
  }
  cat(sprintf("  %-9s runs=%d\n", sb, sum(map_chr(rows, "subject") == sb)))
}

cal <- bind_rows(rows)
write_csv(cal, file.path(dir_out, "resp_calibration.csv"))

cat(sprintf("\nruns measured: %d   participants: %d\n",
            nrow(cal), n_distinct(cal$subject)))

# ---- step 1: the amplitude cut ---------------------------------------------
cat("\n== filtered_iqr distribution (Volts) ==\n")
print(round(quantile(cal$filtered_iqr, c(0, .01, .05, .1, .25, .5, .75, .9, 1),
                     na.rm = TRUE), 5))
cat("\nby task:\n")
print(as.data.frame(cal %>% group_by(task) %>%
  summarise(n = n(), p05 = quantile(filtered_iqr, .05, na.rm = TRUE),
            med = median(filtered_iqr, na.rm = TRUE),
            p95 = quantile(filtered_iqr, .95, na.rm = TRUE), .groups = "drop")))

# ---- step 2: polarity, spontaneous runs only -------------------------------
cat("\n== median duty cycle by task (>0.5 under 'up' polarity means inverted) ==\n")
print(as.data.frame(cal %>% group_by(task) %>%
  summarise(n = sum(is.finite(duty_cycle_med)),
            duty_med = median(duty_cycle_med, na.rm = TRUE),
            frac_over_half = mean(duty_cycle_med > 0.5, na.rm = TRUE),
            .groups = "drop")))
spont <- cal %>% filter(task %in% c("localizer", "authenticity"))
cat(sprintf("\nspontaneous runs (localizer+authenticity): n=%d  median duty=%.3f\n",
            nrow(spont), median(spont$duty_cycle_med, na.rm = TRUE)))
print(polarity_verdict(median(spont$duty_cycle_med, na.rm = TRUE)))

# ---- step 3: rate ceiling vs paradigm --------------------------------------
cat("\n== breath rate (brpm) vs paradigm ==\n")
print(as.data.frame(cal %>% group_by(task) %>%
  summarise(n = n(), p05 = quantile(rate_bt_brpm, .05, na.rm = TRUE),
            med = median(rate_bt_brpm, na.rm = TRUE),
            p95 = quantile(rate_bt_brpm, .95, na.rm = TRUE),
            max = max(rate_bt_brpm, na.rm = TRUE), .groups = "drop")))

png(file.path(dir_fig, "resp_calibration.png"), width = 1400, height = 900, res = 110)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
hist(log10(cal$filtered_iqr), breaks = 40, col = "grey80", border = "white",
     main = "filtered_iqr (log10 Volts)", xlab = "log10 IQR")
hist(cal$duty_cycle_med, breaks = 40, col = "grey80", border = "white",
     main = "median duty cycle", xlab = "duty cycle"); abline(v = 0.5, col = "red")
hist(cal$rate_bt_brpm, breaks = 40, col = "grey80", border = "white",
     main = "breath rate", xlab = "brpm")
plot(log10(cal$filtered_iqr), cal$rate_bt_brpm, pch = 16, col = "#00000060",
     main = "amplitude vs rate", xlab = "log10 IQR", ylab = "brpm")
par(op); dev.off()

cat(sprintf("\nRun ended: %s\n", Sys.time()))
