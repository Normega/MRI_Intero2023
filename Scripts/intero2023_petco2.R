# intero2023_petco2.R
# Extracts PetCO2 per breath and per volume for every canonical BOLD run,
# validates the capnogram over the whole run, and caches cleaned RDS files.
#
# Produces:
#   <cache>/sub-XXX_petco2.rds     per-breath series + per-volume regressor + QC
#   Repo/Data/petco2_runs.csv      one row per run: quality and validity
#   Repo/Data/qc_ledger.csv        appended with capnogram_valid / petco2_usable
#
# The per-volume series is the GLM regressor: linear interpolation between
# consecutive end-tidal values, sampled at TR. Volumes outside the detected
# range stay NA rather than being extrapolated, so a capnogram that dies partway
# shows a gap instead of a flat line.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root  <- "I:/Shared drives/MRI_Intero2023"
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")
TR <- 2

rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_logs.R", "intero_physio.R",
            "intero_petco2.R", "intero_qc.R")) source(file.path(dir_R, f))

cat(sprintf("Run started: %s\n", Sys.time()))

blocks <- read_csv(file.path(dir_out, "physio_blocks.csv"),
                   show_col_types = FALSE, progress = FALSE) %>%
  filter(kind == "bold", !is.na(bids_entity), is.na(dup_slot) | !dup_slot)

# cache filenames come in two vintages; normalise to a subject key
cache_files <- list.files(dir_cache, pattern = "_physio[.]rds$", full.names = TRUE)
cache_key <- tibble(path = cache_files,
                    subject = sub("(_part-[^_]*)?_physio[.]rds$", "",
                                  basename(cache_files)),
                    part = ifelse(grepl("_part-", basename(cache_files)),
                                  sub("^.*_part-([^_]*)_physio[.]rds$", "\\1",
                                      basename(cache_files)), "-"))

runs <- list(); qc <- list()

for (sb in sort(unique(blocks$subject))) {
  b_i <- blocks %>% filter(subject == sb)
  cf  <- cache_key %>% filter(subject == sb)
  if (!nrow(cf)) next
  per_run <- list(); series <- list()

  for (j in seq_len(nrow(b_i))) {
    bp <- b_i$part[j]; if (is.na(bp)) bp <- "-"
    cfp <- cf %>% filter(part == bp)
    if (!nrow(cfp)) cfp <- cf[1, ]
    r <- try(readRDS(cfp$path[1]), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r$co2)) next
    fs  <- r$meta$fs
    co2 <- as.numeric(r$co2$signal)
    s <- max(1L, round(b_i$onset_s[j] * fs))
    e <- min(length(co2), round(b_i$end_s[j] * fs))
    if (e - s < fs * 30) next
    seg <- co2[s:e]; dur <- (e - s) / fs

    pk  <- petco2_flag(petco2_detect(seg, fs))
    val <- capnogram_validity(pk, dur)
    q   <- petco2_quality(pk, dur)
    nvol <- b_i$n_volumes[j]
    if (is.na(nvol)) nvol <- floor(dur / TR)
    vols <- petco2_to_volumes(pk, nvol, tr = TR)

    ent <- b_i$bids_entity[j]
    per_run[[length(per_run)+1]] <- bind_cols(
      tibble(subject = sb, bids_entity = ent, part = bp,
             dur_s = dur, n_volumes = nvol,
             n_vol_covered = sum(!is.na(vols))), q,
      val %>% select(-any_of(c("baseline_med", "rate_brpm"))))
    series[[ent]] <- list(breaths = pk, volumes = vols, tr = TR,
                          onset_s = b_i$onset_s[j], fs = fs)

    qc[[length(qc)+1]] <- qc_entry(sb, "physio", "capnogram_valid",
      if (isTRUE(val$capnogram_valid)) "good" else "fail",
      entity = ent, value = val$baseline_med,
      note = if (is.na(val$fail_reason))
               sprintf("baseline=%.2f swing=%.2f rate=%.1f", val$baseline_med,
                       val$swing_med, val$rate_brpm) else val$fail_reason,
      source = "whole-run capnogram morphology")
    qc[[length(qc)+1]] <- qc_entry(sb, "physio", "petco2_coverage",
      if (!nvol) "fail" else if (sum(!is.na(vols))/nvol > 0.9) "good"
      else if (sum(!is.na(vols))/nvol > 0.5) "degraded" else "fail",
      entity = ent, value = sum(!is.na(vols))/max(nvol,1),
      note = sprintf("%d of %d volumes have PetCO2", sum(!is.na(vols)), nvol),
      source = "interpolated end-tidal series")
  }

  if (!length(per_run)) next
  pr <- bind_rows(per_run)
  runs[[length(runs)+1]] <- pr
  saveRDS(list(subject = sb, runs = pr, series = series,
               built = as.character(Sys.Date()),
               units = "%CO2 (x~7.1 for mmHg at sea level)"),
          file.path(dir_cache, paste0(sb, "_petco2.rds")), compress = "xz")
  cat(sprintf("  %-9s runs=%d valid=%d\n", sb, nrow(pr), sum(pr$capnogram_valid)))
}

rt <- bind_rows(runs)
write_csv(rt, file.path(dir_out, "petco2_runs.csv"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), bind_rows(qc)),
         file.path(dir_out, "qc_ledger.csv"))

cat(sprintf("\nruns processed: %d   capnogram_valid: %d (%.0f%%)\n",
            nrow(rt), sum(rt$capnogram_valid), 100*mean(rt$capnogram_valid)))
cat("\n-- valid runs per participant --\n")
ps <- rt %>% group_by(subject) %>%
  summarise(n = n(), valid = sum(capnogram_valid),
            breathing_valid = sum(capnogram_valid & grepl("breathing", bids_entity)),
            .groups = "drop")
cat(sprintf("participants with >=1 valid run : %d\n", sum(ps$valid > 0)))
cat(sprintf("participants with BOTH breathing runs valid: %d\n",
            sum(ps$breathing_valid == 2)))
cat("\n-- fail reasons --\n")
print(as.data.frame(rt %>% filter(!capnogram_valid) %>% count(fail_reason, sort = TRUE)))
cat(sprintf("\nRun ended: %s\n", Sys.time()))
