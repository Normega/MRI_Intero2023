# intero2023_resp.R
# Per-epoch breathing rate from the belt, pacer adherence, and the
# physiological ground truth of rate change.
#
# Produces:
#   Repo/Data/resp_epochs.csv    one row per breathing epoch, measured + instructed
#   Repo/Data/resp_runs.csv      one row per run, belt gate and coverage
#   Repo/Data/qc_ledger.csv      appended with belt_usable / resp_epoch_coverage
#   <cache>/sub-XXX_resp.rds     breath table per run, cleaned
#
# It does NOT rewrite the event TSVs. That merge is a separate step
# (intero2023_events_physio.R) so the classification thresholds can be
# inspected and revised without regenerating the events themselves.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root  <- "I:/Shared drives/MRI_Intero2023"
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_ev    <- file.path(dir_out, "events")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")

rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R", "intero_resp.R")) source(file.path(dir_R, f))

cat(sprintf("Run started: %s\n", Sys.time()))

blocks <- read_csv(file.path(dir_out, "physio_blocks.csv"),
                   show_col_types = FALSE, progress = FALSE) %>%
  filter(kind == "bold", !is.na(bids_entity), is.na(dup_slot) | !dup_slot,
         grepl("task-breathing", bids_entity))

cache_files <- list.files(dir_cache, pattern = "_physio[.]rds$", full.names = TRUE)
cache_key <- tibble(path = cache_files,
                    subject = sub("(_part-[^_]*)?_physio[.]rds$", "", basename(cache_files)),
                    part = ifelse(grepl("_part-", basename(cache_files)),
                                  sub("^.*_part-([^_]*)_physio[.]rds$", "\\1",
                                      basename(cache_files)), "-"))

ep_all <- list(); run_all <- list(); qc <- list()

for (sb in sort(unique(blocks$subject))) {
  b_i <- blocks %>% filter(subject == sb)
  cf  <- cache_key %>% filter(subject == sb)
  if (!nrow(cf)) next
  bts <- list()

  for (j in seq_len(nrow(b_i))) {
    ent <- b_i$bids_entity[j]
    bp <- b_i$part[j]; if (is.na(bp)) bp <- "-"
    cfp <- cf %>% filter(part == bp); if (!nrow(cfp)) cfp <- cf[1, ]
    r <- try(readRDS(cfp$path[1]), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r$belt)) next

    fs <- r$meta$fs; belt <- as.numeric(r$belt$signal)
    s <- max(1L, round(b_i$onset_s[j] * fs))
    e <- min(length(belt), round(b_i$end_s[j] * fs))
    if (e - s < fs * 60) next
    a <- belt_analyse(belt[s:e], fs, id = paste(sb, ent, sep = "_"))
    if (is.null(a)) next
    grade <- belt_grade(a$quality)

    evf <- file.path(dir_ev, sprintf("%s_%s_events.tsv", sb, ent))
    ptf <- file.path(dir_ev, sprintf("%s_%s_pacertrace.tsv", sb, ent))
    if (!file.exists(evf)) next
    ev <- suppressWarnings(read_tsv(evf, show_col_types = FALSE, progress = FALSE,
                                    col_types = cols(.default = col_character()))) %>%
      mutate(onset = as.numeric(onset), duration = as.numeric(duration),
             epoch_index = suppressWarnings(as.integer(epoch_index))) %>%
      filter(trial_type %in% c("paced_breathing", "free_breathing"))
    if (!nrow(ev)) next
    ins <- if (file.exists(ptf))
      pacer_instructed(suppressWarnings(read_tsv(ptf, show_col_types = FALSE,
        progress = FALSE, col_types = cols(.default = col_character())))) else NULL

    er <- epoch_resp(a$bt, ev, ins)
    ad <- run_adherence(a$bt, ev, ins)
    er <- er %>% mutate(subject = sb, bids_entity = ent, belt_grade = grade,
                        adhere_r = ad$adhere_r, adhere_gain = ad$adhere_gain,
                        adhere_mae_run_s = ad$adhere_mae_s,
                        period_ratio_oct = ad$period_ratio_oct,
                        .before = 1) %>%
      left_join(ev %>% select(epoch_index, trial_type, change_direction,
                              salience, response, sdt_class) %>%
                  mutate(change_amount = NA_character_) %>% select(-change_amount),
                by = c("epoch_index", "trial_type")) %>%
      mutate(measured = resp_measured(n_breaths, coverage_frac, grade))
    ep_all[[length(ep_all) + 1]] <- er
    bts[[ent]] <- strip_phase(a$bt)

    run_all[[length(run_all) + 1]] <- tibble(
      subject = sb, bids_entity = ent, part = bp, belt_grade = grade,
      filtered_iqr = a$quality$filtered_iqr %||% NA_real_,
      band_ratio   = a$quality$band_ratio %||% NA_real_,
      n_breaths_run = nrow(a$bt),
      n_epochs = nrow(er), n_measured = sum(er$measured, na.rm = TRUE)) %>%
      bind_cols(ad)

    qc[[length(qc) + 1]] <- qc_entry(sb, "physio", "pacer_adherence",
      if (!is.finite(ad$adhere_r)) "fail"
      else if (ad$adhere_r >= 0.5) "good"
      else if (ad$adhere_r >= ADHERE_MIN_R) "info" else "degraded",
      entity = ent, value = ad$adhere_r,
      note = sprintf("gain=%.2f mae=%.2fs n=%d period_ratio=%.2f",
                     ad$adhere_gain, ad$adhere_mae_s, ad$n_paired, ad$period_ratio_oct),
      source = "measured vs instructed breath period, pooled over paced epochs")

    qc[[length(qc) + 1]] <- qc_entry(sb, "physio", "belt_usable",
      if (grade == "good") "good" else if (grade == "degraded") "degraded" else "fail",
      entity = ent, value = a$quality$filtered_iqr %||% NA_real_,
      note = sprintf("filtered_iqr=%.5f band_ratio=%.3f",
                     a$quality$filtered_iqr %||% NA_real_, a$quality$band_ratio %||% NA_real_),
      source = "resp_quality, cut calibrated on this sample")
    qc[[length(qc) + 1]] <- qc_entry(sb, "physio", "resp_epoch_coverage",
      if (nrow(er) == 0) "fail"
      else if (mean(er$measured, na.rm = TRUE) > 0.8) "good"
      else if (mean(er$measured, na.rm = TRUE) > 0.5) "degraded" else "fail",
      entity = ent, value = mean(er$measured, na.rm = TRUE),
      note = sprintf("%d of %d epochs measurable", sum(er$measured, na.rm = TRUE), nrow(er)),
      source = "epoch_features coverage + breath count")
  }

  if (length(bts))
    saveRDS(list(subject = sb, breaths = bts, built = as.character(Sys.Date()),
                 note = "no phase-timing columns: belt polarity unresolved"),
            file.path(dir_cache, paste0(sb, "_resp.rds")), compress = "xz")
  cat(sprintf("  %-9s runs=%d\n", sb, length(bts)))
}

ep <- bind_rows(ep_all); rn <- bind_rows(run_all)

# ---- the No Change null, which sets the classification thresholds ----------
nc <- ep %>% filter(trial_type == "paced_breathing", change_direction == "No Change",
                    measured, is.finite(actual_ratio))
cat(sprintf("\n== No Change null: n=%d measured epochs ==\n", nrow(nc)))
print(round(quantile(nc$actual_ratio, c(.05,.25,.5,.75,.95), na.rm = TRUE), 3))
lo <- unname(quantile(nc$actual_ratio, .05, na.rm = TRUE))
hi <- unname(quantile(nc$actual_ratio, .95, na.rm = TRUE))
cat(sprintf("thresholds: Faster if ratio < %.3f, Slower if > %.3f\n", lo, hi))

ep$actual_change_physio <- classify_actual_change(ep$actual_ratio, ep$measured,
                                                  lo = lo, hi = hi)

# adherence_verified is RUN-level: it says the participant demonstrably
# followed the pacer during this run, so the instructed change can be trusted
# as the stimulus for this epoch. It is NOT "this epoch's measured change
# matched its instructed change" -- that per-epoch agreement is only estimable
# at AUC ~0.68 and would encode measurement noise as behaviour.
ep$adherence_verified <- ep$measured &
  is.finite(ep$adhere_r) & ep$adhere_r >= ADHERE_MIN_R &
  ep$belt_grade %in% c("good", "degraded")

write_csv(ep, file.path(dir_out, "resp_epochs.csv"))
write_csv(rn, file.path(dir_out, "resp_runs.csv"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), bind_rows(qc)),
         file.path(dir_out, "qc_ledger.csv"))

cat(sprintf("\nepochs=%d  runs=%d  participants=%d\n",
            nrow(ep), nrow(rn), n_distinct(ep$subject)))
cat("\n-- belt grade by run --\n"); print(table(rn$belt_grade))
cat("\n-- epochs measurable --\n"); print(table(ep$trial_type, ep$measured))
cat("\n-- instructed vs actual, paced measured epochs --\n")
pm <- ep %>% filter(trial_type == "paced_breathing", measured)
print(as.data.frame(pm %>% group_by(change_direction, salience) %>% summarise(
  n = n(), instructed = round(median(instructed_ratio, na.rm = TRUE), 3),
  actual = round(median(actual_ratio, na.rm = TRUE), 3),
  mae_s = round(median(adherence_mae_s, na.rm = TRUE), 2),
  r = round(median(adherence_r, na.rm = TRUE), 2), .groups = "drop")))
cat("\n-- run-level pacer adherence by belt grade --\n")
# NB: the summary columns must not reuse a source column's name. summarise()
# evaluates its arguments in order against a mask that already holds the
# earlier results, so naming the median `adhere_r` would make the next
# expression read the median rather than the column, and frac_tracking would
# collapse to 0 or 1.
print(as.data.frame(rn %>% group_by(belt_grade) %>% summarise(
  n = n(), n_finite = sum(is.finite(adhere_r)),
  adhere_r_med = round(median(adhere_r, na.rm = TRUE), 3),
  period_ratio_med = round(median(period_ratio_oct, na.rm = TRUE), 3),
  frac_tracking = round(mean(adhere_r >= ADHERE_MIN_R, na.rm = TRUE), 2),
  .groups = "drop")))

usable <- rn %>% filter(belt_grade %in% c("good", "degraded"),
                        is.finite(adhere_r), adhere_r >= ADHERE_MIN_R)
ps <- usable %>% count(subject)
cat(sprintf("\nRUNS with a working belt AND verified pacer tracking: %d of %d\n",
            nrow(usable), nrow(rn)))
cat(sprintf("PARTICIPANTS: both breathing runs %d, at least one %d\n",
            sum(ps$n == 2), nrow(ps)))
cat(sprintf("\nepochs with adherence_verified: %d of %d paced (%.0f%%)\n",
            sum(ep$adherence_verified & ep$trial_type == "paced_breathing", na.rm = TRUE),
            sum(ep$trial_type == "paced_breathing"),
            100 * mean(ep$adherence_verified[ep$trial_type == "paced_breathing"], na.rm = TRUE)))
cat(sprintf("\nRun ended: %s\n", Sys.time()))
