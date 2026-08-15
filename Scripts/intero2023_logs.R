# intero2023_logs.R
# Builds the log-derived run inventory and QC gate for Intero2023 / breathMRI2026.
#
# Produces (in Repo/Data/):
#   log_run_attempts.csv  every run attempt found in every log file, nothing dropped
#   log_runs.csv          the run of record for each participant/session/run
#   log_sessions.csv      one row per participant-session
#   log_qc_flags.csv      reconciliation against the imaging inventory
#
# This is the gate the physio and imaging work depends on: a run is analysable
# only if the log gives it a scanner trigger to anchor onsets to.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(purrr); library(tibble); library(tidyr)
})

# ---- paths ----------------------------------------------------------------

dir_root  <- "I:/Shared drives/MRI_Intero2023"
dir_logs  <- file.path(dir_root, "Data", "Logs")
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_R     <- file.path(dir_root, "Repo", "R")
f_invent  <- file.path(dir_root, "Prompts", "inventory.tsv")

source(file.path(dir_R, "intero_ids.R"))
source(file.path(dir_R, "intero_logs.R"))

dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Run started: %s\n", Sys.time()))

# ---- 1. scan --------------------------------------------------------------

scan <- scan_log_dir(dir_logs)
files <- scan$files
attempts <- scan$attempts

cat(sprintf("Log files scanned      : %d\n", nrow(files)))
cat(sprintf("  byte-identical dupes : %d\n", sum(files$is_duplicate)))
cat(sprintf("  empty / no data rows : %d\n",
            sum(!files$file %in% attempts$file)))
cat(sprintf("  pilot/QA files       : %d\n", sum(files$is_pilot)))
cat(sprintf("Run attempts found     : %d (%d scanner runs)\n",
            nrow(attempts), sum(attempts$is_scanner_run)))

non_std <- attempts %>% filter(!is_scanner_run) %>% count(log_run)
if (nrow(non_std)) {
  cat("Non-scanner runs (excluded from the 5-run protocol):\n")
  for (i in seq_len(nrow(non_std)))
    cat(sprintf("  %-18s %d\n", non_std$log_run[i], non_std$n[i]))
}

# ---- 2. select the run of record ------------------------------------------

# Duplicate files would otherwise inflate the attempt count for a run.
selected <- attempts %>%
  filter(!is_duplicate) %>%
  select_runs()

runs <- selected %>% filter(is_selected)

cat(sprintf("\nRuns of record         : %d\n", nrow(runs)))
cat(sprintf("  with scanner trigger : %d\n", sum(runs$has_trigger)))
cat(sprintf("  WITHOUT trigger      : %d  <- not alignable\n",
            sum(!runs$has_trigger)))
cat(sprintf("  superseded attempts  : %d\n", sum(!selected$is_selected)))

sessions <- summarise_sessions(selected)
cat(sprintf("Participant-sessions   : %d (%d real, %d pilot)\n",
            nrow(sessions), sum(!sessions$is_pilot), sum(sessions$is_pilot)))

multi <- sessions %>% filter(!is_pilot) %>% count(pid) %>% filter(n > 1)
if (nrow(multi)) {
  cat("Participants with >1 session date:\n")
  for (i in seq_len(nrow(multi)))
    cat(sprintf("  %s (%d dates)\n", pid_lab(multi$pid[i]), multi$n[i]))
}

# ---- 3. reconcile against the imaging inventory ---------------------------

# inventory.tsv is the Linux box report: the dataset of record. Volume counts
# come from the NIfTI headers, so run duration x TR is the independent check
# on the log timing model.
TR <- 2.0

qc <- tibble()
if (file.exists(f_invent)) {
  inv <- read_tsv(f_invent, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      pid           = pid_num(.data$subject),
      protocol_slot = as.integer(.data$protocol_run),
      img_dur_s     = .data$n_volumes * TR
    ) %>%
    select(pid, protocol_slot, n_volumes, img_dur_s, fmriprep_status,
           sdc, fd_mean, fd_max, `fd_gt0.5`)

  qc <- runs %>%
    filter(!is_pilot) %>%
    select(pid, subject, session_date, log_run, protocol_slot, bids_entity,
           protocol_version, has_trigger, completed, n_trials, dur_s,
           resp_rate, file) %>%
    full_join(inv, by = c("pid", "protocol_slot")) %>%
    mutate(
      subject     = coalesce(subject, pid_bids(pid)),
      bids_entity = coalesce(bids_entity,
                             bids_entities(INTERO_RUNS$log_run[
                               match(protocol_slot, INTERO_RUNS$protocol_slot)])),
      has_log     = !is.na(file),
      has_imaging = !is.na(n_volumes),
      # Log duration is measured from the first trigger to the last logged
      # event; the acquisition keeps running a little past that, so the log is
      # expected to be slightly SHORTER. A large gap either way is a flag.
      dur_diff_s  = dur_s - img_dur_s
    ) %>%
    mutate(
      flag = case_when(
        has_log & !has_imaging & has_trigger ~ "log_only_triggered",
        has_log & !has_imaging               ~ "log_only",
        !has_log & has_imaging               ~ "imaging_only",
        has_log & has_imaging & !has_trigger ~ "no_trigger_not_alignable",
        !is.na(dur_diff_s) & dur_diff_s > 5  ~ "log_longer_than_acquisition",
        !is.na(dur_diff_s) & dur_diff_s < -60 ~ "log_much_shorter",
        TRUE                                 ~ "ok"
      )
    ) %>%
    arrange(pid, protocol_slot)

  cat("\n-- reconciliation vs imaging inventory --\n")
  tally <- qc %>% count(flag, sort = TRUE)
  for (i in seq_len(nrow(tally)))
    cat(sprintf("  %-28s %d\n", tally$flag[i], tally$n[i]))

  bad <- qc %>% filter(flag != "ok")
  if (nrow(bad)) {
    cat("\n  flagged runs:\n")
    for (i in seq_len(nrow(bad))) {
      cat(sprintf("    %-9s %-24s %-28s log_dur=%-8s vols=%s\n",
                  bad$subject[i], bad$bids_entity[i], bad$flag[i],
                  ifelse(is.na(bad$dur_s[i]), "-", sprintf("%.0fs", bad$dur_s[i])),
                  ifelse(is.na(bad$n_volumes[i]), "-", bad$n_volumes[i])))
    }
  }

  ok <- qc %>% filter(flag == "ok", !is.na(dur_diff_s))
  if (nrow(ok)) {
    cat(sprintf("\n  log-vs-acquisition duration difference (n=%d matched runs):\n",
                nrow(ok)))
    cat(sprintf("    median %.1fs   IQR [%.1f, %.1f]   range [%.1f, %.1f]\n",
                median(ok$dur_diff_s),
                quantile(ok$dur_diff_s, .25), quantile(ok$dur_diff_s, .75),
                min(ok$dur_diff_s), max(ok$dur_diff_s)))
  }
} else {
  cat(sprintf("\n! inventory not found at %s -- skipping reconciliation\n", f_invent))
}

# ---- 4. write -------------------------------------------------------------

stamp <- function(df) {
  df %>% mutate(log_script_version = "1.0",
                log_date_built     = as.character(Sys.Date()))
}

write_csv(stamp(selected), file.path(dir_out, "log_run_attempts.csv"))
write_csv(stamp(runs),     file.path(dir_out, "log_runs.csv"))
write_csv(stamp(sessions), file.path(dir_out, "log_sessions.csv"))
if (nrow(qc)) write_csv(stamp(qc), file.path(dir_out, "log_qc_flags.csv"))

cat("\n-- written --\n")
for (f in c("log_run_attempts.csv", "log_runs.csv", "log_sessions.csv",
            "log_qc_flags.csv")) {
  p <- file.path(dir_out, f)
  if (file.exists(p)) cat(sprintf("  %-24s %d rows\n", f, nrow(read_csv(p, show_col_types = FALSE, progress = FALSE))))
}
cat(sprintf("Run ended: %s\n", Sys.time()))
