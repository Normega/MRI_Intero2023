# intero2023_events.R
# Generates BIDS events files + the pacer trace + a behavioural trial table
# for the Intero2023 breathing runs.
#
# Produces:
#   Repo/Data/events/sub-XXX_task-breathing_run-N_events.tsv   BIDS, per run
#   Repo/Data/events/sub-XXX_task-breathing_run-N_pacertrace.tsv
#   Repo/Data/breathing_trials.csv    one row per breathing epoch, all subjects
#   Repo/Data/events_coverage.csv     what was produced and what was refused
#
# Naming follows the LINUX BOX copy (task-breathing_run-1), which is the
# dataset of record. The Google Drive Data/BIDS copy uses task-breath1, has no
# fieldmaps, and has shifted labels for sub-082/104 -- do not pair with it.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(purrr); library(tibble)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_logs <- file.path(dir_root, "Data", "Logs")
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_ev   <- file.path(dir_out, "events")
dir_R    <- file.path(dir_root, "Repo", "R")

for (f in c("intero_ids.R", "intero_logs.R", "intero_events.R"))
  source(file.path(dir_R, f))

dir.create(dir_ev, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Run started: %s\n", Sys.time()))

runs <- read_csv(file.path(dir_out, "log_runs.csv"),
                 show_col_types = FALSE, progress = FALSE) %>%
  filter(!is_pilot, has_trigger, log_run %in% c("ChangeDetection1", "ChangeDetection2"))

cat(sprintf("Breathing runs with a trigger: %d\n", nrow(runs)))

log_cache <- new.env(parent = emptyenv())
get_log <- function(fn) {
  if (!exists(fn, envir = log_cache))
    assign(fn, read_log_file(file.path(dir_logs, fn)), envir = log_cache)
  get(fn, envir = log_cache)
}

cov <- list(); trials <- list()

for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]
  d <- get_log(r$file) %>% filter(Run == r$log_run)
  res <- extract_breathing_run(d)

  n_ev <- nrow(res$events)
  # Cycle count is NOT fixed by protocol version: P51 ran 5 per breathing run,
  # P053-P059 ran 7, and v2 runs 7. So completeness is judged by balance --
  # every paced trial should be followed by a free period -- rather than
  # against a hardcoded expectation that would misreport whole cohorts.
  n_p <- if (n_ev) sum(res$events$trial_type == "paced_breathing") else 0L
  n_f <- if (n_ev) sum(res$events$trial_type == "free_breathing")  else 0L
  status <- if (n_ev == 0) "refused_unparsed"
            else if (n_p == n_f && n_p >= 5L) "complete"
            else "partial"

  cov[[length(cov) + 1]] <- tibble(
    subject = r$subject, bids_entity = r$bids_entity,
    protocol_version = r$protocol_version, file = r$file,
    n_events = n_ev,
    n_paced = n_p, n_free = n_f, n_cycles = n_p,
    n_pacer = nrow(res$pacer), status = status
  )

  if (n_ev == 0) next

  base <- paste0(r$subject, "_", r$bids_entity)
  # BIDS wants onset/duration/trial_type first and n/a for missing values.
  ev_out <- res$events %>%
    relocate(onset, duration, trial_type) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "n/a", as.character(.x))))
  write_tsv(ev_out, file.path(dir_ev, paste0(base, "_events.tsv")), na = "n/a")
  write_tsv(res$pacer, file.path(dir_ev, paste0(base, "_pacertrace.tsv")), na = "n/a")

  trials[[length(trials) + 1]] <- res$events %>%
    filter(trial_type %in% c("paced_breathing", "free_breathing")) %>%
    mutate(subject = r$subject, bids_entity = r$bids_entity,
           protocol_version = r$protocol_version, .before = 1)
}

coverage <- bind_rows(cov)
trial_tbl <- bind_rows(trials)

write_csv(coverage, file.path(dir_out, "events_coverage.csv"))
write_csv(trial_tbl, file.path(dir_out, "breathing_trials.csv"))

cat("\n-- coverage --\n")
cov_tally <- count(coverage, status, sort = TRUE)
for (i in seq_len(nrow(cov_tally)))
  cat(sprintf("  %-28s %d\n", cov_tally$status[i], cov_tally$n[i]))

part <- coverage %>% filter(status == "partial")
if (nrow(part)) {
  cat("\n  partial runs:\n")
  for (i in seq_len(nrow(part)))
    cat(sprintf("    %-9s %-24s paced=%d free=%d\n",
                part$subject[i], part$bids_entity[i], part$n_paced[i], part$n_free[i]))
}

cat(sprintf("\nBreathing epochs written: %d (%d paced, %d free)\n",
            nrow(trial_tbl),
            sum(trial_tbl$trial_type == "paced_breathing"),
            sum(trial_tbl$trial_type == "free_breathing")))

paced <- trial_tbl %>% filter(trial_type == "paced_breathing")
cat("\n-- signal detection, paced trials --\n")
sdt <- paced %>% count(sdt_class, sort = TRUE)
for (i in seq_len(nrow(sdt)))
  cat(sprintf("  %-20s %4d (%.1f%%)\n", sdt$sdt_class[i], sdt$n[i],
              100 * sdt$n[i] / nrow(paced)))

cat("\n-- detection by salience (change trials only) --\n")
det <- paced %>%
  filter(change_direction != "No Change", !is.na(salience)) %>%
  group_by(salience) %>%
  summarise(n = n(), hit_rate = mean(sdt_class == "hit"),
            miss_rate = mean(sdt_class == "miss"),
            mean_conf = mean(confidence, na.rm = TRUE), .groups = "drop")
print(as.data.frame(det))

cat("\n-- confidence by outcome (the awareness axis) --\n")
print(as.data.frame(
  paced %>% group_by(sdt_class) %>%
    summarise(n = n(), mean_conf = mean(confidence, na.rm = TRUE),
              mean_arousal = mean(arousal, na.rm = TRUE), .groups = "drop")))

cat(sprintf("\nRun ended: %s\n", Sys.time()))
