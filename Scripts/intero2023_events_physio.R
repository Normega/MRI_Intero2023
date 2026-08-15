# intero2023_events_physio.R
# Merges the belt-derived respiration measures into the event files, filling
# the two columns intero_events.R reserved as NA.
#
# Kept separate from intero2023_resp.R so the thresholds behind the merged
# columns can be revised without re-reading any physiology.
#
# WHAT GOES IN, AND HOW MUCH TO TRUST IT
#
#   adherence_verified   RUN-level. The participant demonstrably tracked the
#                        pacer during this run, so the INSTRUCTED change is a
#                        valid stimulus for this epoch and sdt_class computed
#                        against it means what it says. Well estimated: pooled
#                        over ~23 paired breaths per run.
#   actual_period_ratio  EPOCH-level, continuous. Measured breath period at
#                        the end of the epoch relative to its start. This is
#                        the recommended regressor.
#   actual_change_physio EPOCH-level, categorical, and DELIBERATELY
#                        CONSERVATIVE. Thresholded at the 5th/95th percentile
#                        of the No Change null, so it fires only on the
#                        largest changes. An epoch holds 3-4 breaths and five
#                        different estimators all separate Faster from Slower
#                        at AUC 0.66-0.69, so this column cannot carry a
#                        per-trial analysis on its own -- prefer
#                        actual_period_ratio, or condition on
#                        adherence_verified and use the instructed change.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_ev   <- file.path(dir_out, "events")
dir_R    <- file.path(dir_root, "Repo", "R")
for (f in c("intero_ids.R", "intero_qc.R")) source(file.path(dir_R, f))

cat(sprintf("Run started: %s\n", Sys.time()))

ep <- read_csv(file.path(dir_out, "resp_epochs.csv"), show_col_types = FALSE,
               progress = FALSE) %>%
  transmute(subject, bids_entity, epoch_index, trial_type,
            resp_measured        = measured,
            resp_n_breaths       = n_breaths,
            resp_coverage        = round(coverage_frac, 3),
            resp_rate_brpm       = round(rate_med_brpm, 2),
            actual_period_ratio  = round(actual_ratio, 4),
            instructed_period_ratio = round(instructed_ratio, 4),
            adhere_r             = round(adhere_r, 3),
            belt_grade,
            adherence_verified,
            actual_change_physio)

# ---- event TSVs ------------------------------------------------------------
n_files <- 0L; n_filled <- 0L
for (f in list.files(dir_ev, pattern = "_events[.]tsv$", full.names = TRUE)) {
  bn  <- basename(f)
  sb  <- str_match(bn, "^(sub-[0-9]+)_")[, 2]
  ent <- str_match(bn, "^sub-[0-9]+_(task-[a-z]+(?:_run-[0-9]+)?)_events")[, 2]
  if (is.na(sb) || is.na(ent)) next
  e_i <- ep %>% filter(subject == sb, bids_entity == ent)
  if (!nrow(e_i)) next

  d <- suppressWarnings(read_tsv(f, show_col_types = FALSE, progress = FALSE,
                                 col_types = cols(.default = col_character())))
  d$.k <- suppressWarnings(as.integer(d$epoch_index))

  # Join on (epoch_index, trial_type): epoch_index alone is not unique because
  # a paced and a free epoch share an index within each cycle.
  j <- e_i %>% select(-subject, -bids_entity) %>% rename(.k = epoch_index)
  d2 <- d %>% select(-any_of(c("adherence_verified", "actual_change_physio"))) %>%
    left_join(j, by = c(".k", "trial_type"))

  # Rating and valence rows are not breathing epochs and must stay n/a rather
  # than inherit a neighbouring epoch's physiology.
  is_ep <- d2$trial_type %in% c("paced_breathing", "free_breathing")
  newc <- c("resp_measured", "resp_n_breaths", "resp_coverage", "resp_rate_brpm",
            "actual_period_ratio", "instructed_period_ratio", "adhere_r",
            "belt_grade", "adherence_verified", "actual_change_physio")
  for (cc in newc) {
    v <- as.character(d2[[cc]])
    v[!is_ep | is.na(v)] <- "n/a"
    d2[[cc]] <- v
  }
  d2$.k <- NULL
  write_tsv(d2, f, na = "n/a")
  n_files <- n_files + 1L
  n_filled <- n_filled + sum(d2$adherence_verified != "n/a")
}
cat(sprintf("event files updated: %d  (%d breathing epochs carry physio)\n",
            n_files, n_filled))

# ---- the pooled trial table ------------------------------------------------
bt <- read_csv(file.path(dir_out, "breathing_trials.csv"), show_col_types = FALSE,
               progress = FALSE) %>%
  select(-any_of(c("adherence_verified", "actual_change_physio"))) %>%
  left_join(ep, by = c("subject", "bids_entity", "epoch_index", "trial_type"))
write_csv(bt, file.path(dir_out, "breathing_trials.csv"))

cat(sprintf("\nbreathing_trials.csv: %d rows, %d with respiration measured\n",
            nrow(bt), sum(bt$resp_measured %in% TRUE)))

paced <- bt %>% filter(trial_type == "paced_breathing")
cat("\n-- the analysable set for 'change with and without awareness' --\n")
av <- paced %>% filter(adherence_verified %in% TRUE, !is.na(response))
cat(sprintf("paced epochs with verified adherence and a response: %d (%d participants)\n",
            nrow(av), n_distinct(av$subject)))
cat("\nsdt_class on those epochs (instructed change is valid here):\n")
print(table(av$sdt_class, useNA = "ifany"))
cat("\nby instructed change x salience:\n")
print(as.data.frame(av %>% count(change_direction, salience, sdt_class) %>%
        tidyr::pivot_wider(names_from = sdt_class, values_from = n, values_fill = 0)))
cat(sprintf("\nRun ended: %s\n", Sys.time()))
