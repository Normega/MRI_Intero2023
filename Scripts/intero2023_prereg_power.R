# intero2023_prereg_power.R
# Establishes the analysable sample and per-cell trial counts that the
# preregistration commits to. Run BEFORE any imaging analysis.
#
# A first-level contrast is only estimable for a participant who has at least
# one trial in each of its cells, so "N = 54" is meaningless here -- the real N
# differs per contrast and per physiological covariate set. This script reports
# each one explicitly so the preregistration states numbers it can be held to.
#
# Produces: Repo/Data/prereg_sample.csv
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_reg  <- file.path(dir_out, "regressors")
dir_ready<- file.path(dir_root, "Data", "Ready")

EXCLUDE_LABEL <- c("sub-082", "sub-104")   # shifted run labels
MIN_TRIALS <- 3L                            # per cell, per participant

b <- read_csv(file.path(dir_out, "breathing_trials.csv"), show_col_types = FALSE,
              progress = FALSE) %>%
  mutate(confidence = suppressWarnings(as.numeric(confidence)))

TASK_MAP <- c(breath1 = "task-breathing_run-1", breath2 = "task-breathing_run-2")
img <- tibble(path = list.files(dir_ready, pattern = "desc-confounds_timeseries[.]tsv$",
                                recursive = TRUE, full.names = TRUE)) %>%
  mutate(bn = basename(path),
         subject = str_match(bn, "^(sub-[0-9]+)_")[, 2],
         rtask = str_match(bn, "_task-([a-z0-9]+)_")[, 2],
         bids_entity = unname(TASK_MAP[rtask])) %>%
  dplyr::filter(!is.na(bids_entity)) %>%
  select(subject, bids_entity) %>% distinct() %>% mutate(has_img = TRUE)

reg <- tibble(f = list.files(dir_reg, pattern = "[.]tsv$")) %>%
  mutate(subject = str_match(f, "^(sub-[0-9]+)_")[, 2],
         bids_entity = str_match(f, "^sub-[0-9]+_(task-[a-z]+(?:_run-[0-9]+)?)_")[, 2],
         kind = ifelse(grepl("_petco2[.]tsv$", f), "petco2", "resp")) %>%
  dplyr::filter(!is.na(bids_entity), grepl("breathing", bids_entity)) %>%
  distinct(subject, bids_entity, kind) %>%
  mutate(v = TRUE) %>% tidyr::pivot_wider(names_from = kind, values_from = v, values_fill = FALSE)

runs <- b %>% dplyr::filter(trial_type == "paced_breathing") %>%
  distinct(subject, bids_entity) %>%
  left_join(img, by = c("subject", "bids_entity")) %>%
  left_join(reg, by = c("subject", "bids_entity")) %>%
  mutate(across(c(has_img, resp, petco2), ~ tidyr::replace_na(.x, FALSE)),
         label_ok = !subject %in% EXCLUDE_LABEL)

cat(sprintf("Run started: %s\n\n", Sys.time()))
cat("=== BREATHING RUNS: what exists ===\n")
cat(sprintf("  runs with events            : %d\n", nrow(runs)))
cat(sprintf("  + usable imaging            : %d\n", sum(runs$has_img & runs$label_ok)))
cat(sprintf("  + RVT regressor             : %d\n", sum(runs$has_img & runs$label_ok & runs$resp)))
cat(sprintf("  + PetCO2 regressor          : %d\n", sum(runs$has_img & runs$label_ok & runs$petco2)))
cat(sprintf("  + BOTH physio regressors    : %d\n",
            sum(runs$has_img & runs$label_ok & runs$resp & runs$petco2)))

usable <- runs %>% dplyr::filter(has_img, label_ok)

# ---- per-participant trial counts in each modelled condition ---------------
ch <- b %>% dplyr::filter(trial_type == "paced_breathing", !is.na(response)) %>%
  semi_join(usable, by = c("subject", "bids_entity")) %>%
  mutate(cell = case_when(
    change_direction == "No Change" ~ "nochange_trial",
    change_direction == response    ~ "aware_correct",
    response == "No Change"         ~ "unaware_saidNC",
    TRUE                            ~ "wrong_direction"))
fr <- b %>% dplyr::filter(trial_type == "free_breathing") %>%
  semi_join(usable, by = c("subject", "bids_entity")) %>% mutate(cell = "free")

cnt <- bind_rows(ch %>% select(subject, cell), fr %>% select(subject, cell)) %>%
  count(subject, cell) %>% tidyr::pivot_wider(names_from = cell, values_from = n, values_fill = 0)

cat("\n=== TRIALS PER PARTICIPANT, by modelled condition ===\n")
print(as.data.frame(cnt %>% summarise(across(-subject, list(
  median = ~ median(.x), min = ~ min(.x), max = ~ max(.x))))) %>% t())

cat("\n=== ESTIMABILITY OF EACH PLANNED CONTRAST ===\n")
est <- tibble(
  contrast = c("C1 paced vs free",
               "C2 aware vs unaware (categorical)",
               "C3 confidence parametric (all change trials)",
               "C4 wrong direction vs aware",
               "C5 salience high vs low (paced)"),
  n_participants = c(
    sum(cnt$aware_correct + cnt$unaware_saidNC + cnt$wrong_direction >= MIN_TRIALS &
          cnt$free >= MIN_TRIALS),
    sum(cnt$aware_correct >= MIN_TRIALS & cnt$unaware_saidNC >= MIN_TRIALS),
    sum(cnt$aware_correct + cnt$unaware_saidNC + cnt$wrong_direction >= 2 * MIN_TRIALS),
    sum(cnt$wrong_direction >= MIN_TRIALS & cnt$aware_correct >= MIN_TRIALS),
    NA_integer_))
sal <- ch %>% dplyr::filter(change_direction != "No Change") %>%
  count(subject, salience) %>% tidyr::pivot_wider(names_from = salience, values_from = n, values_fill = 0)
est$n_participants[5] <- sum(sal$high >= MIN_TRIALS & sal$low >= MIN_TRIALS)
print(as.data.frame(est))
cat(sprintf("\n(minimum %d trials per cell; total participants with usable imaging: %d)\n",
            MIN_TRIALS, n_distinct(usable$subject)))

cat("\n=== PHYSIOLOGICAL COVARIATE AVAILABILITY, participant level ===\n")
ps <- usable %>% group_by(subject) %>%
  summarise(n_runs = n(), rvt_runs = sum(resp), co2_runs = sum(petco2), .groups = "drop")
cat(sprintf("  participants with >=1 usable breathing run : %d\n", nrow(ps)))
cat(sprintf("  with RVT on all their runs                 : %d\n", sum(ps$rvt_runs == ps$n_runs)))
cat(sprintf("  with PetCO2 on all their runs              : %d\n", sum(ps$co2_runs == ps$n_runs)))
cat(sprintf("  with BOTH on all their runs                : %d\n",
            sum(ps$rvt_runs == ps$n_runs & ps$co2_runs == ps$n_runs)))

adh <- b %>% dplyr::filter(trial_type == "paced_breathing") %>%
  semi_join(usable, by = c("subject", "bids_entity")) %>%
  group_by(subject) %>% summarise(adh = any(adherence_verified %in% TRUE), .groups = "drop")
cat(sprintf("  with verified pacer adherence on >=1 run   : %d\n", sum(adh$adh)))

out <- cnt %>% left_join(ps, by = "subject") %>% left_join(adh, by = "subject")
write_csv(out, file.path(dir_out, "prereg_sample.csv"))
cat(sprintf("\nwrote %s\n", file.path(dir_out, "prereg_sample.csv")))
cat(sprintf("Run ended: %s\n", Sys.time()))
