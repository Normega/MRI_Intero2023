# intero2023_h3_sets.R
# Freezes the H3 participant sets of prereg v0.7 from the first-level event
# sheets (trial counts only; no imaging data). Added 2026-11-25, pre-freeze.
#
#   in_H3  (primary, N = 45): >= 1 paced_aware AND >= 1 paced_unaware pooled
#          across gate-passing runs. v0.7 fits H3 in one joint model of all the
#          participant's gate-passing runs (h3_runs), so every trial enters.
#          (In these data the set is identical to v0.6's "both cells in one
#          run" rule; the script asserts this.)
#   n_aware_sharedrun / n_unaware_sharedrun: trials in runs that contain both
#          cells -- what the per-run robustness check (v0.6 form) uses.
#   in_S7  (sensitivity S7 and H7b, N = 23): the v0.5 rule, >= 3 aware and
#          >= 3 unaware trials pooled across gate-passing runs.
#
# Output: Data/h3_sets.csv, one row per analysable participant (47).
# Script version: 1.1 (2026-11-25: joint two-run H3 model)
# Author: RADlab

suppressPackageStartupMessages({ library(readr); library(dplyr) })
dir_repo <- Sys.getenv("INTERO_REPO", "/media/normega/FB_Cluster_Data/intero2024/MRI_Intero2023")
fl <- file.path(dir_repo, "FirstLevel")
m <- bind_rows(lapply(c("firstlevel_manifest.csv", "firstlevel_manifest_candidates.csv"), function(f)
  read_csv(file.path(fl, f), col_types = cols(.default = "c")))) %>%
  dplyr::filter(motion_gate_fail == "FALSE", subject != "sub-056") %>% distinct(subject, bids_entity)
runs <- bind_rows(lapply(seq_len(nrow(m)), function(i) {
  e <- read_tsv(file.path(fl, sprintf("%s_%s_desc-model_events.tsv", m$subject[i], m$bids_entity[i])),
                col_types = cols(.default = "c"))
  tibble(subject = m$subject[i], bids_entity = m$bids_entity[i],
         a = sum(e$condition == "paced_aware"), u = sum(e$condition == "paced_unaware"))
})) %>% mutate(contrast_run = a >= 1 & u >= 1)
sets <- runs %>% group_by(subject) %>% summarise(
  n_runs = n(), n_aware_pooled = sum(a), n_unaware_pooled = sum(u),
  h3_runs = paste(sub("task-breathing_", "", bids_entity), collapse = ";"),
  n_aware_sharedrun = sum(a[contrast_run]), n_unaware_sharedrun = sum(u[contrast_run]),
  any_shared_run = any(contrast_run), .groups = "drop") %>%
  mutate(in_H3 = n_aware_pooled >= 1 & n_unaware_pooled >= 1,
         in_S7 = n_aware_pooled >= 3 & n_unaware_pooled >= 3)
stopifnot(nrow(sets) == 47, sum(sets$in_H3) == 45, sum(sets$in_S7) == 23, all(!sets$in_S7 | sets$in_H3),
          identical(sets$in_H3, sets$any_shared_run))
sets <- sets %>% select(-any_shared_run)
sets$h3_runs[!sets$in_H3] <- ""
write_csv(sets, file.path(dir_repo, "Data", "h3_sets.csv"))
cat(sprintf("H3 %d, S7 %d; excluded from H3: %s\n", sum(sets$in_H3), sum(sets$in_S7),
            paste(sets$subject[!sets$in_H3], collapse = ", ")))
