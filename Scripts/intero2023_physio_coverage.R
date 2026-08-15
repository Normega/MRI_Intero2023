# intero2023_physio_coverage.R
# Segments every Biopac .acq into acquisition blocks, measures per-run signal
# coverage, and caches a decimated copy so downstream work never re-reads 44 MB
# of 2 kHz data.
#
# Produces:
#   Repo/Data/physio_blocks.csv     every trigger block found, all participants
#   Repo/Data/physio_coverage.csv   one row per canonical BOLD run
#   <cache>/sub-XXX_physio.rds      decimated CO2 + belt + block table
#
# The RDS cache lives OUTSIDE Repo/ by default: it is a regenerable build
# artefact of roughly 4 MB x 56 participants, which does not belong in a
# shareable code repo. Override with --cache=<dir> if you want it elsewhere.
#
# Usage:  Rscript intero2023_physio_coverage.R [--limit=N] [--cache=DIR]
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  m <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else default
}

dir_root    <- "I:/Shared drives/MRI_Intero2023"
dir_physio  <- file.path(dir_root, "Data", "Physio")
dir_out     <- file.path(dir_root, "Repo", "Data")
dir_R       <- file.path(dir_root, "Repo", "R")
dir_cache   <- getarg("cache", file.path(dir_root, "Data", "Physio", "cache"))
f_invent    <- file.path(dir_root, "Prompts", "inventory.tsv")
limit       <- as.integer(getarg("limit", NA))

reticulate::py_require("bioread==2025.5.2")
rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_logs.R", "intero_physio.R"))
  source(file.path(dir_R, f))

dir.create(dir_cache, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Run started: %s\ncache: %s\n", Sys.time(), dir_cache))

# expected volume counts per participant/run, from the dataset of record
inv <- read_tsv(f_invent, show_col_types = FALSE, progress = FALSE) %>%
  transmute(pid = pid_num(subject),
            protocol_slot = as.integer(protocol_run),
            n_volumes = as.integer(n_volumes))

acqs <- list.files(dir_physio, pattern = "^Intero_P\\d+[.]acq$", full.names = TRUE)
pids <- pid_num(basename(acqs))
keep <- !is_pilot(pids) & !is.na(pids)
acqs <- acqs[keep]; pids <- pids[keep]
ord  <- order(pids); acqs <- acqs[ord]; pids <- pids[ord]
if (!is.na(limit)) { acqs <- head(acqs, limit); pids <- head(pids, limit) }

cat(sprintf("participants with .acq (excl. pilot): %d\n", length(acqs)))

all_blocks <- list()

for (i in seq_along(acqs)) {
  pid <- pids[i]; sub <- pid_bids(pid)
  cache_f <- file.path(dir_cache, paste0(sub, "_physio.rds"))
  res <- try(read_acq_session(acqs[i], id = sub), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat(sprintf("  %-9s READ FAILED: %s\n", sub,
                conditionMessage(attr(res, "condition"))))
    all_blocks[[length(all_blocks) + 1]] <-
      tibble(pid = pid, subject = sub, block = NA_integer_, kind = "read_failed")
    next
  }
  exp_i <- inv %>% filter(pid == !!pid) %>% select(protocol_slot, n_volumes)
  res$blocks <- label_bold_blocks(res$blocks, exp_i)

  saveRDS(res, cache_f, compress = "xz")
  all_blocks[[length(all_blocks) + 1]] <-
    res$blocks %>% mutate(pid = pid, subject = sub, .before = 1)

  nb <- sum(res$blocks$kind == "bold"); nf <- sum(res$blocks$kind == "fieldmap")
  nx <- sum(res$blocks$kind == "fragment")
  cat(sprintf("  %-9s %5.1f min  bold=%d fmap=%d frag=%d  matched=%d/%d  cache=%.1fMB\n",
              sub, res$meta$duration_s / 60, nb, nf, nx,
              sum(res$blocks$pulse_match, na.rm = TRUE), nrow(exp_i),
              file.size(cache_f) / 1e6))
}

blocks <- bind_rows(all_blocks)
write_csv(blocks, file.path(dir_out, "physio_blocks.csv"))

# ---- per-run coverage -----------------------------------------------------
coverage <- blocks %>%
  filter(kind == "bold", !is.na(bids_entity)) %>%
  transmute(
    pid, subject, protocol_slot, bids_entity, n_pulses, n_volumes,
    onset_s, dur_s, co2_p99, co2_med, belt_sd, belt_iqr, belt_rng,
    # CO2 above ~2 %CO2 at the 99th percentile means the mask was on; a removed
    # mask collapses to ~0.5. Deliberately a wide margin, not a tuned cut.
    co2_present  = co2_p99 > 2,
    # No amplitude cut for the belt -- gain is not comparable between sessions.
    # Flatline is the one thing that IS comparable, so only that is flagged.
    belt_flatline = belt_iqr < 1e-3
  ) %>%
  arrange(pid, protocol_slot)

write_csv(coverage, file.path(dir_out, "physio_coverage.csv"))

cat("\n-- block classification --\n")
print(as.data.frame(count(blocks, kind)))

cat("\n-- canonical BOLD runs with physio, by task --\n")
print(as.data.frame(
  coverage %>% group_by(bids_entity) %>%
    summarise(n = n(),
              co2_ok = sum(co2_present), co2_lost = sum(!co2_present),
              belt_flat = sum(belt_flatline),
              med_belt_iqr = median(belt_iqr, na.rm = TRUE), .groups = "drop")))

# Two different failures look the same per-run and must not be conflated:
# a capnograph that never recorded at all (every run dead -- the "missing cap"
# hardware fault the RA logged around P068) versus a mask removed for part of
# the session (some runs good, some dead).
sess <- coverage %>%
  group_by(pid, subject) %>%
  summarise(n_runs = n(), n_co2 = sum(co2_present), .groups = "drop") %>%
  mutate(co2_status = case_when(n_co2 == 0        ~ "capnograph_dead_all_session",
                                n_co2 < n_runs    ~ "mask_off_some_runs",
                                TRUE              ~ "ok"))
write_csv(sess, file.path(dir_out, "physio_session_status.csv"))

cat("\n-- session-level CO2 status --\n")
print(as.data.frame(count(sess, co2_status)))

cat("\n-- participants with NO usable CO2 at all --\n")
dead <- sess %>% filter(co2_status == "capnograph_dead_all_session")
cat(if (nrow(dead)) paste("  ", dead$subject, collapse = "\n") else "  none", "\n")

cat("\n-- runs lost mid-session (mask removed) --\n")
partial <- sess %>% filter(co2_status == "mask_off_some_runs") %>% pull(pid)
lost <- coverage %>% filter(!co2_present, pid %in% partial)
if (nrow(lost)) {
  for (i in seq_len(nrow(lost)))
    cat(sprintf("  %-9s %-24s co2_p99=%.2f belt_iqr=%.4f\n",
                lost$subject[i], lost$bids_entity[i], lost$co2_p99[i], lost$belt_iqr[i]))
} else cat("  none\n")

cat(sprintf("\nblocks: %d rows -> physio_blocks.csv\ncoverage: %d runs -> physio_coverage.csv\n",
            nrow(blocks), nrow(coverage)))
cat(sprintf("Run ended: %s\n", Sys.time()))
