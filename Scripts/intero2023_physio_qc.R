# intero2023_physio_qc.R
# Multi-part .acq ingestion, wallclock run matching, and the three-way QC gate.
# Supersedes intero2023_physio_coverage.R.
#
# THE THREE CHECKS, and why all three are needed:
#   volume_count_vs_mode  n_volumes vs the cohort mode FOR ITS PROTOCOL VERSION.
#     v1 and v2 have different run lengths (breathing 265 vs 242), so comparing
#     across the split invents deviations and hides real ones. Catches conversion
#     truncation (sub-065: 167 of 196) and scanner interruptions (sub-094).
#   pulse_vs_volume       acq trigger pulses vs NIfTI volumes. Confirms the
#     physio saw the acquisition that survived into BIDS.
#   acq_log_time_match    acq block start vs the log's first-trigger wallclock.
#     Confirms the TASK that ran is the one on those images.
#
# A run needs all three. sub-056 passes pulse_vs_volume and fails
# acq_log_time_match -- the signature of images and events coming from different
# acquisitions, which is the failure most likely to reach results unnoticed.
#
# Produces: Repo/Data/{physio_blocks,physio_coverage,qc_ledger}.csv
#
# Script version: 2.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_phys <- file.path(dir_root, "Data", "Physio")
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_R    <- file.path(dir_root, "Repo", "R")
dir_cache<- file.path(dir_phys, "cache")

reticulate::py_require("bioread==2025.5.2")
rk <- "I:/Shared drives/TechnicalReference/respkit/R"
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_logs.R", "intero_physio.R", "intero_qc.R"))
  source(file.path(dir_R, f))

dir.create(dir_cache, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Run started: %s\n", Sys.time()))

# ---- reference tables -----------------------------------------------------
inv <- read_tsv(file.path(dir_root, "Prompts", "inventory.tsv"),
                show_col_types = FALSE, progress = FALSE) %>%
  mutate(pid = pid_num(subject),
         protocol_slot = as.integer(protocol_run),
         n_volumes = as.integer(n_volumes),
         version = protocol_version(pid))

# cohort mode per (slot, version) -- the v1/v2 split is essential here
modes <- inv %>% group_by(protocol_slot, version) %>%
  summarise(mode_vol = as.integer(names(sort(table(n_volumes), decreasing = TRUE))[1]),
            n_cohort = n(), .groups = "drop")

logs <- read_csv(file.path(dir_out, "log_runs.csv"),
                 show_col_types = FALSE, progress = FALSE)

idx <- acq_file_index(dir_phys) %>% filter(!is_pilot(pid))
cat(sprintf("acq files: %d, participants: %d, multi-part: %s\n",
            nrow(idx), n_distinct(idx$pid),
            paste(sort(unique(idx$pid[duplicated(idx$pid)])), collapse = ",")))

all_blocks <- list(); qc <- list()

for (P in sort(unique(idx$pid))) {
  sub <- pid_bids(P); parts <- idx %>% filter(pid == P)
  exp_i <- inv %>% filter(pid == P) %>% select(protocol_slot, n_volumes, version)
  lr <- logs %>% filter(pid == P, has_trigger, !is.na(protocol_slot)) %>%
        transmute(protocol_slot, trigger_time = as.POSIXct(trigger_time, tz = "UTC"))

  bl <- map_dfr(seq_len(nrow(parts)), function(i) {
    r <- try(read_acq_session(parts$path[i], id = sub), silent = TRUE)
    if (inherits(r, "try-error")) return(tibble())
    t0 <- acq_start_time(parts$path[i])
    if (!is.null(r$co2))
      saveRDS(r, file.path(dir_cache, sprintf("%s_part-%s_physio.rds", sub, parts$part[i])),
              compress = "xz")
    if (!nrow(r$blocks)) return(tibble())
    r$blocks %>% mutate(part = parts$part[i], abs_time = t0 + onset_s, .before = 1)
  })

  if (!nrow(bl)) {
    # Recorded explicitly. A participant with no trigger blocks must still
    # appear, or the table stops reconciling against the participant list.
    all_blocks[[length(all_blocks)+1]] <- tibble(subject = sub, pid = P,
                                                 kind = "no_trigger_blocks")
    qc[[length(qc)+1]] <- qc_entry(sub, "physio", "acq_trigger_present", "fail",
      scope = "session", note = "no trigger blocks detected in any .acq part",
      source = paste(parts$file, collapse = ";"))
    cat(sprintf("  %-9s NO TRIGGER BLOCKS (%d part(s))\n", sub, nrow(parts)))
    next
  }

  bl <- bl %>% arrange(abs_time)
  bold <- bl %>% filter(kind == "bold")
  m <- match_blocks_by_time(bold, lr, exp_i %>% select(protocol_slot, n_volumes))
  m$grade_time <- grade_time_match(m$dt_s, m$margin_s, m$pulse_delta)

  # join cohort mode for the volume-count check
  m <- m %>%
    left_join(exp_i %>% select(protocol_slot, version), by = "protocol_slot") %>%
    left_join(modes, by = c("protocol_slot", "version")) %>%
    mutate(vol_delta = n_volumes - mode_vol,
           vol_pct   = 100 * vol_delta / mode_vol)

  all_blocks[[length(all_blocks)+1]] <-
    bind_rows(m, bl %>% filter(kind != "bold")) %>%
    mutate(subject = sub, pid = P, .before = 1)

  for (i in seq_len(nrow(m))) {
    ent <- m$bids_entity[i]; if (is.na(ent)) ent <- sprintf("block%d", i)
    qc[[length(qc)+1]] <- qc_entry(sub, "physio", "acq_log_time_match",
      m$grade_time[i], entity = ent, value = m$dt_s[i],
      note = sprintf("margin=%.0fs part=%s", m$margin_s[i], m$part[i]),
      source = "acq header vs log trigger_time")
    if (!is.na(m$pulse_delta[i]))
      qc[[length(qc)+1]] <- qc_entry(sub, "physio", "pulse_vs_volume",
        if (m$pulse_delta[i] <= 1) "exact" else if (m$pulse_delta[i] <= QC_PULSE_TOL) "good" else "degraded",
        entity = ent, value = m$pulse_delta[i],
        note = sprintf("%d pulses vs %d volumes", m$n_pulses[i], m$n_volumes[i]),
        source = "acq trigger count vs NIfTI header")
    if (!is.na(m$vol_delta[i]))
      qc[[length(qc)+1]] <- qc_entry(sub, "imaging", "volume_count_vs_mode",
        if (m$vol_delta[i] == 0) "exact" else if (abs(m$vol_pct[i]) < 2) "good"
        else if (abs(m$vol_pct[i]) < 10) "degraded" else "fail",
        entity = ent, value = m$vol_delta[i],
        note = sprintf("%d vs %s mode %d (%.1f%%)", m$n_volumes[i],
                       m$version[i], m$mode_vol[i], m$vol_pct[i]),
        source = "inventory.tsv vs cohort mode")
    if (!is.na(m$co2_p99[i]))
      qc[[length(qc)+1]] <- qc_entry(sub, "physio", "co2_present",
        if (m$co2_p99[i] > 2) "good" else "fail", entity = ent,
        value = m$co2_p99[i], note = "p99 %CO2; <2 means mask off or capnograph dead",
        source = "CO2100C channel")
    if (!is.na(m$belt_iqr[i]))
      qc[[length(qc)+1]] <- qc_entry(sub, "physio", "belt_flatline",
        if (m$belt_iqr[i] < 1e-3) "fail" else "info", entity = ent,
        value = m$belt_iqr[i], note = "IQR volts; amplitude not comparable across sessions",
        source = "DA100C channel")
  }
  cat(sprintf("  %-9s parts=%-7s bold=%d matched=%d exact=%d\n", sub,
              paste(parts$part, collapse = "/"), nrow(bold),
              sum(!is.na(m$bids_entity)), sum(m$grade_time == "exact")))
}

blocks <- bind_rows(all_blocks)
ledger <- bind_rows(qc)
write_csv(blocks, file.path(dir_out, "physio_blocks.csv"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), ledger),
         file.path(dir_out, "qc_ledger.csv"))

cov <- blocks %>% filter(kind == "bold", !is.na(bids_entity), !dup_slot) %>%
  select(subject, pid, protocol_slot, bids_entity, n_pulses, n_volumes,
         dt_s, margin_s, pulse_delta, vol_delta, grade_time,
         co2_p99, belt_iqr, part)
write_csv(cov, file.path(dir_out, "physio_coverage.csv"))

cat("\n-- QC ledger by check x grade --\n")
print(as.data.frame(ledger %>% count(check, grade) %>% tidyr::pivot_wider(
  names_from = grade, values_from = n, values_fill = 0)))
cat(sprintf("\nblocks=%d  coverage runs=%d  qc rows=%d  subjects=%d\n",
            nrow(blocks), nrow(cov), nrow(ledger), n_distinct(blocks$subject)))

cat("\n-- runs failing acq_log_time_match --\n")
bad <- ledger %>% filter(check == "acq_log_time_match", grade %in% c("fail","ambiguous"))
if (nrow(bad)) for (i in seq_len(nrow(bad)))
  cat(sprintf("  %-9s %-24s dt=%.0fs %s\n", bad$subject[i], bad$entity[i],
              bad$value[i], bad$note[i])) else cat("  none\n")
cat(sprintf("\nRun ended: %s\n", Sys.time()))
