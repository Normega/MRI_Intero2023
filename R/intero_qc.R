# intero_qc.R
# A single append-only QC ledger for the whole project.
#
# Long format, one row per assessment, so any stage -- log parsing, physio,
# events, imaging -- can add checks without changing a schema. A wide
# one-row-per-participant table would have to be renegotiated every time a new
# check appears, and would hide the fact that most checks are per RUN, not per
# participant.
#
#   subject      sub-058
#   scope        session | run
#   entity       task-breathing_run-1, or NA for session-level
#   stage        logs | physio | events | imaging
#   check        acq_log_time_match
#   value        0.4                     numeric where meaningful
#   grade        exact | good | degraded | ambiguous | fail | info
#   note         free text
#   source       what produced the number
#   date_checked
#
# `grade` is deliberately coarse and ordered. Nothing here decides inclusion --
# it records what was measured so exclusion stays a downstream, reviewable call.
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(readr) })

QC_GRADES <- c("exact", "good", "info", "degraded", "ambiguous", "fail")

#' One QC row. Vectorised over its arguments.
qc_entry <- function(subject, stage, check, grade,
                     scope = "run", entity = NA_character_,
                     value = NA_real_, note = NA_character_,
                     source = NA_character_) {
  stopifnot(all(grade %in% QC_GRADES))
  tibble(subject = subject, scope = scope, entity = entity, stage = stage,
         check = check, value = as.numeric(value), grade = grade,
         note = note, source = source,
         date_checked = as.character(Sys.Date()))
}

#' Merge new rows into the ledger, replacing any prior entry for the same
#' (subject, entity, stage, check) so re-running a stage updates rather than
#' duplicates.
qc_merge <- function(ledger, new) {
  if (is.null(ledger) || !nrow(ledger)) return(new)
  key <- function(d) paste(d$subject, d$entity, d$stage, d$check, sep = "\r")
  ledger[!key(ledger) %in% key(new), , drop = FALSE] %>% bind_rows(new)
}

#' Read the ledger with column types PINNED.
#'
#' Type inference is not safe here: a column that happens to be all-NA in the
#' stored file (entity on a session-scoped run, or note) is guessed as logical,
#' and binding logical to character fails in vctrs. Pinning also keeps `value`
#' numeric when a stage happens to write only NAs.
qc_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE, progress = FALSE,
           col_types = cols(
             subject = col_character(), scope = col_character(),
             entity = col_character(),  stage = col_character(),
             check = col_character(),   value = col_double(),
             grade = col_character(),   note = col_character(),
             source = col_character(),  date_checked = col_character()))
}

qc_write <- function(ledger, path) {
  ledger %>%
    mutate(grade = factor(.data$grade, levels = QC_GRADES)) %>%
    arrange(.data$subject, .data$stage, .data$entity, .data$check) %>%
    write_csv(path)
  invisible(ledger)
}

#' Grade an acq-block <-> log-run time match.
#'
#' Two numbers matter and they mean different things:
#'   dt_s     how close the best candidate is. Small dt is necessary.
#'   margin_s how much closer the best is than the runner-up. This is what makes
#'            the match *identifiable*. A 0.5 s match whose runner-up is 600 s
#'            away is certain; the same 0.5 s with a runner-up at 3 s is a coin
#'            flip, and pulse counts cannot break the tie when two runs of the
#'            same task have identical lengths (loc1 and loc2 are both 190).
#' pulse_delta corroborates but must not veto: v1 acquisitions were clipped by
#'   up to ~20 volumes, so a large delta with a tight time match is a clipped
#'   run, not a mismatch.
grade_time_match <- function(dt_s, margin_s, pulse_delta) {
  dplyr::case_when(
    is.na(dt_s)                              ~ "fail",
    dt_s > QC_TIME_TOL_S                     ~ "fail",
    !is.na(margin_s) & margin_s < QC_MARGIN_MIN_S ~ "ambiguous",
    dt_s <= 2 & !is.na(pulse_delta) & pulse_delta == 0 ~ "exact",
    !is.na(pulse_delta) & pulse_delta > QC_PULSE_TOL   ~ "degraded",
    TRUE                                     ~ "good"
  )
}

QC_TIME_TOL_S    <- 15   # beyond this the block is not that run
QC_MARGIN_MIN_S  <- 30   # runner-up must be at least this much worse
QC_PULSE_TOL     <- 5    # pulses; larger implies a clipped acquisition
