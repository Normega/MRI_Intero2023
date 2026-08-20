# intero_ids.R
# Canonical participant-ID handling for Intero2023 / breathMRI2026.
#
# Every source in this project spells the participant differently:
#   PsychoPy logs      P50, P130          (Pnum column: 50, 130)
#   run_notes/screener P052 ... but also a stray unpadded P90
#   BIDS / fMRIPrep    sub-052
#   Biopac physio      Intero_P52.acq
#   Raw MRI dirs       P052_20991231   (synthetic example; not a real session date)
#
# Everything joins through pid_num() (an integer) and is displayed as
# pid_bids() ("sub-052"). Nothing downstream should parse an ID by hand.
#
# Author: RADlab

# ---- core normalisation ---------------------------------------------------

#' Extract the canonical integer participant number from any project ID spelling.
#'
#' Accepts "P50", "P050", "sub-050", "Intero_P52.acq", "P052_20991231", 50, "50".
#' Returns NA_integer_ for anything with no recognisable participant number.
pid_num <- function(x) {
  x <- as.character(x)
  out <- rep(NA_integer_, length(x))
  if (!length(x)) return(out)

  # Bare numbers ("50", 130) are already the number.
  bare <- !is.na(x) & grepl("^\\s*\\d+\\s*$", x)
  out[bare] <- suppressWarnings(as.integer(trimws(x[bare])))

  # Everything else must carry an explicit P### or sub-### marker. We do NOT
  # fall back to "first number anywhere" -- filenames here are full of dates
  # and timestamps that would match (P052_20991231 would yield 20991231).
  #
  # No \b before P: "_" is a word character, so a word boundary would refuse
  # to match "Intero_P52.acq", which is exactly how the physio files are named.
  idx <- which(!bare & !is.na(x))
  if (length(idx)) {
    # regexec (not regexpr) so the digits come back as a capture group;
    # regexpr + regmatches silently drops non-matching elements and
    # misaligns the result against the input vector.
    m <- regexec("(?:sub-|P)0*([0-9]+)", x[idx], ignore.case = TRUE)
    got <- regmatches(x[idx], m)
    out[idx] <- vapply(got, function(z) {
      if (length(z) >= 2L) suppressWarnings(as.integer(z[2])) else NA_integer_
    }, integer(1))
  }
  out
}

#' Canonical BIDS-style label, e.g. 52 -> "sub-052".
pid_bids <- function(x) {
  n <- pid_num(x)
  ifelse(is.na(n), NA_character_, sprintf("sub-%03d", n))
}

#' Canonical lab-style label, e.g. 52 -> "P052".
pid_lab <- function(x) {
  n <- pid_num(x)
  ifelse(is.na(n), NA_character_, sprintf("P%03d", n))
}

# ---- study roles ----------------------------------------------------------

# Confirmed by PI 2026-08-05. P50 and P999 are pilot/QA; P00 and P01 are
# pre-study test records. These are flagged, never silently dropped, so that
# participant-flow counts stay reportable.
INTERO_PILOT_IDS <- c(0L, 1L, 50L, 999L)

#' TRUE for pilot / QA records that are not real participants.
is_pilot <- function(x) {
  n <- pid_num(x)
  !is.na(n) & n %in% INTERO_PILOT_IDS
}

# The task was revised after the early cohort: from P60 onward the breathing
# runs log explicit ChangeTrial rows and the paced rate was fixed at 12 brpm.
# P51-P59 ran the earlier version and are retained for now but must stay
# separable -- see intero2023-pilot-ids.
INTERO_V2_FIRST_ID <- 60L

#' Protocol version: "v1" for the early cohort, "v2" from P60 on.
protocol_version <- function(x) {
  n <- pid_num(x)
  ifelse(is.na(n), NA_character_,
         ifelse(n >= INTERO_V2_FIRST_ID, "v2", "v1"))
}

# ---- run naming -----------------------------------------------------------

# The dataset of record is the Linux box copy, which uses paradigm labels with
# BIDS run- entities. The Google Drive Data/BIDS copy uses loc1/breath1-style
# names, is missing all fieldmaps, and has shifted labels for sub-082/sub-104.
# Emit the Linux spelling: anything named task-breath1_events.tsv pairs with
# nothing. See intero2023-bids-mislabelled-runs.
#
# protocol_slot is the Run<N> in the scanner's ProtocolName and is the only
# trustworthy run identifier.
INTERO_RUNS <- data.frame(
  log_run       = c("Localizer1", "Localizer2", "ChangeDetection1",
                    "ChangeDetection2", "Authenticity"),
  protocol_slot = c(1L, 2L, 3L, 4L, 5L),
  bids_task     = c("localizer", "localizer", "breathing",
                    "breathing", "authenticity"),
  bids_run      = c("1", "2", "1", "2", NA_character_),
  stringsAsFactors = FALSE
)

#' BIDS entity string for a log Run label, e.g. "Localizer2" -> "task-localizer_run-2".
bids_entities <- function(log_run) {
  i <- match(log_run, INTERO_RUNS$log_run)
  task <- INTERO_RUNS$bids_task[i]
  run  <- INTERO_RUNS$bids_run[i]
  ifelse(is.na(task), NA_character_,
         ifelse(is.na(run),
                paste0("task-", task),
                paste0("task-", task, "_run-", run)))
}
