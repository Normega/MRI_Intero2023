# intero_logs.R
# Parsing and session assembly for the Intero2023 PsychoPy task logs.
#
# ---------------------------------------------------------------------------
# WHAT THESE LOGS ACTUALLY ARE  -- read this before changing anything.
#
# A log row is a *state snapshot at an event*, not a trial record. Two
# consequences, both verified against Data/Logs on 2026-08-05:
#
#   1. RESPONSES LAG BY ONE ROW. The value in Response/RT/Accuracy/Confidence/
#      Valence/Arousal/Authenticity is whatever had been recorded at the moment
#      the row was written -- i.e. it belongs to the PREVIOUS logged event.
#      Verified on P130's Authenticity run: 16 stimuli from Values/P130.txt,
#      the first stimulus logs a blank rating, and the last stimulus's rating
#      lands on the trailing Completion row.
#
#   2. VALUES GO STALE. When no new response occurs, the previous value simply
#      persists into later rows (Fixation, Free, Completion, the next trial).
#      A non-empty Response does NOT mean a response happened. In P130's
#      breathing run RT=0.5939914999980829 repeats across six later rows.
#
# Naively reading these columns fabricates data. Use shift_responses().
#
# TIMING. At the PulseReceived row, RunStartTime is rewritten to the first
# scanner-trigger wallclock and RunTime resets to zero. So RunTime is already
# time-since-first-TR and is the correct onset clock. Only the FIRST trigger of
# each run is logged, so the logs alone cannot detect dropped volumes -- that
# needs the Biopac trigger channel and the NIfTI volume counts.
#
# Author: RADlab

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

# ---- small helpers --------------------------------------------------------

#' Parse PsychoPy's "H:MM:SS.ffffff" elapsed-time strings to seconds.
#'
#' Returns NA for blanks and for anything that is not H:M:S. Deliberately
#' strict: a silently-coerced 0 here would look like a valid onset.
parse_hms <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_real_, length(x))
  ok <- grepl("^\\d+:\\d{1,2}:\\d{1,2}(\\.\\d+)?$", x)
  if (any(ok)) {
    p <- str_split_fixed(x[ok], ":", 3)
    out[ok] <- as.numeric(p[, 1]) * 3600 + as.numeric(p[, 2]) * 60 + as.numeric(p[, 3])
  }
  out
}

#' Parse the "2024-06-03 11:02:37.247504" wallclock stamps.
parse_wallclock <- function(x) {
  as.POSIXct(trimws(as.character(x)), format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- filename metadata ----------------------------------------------------

LOG_FILE_RE <- paste0(
  "^Log\\.Intero2023\\.",
  "(P\\d+)\\.",              # participant
  "Session(\\w+)_",          # Test / Train
  "(\\d{4}-\\d{2}-\\d{2})_", # date
  "(\\d{1,2})h(\\d{2})\\.(\\d{2})\\.(\\d+)"
)

#' Pull participant / session / launch-time out of a log filename.
#'
#' Tolerates the " (1)" suffix Drive adds to duplicated downloads.
parse_log_filename <- function(fn) {
  base <- basename(fn)
  m <- str_match(base, LOG_FILE_RE)
  tibble(
    file          = base,
    pid           = pid_num(m[, 2]),
    session_type  = m[, 3],
    session_date  = m[, 4],
    launch_time   = ifelse(is.na(m[, 5]), NA_character_,
                           sprintf("%02d:%s:%s", as.integer(m[, 5]), m[, 6], m[, 7])),
    is_dup_suffix = grepl("\\(\\d+\\)\\.csv$", base),
    parsed_ok     = !is.na(m[, 2])
  )
}

# ---- reading --------------------------------------------------------------

# Every column is read as text. Type guessing on these files is actively
# harmful: Response holds both integers and words ("No Change", "Slower"),
# and RespiratoryRate appears as both "12" and "12.0".
LOG_COLS <- c("ExpName", "Pnum", "Session", "RespiratoryRate", "ReferenceTime",
              "FirstRun", "StartTime", "ExptStartTime", "RunStartTime",
              "BlockStartTime", "TotalTime", "RunTime", "BlockTime", "Run",
              "Block", "Description", "PulseTime", "Pause", "TrialNum",
              "StimOnset", "Direction", "Stimulus", "RT", "Valence", "Arousal",
              "TrialType", "ChangeAmount", "Salience", "Response", "Accuracy",
              "Confidence", "Authenticity", "AvgACC", "Threshold")

#' Read one raw log CSV. Returns a zero-row tibble for empty/headerless files.
read_log_file <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(tibble())
  out <- try(suppressWarnings(readr::read_csv(
    path,
    col_types   = cols(.default = col_character()),
    name_repair = "unique_quiet",
    progress    = FALSE
  )), silent = TRUE)
  if (inherits(out, "try-error") || nrow(out) == 0) return(tibble())

  # Trailing unnamed column from the writer's stray comma.
  out <- out[, names(out) %in% LOG_COLS, drop = FALSE]
  for (cn in setdiff(LOG_COLS, names(out))) out[[cn]] <- NA_character_
  out <- out[, LOG_COLS, drop = FALSE]

  out %>%
    mutate(
      .row       = row_number(),
      run_time_s = parse_hms(.data$RunTime),
      total_s    = parse_hms(.data$TotalTime),
      run_start  = parse_wallclock(.data$RunStartTime),
      trial_num  = suppressWarnings(as.integer(.data$TrialNum))
    )
}

# ---- the off-by-one / staleness fix ---------------------------------------

#' Recover per-event responses from the lagged, carry-forward response columns.
#'
#' The response to the event on row i is written on row i+1. It is only a
#' genuine response if row i+1 actually recorded something new -- otherwise the
#' previous value is just persisting.
#'
#' Freshness test, in order of preference:
#'   * if RT is populated for this run, a new response occurred when RT changed;
#'   * otherwise (rating runs log no RT) fall back to "the target column is
#'     non-empty", and separately flag consecutive-identical values, which are
#'     ambiguous between a real repeat and a missed response.
#'
#' @param df   rows of a single run, in log order
#' @param cols response-bearing columns to recover
#' @return df with `<col>_resp` columns added, plus `resp_fresh` and
#'   `resp_repeat` diagnostics.
shift_responses <- function(df, cols = c("Response", "RT", "Accuracy",
                                         "Confidence", "Valence", "Arousal",
                                         "Authenticity")) {
  if (nrow(df) == 0) return(df)
  cols <- intersect(cols, names(df))

  nxt <- function(x) c(x[-1], NA)              # value from the following row
  blank <- function(x) is.na(x) | trimws(x) == ""

  rt_raw    <- if ("RT" %in% names(df)) suppressWarnings(as.numeric(df$RT)) else rep(NA_real_, nrow(df))
  rt_usable <- any(is.finite(rt_raw))

  if (rt_usable) {
    # A new response was recorded on row i when RT differs from row i-1.
    prev_rt <- c(NA_real_, rt_raw[-length(rt_raw)])
    new_here <- is.finite(rt_raw) & (!is.finite(prev_rt) | rt_raw != prev_rt)
  } else {
    new_here <- rep(NA, nrow(df))
  }

  for (cn in cols) {
    v_next <- nxt(df[[cn]])
    if (rt_usable) {
      fresh_next <- nxt(new_here)
      fresh_next[is.na(fresh_next)] <- FALSE
    } else {
      fresh_next <- !blank(v_next)
    }
    df[[paste0(cn, "_resp")]] <- ifelse(fresh_next, v_next, NA_character_)
  }

  # Diagnostics, kept alongside rather than applied, so the caller can judge.
  if (rt_usable) {
    df$resp_fresh <- nxt(new_here)
    df$resp_fresh[is.na(df$resp_fresh)] <- FALSE
    df$resp_repeat <- FALSE
  } else {
    tgt <- if ("Authenticity" %in% cols) "Authenticity" else cols[1]
    v_next <- nxt(df[[tgt]])
    df$resp_fresh  <- !blank(v_next)
    df$resp_repeat <- !blank(v_next) & !is.na(df[[tgt]]) &
      trimws(v_next) == trimws(df[[tgt]] %||% "")
    df$resp_repeat[is.na(df$resp_repeat)] <- FALSE
  }
  df
}

# ---- per-run extraction ---------------------------------------------------

TRIAL_DESCS <- c("Trial", "ChangeTrial")
DONE_DESCS  <- c("Completion", "BlockComplete")

#' Summarise every scanner run present in one log file.
#'
#' One row per (file, Run). StartUp and the pilot-only ChangeThreshold
#' staircase are excluded from the scanner-run set but reported in `other_runs`.
extract_runs <- function(df, file) {
  if (nrow(df) == 0) {
    return(tibble(file = file, log_run = character(0)))
  }
  d <- df %>% filter(!is.na(.data$Run), .data$Run != "", .data$Run != "StartUp")
  if (nrow(d) == 0) return(tibble(file = file, log_run = character(0)))

  d %>%
    group_by(log_run = .data$Run) %>%
    summarise(
      n_rows       = n(),
      n_trigger    = sum(.data$Description == "PulseReceived", na.rm = TRUE),
      has_trigger  = .data$n_trigger > 0,
      trigger_time = suppressWarnings(
        min(.data$run_start[.data$Description == "PulseReceived"], na.rm = TRUE)),
      n_trials     = sum(.data$Description %in% TRIAL_DESCS, na.rm = TRUE),
      completed    = any(.data$Description %in% DONE_DESCS, na.rm = TRUE),
      dur_s        = suppressWarnings(max(.data$run_time_s, na.rm = TRUE)),
      resp_rate    = suppressWarnings(
        as.numeric(first(na.omit(.data$RespiratoryRate)))),
      .groups = "drop"
    ) %>%
    mutate(
      trigger_time = as.POSIXct(ifelse(is.finite(.data$trigger_time),
                                       .data$trigger_time, NA), tz = "UTC",
                                origin = "1970-01-01"),
      dur_s = ifelse(is.finite(.data$dur_s), .data$dur_s, NA_real_),
      file  = file
    ) %>%
    left_join(INTERO_RUNS, by = "log_run") %>%
    relocate(file, .before = 1)
}

# ---- directory scan -------------------------------------------------------

#' Build the full per-attempt run inventory across a directory of logs.
#'
#' Nothing is dropped here. Empty files, duplicates and pilot records are all
#' represented and flagged, so counts stay reconcilable with the raw directory.
scan_log_dir <- function(dir) {
  files <- sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
  stopifnot(length(files) > 0)

  meta <- parse_log_filename(files) %>%
    mutate(
      path     = files,
      size     = file.size(files),
      md5      = unname(tools::md5sum(files)),
      is_pilot = is_pilot(.data$pid)
    )

  # Byte-identical copies (Drive's " (1)" downloads). Keep the first
  # lexicographically, mark the rest.
  #
  # Zero-byte files are excluded from this test: every empty file has the same
  # md5 (d41d8cd9...), so hashing alone would declare all 18 of them copies of
  # one another. They are not duplicates of anything -- they are launches that
  # died before writing a header -- and are marked is_empty instead.
  # Keep the canonically-named file and mark the " (1)" copy as the duplicate.
  # Ordering on the filename alone gets this backwards: a space (0x20) sorts
  # before a dot (0x2E), so "...340 (1).csv" precedes "...340.csv" and the
  # original would be the one discarded.
  meta <- meta %>%
    mutate(is_empty = .data$size == 0) %>%
    group_by(.data$md5) %>%
    arrange(.data$is_dup_suffix, .data$file, .by_group = TRUE) %>%
    mutate(dup_of = ifelse(.data$is_empty | row_number() == 1,
                           NA_character_, first(.data$file))) %>%
    ungroup() %>%
    mutate(is_duplicate = !is.na(.data$dup_of)) %>%
    arrange(.data$file)

  runs <- map2_dfr(meta$path, meta$file, function(p, f) {
    extract_runs(read_log_file(p), f)
  })

  meta_join <- meta %>%
    select(file, pid, session_type, session_date, launch_time,
           is_pilot, is_empty, is_duplicate, dup_of, size, md5)

  attempts <- runs %>%
    left_join(meta_join, by = "file") %>%
    mutate(
      subject          = pid_bids(.data$pid),
      protocol_version = protocol_version(.data$pid),
      bids_entity      = bids_entities(.data$log_run),
      is_scanner_run   = !is.na(.data$protocol_slot)
    ) %>%
    arrange(.data$pid, .data$session_date, .data$launch_time, .data$protocol_slot)

  list(files = meta, attempts = attempts)
}

# ---- session assembly -----------------------------------------------------

#' Choose, for each (participant, session, run), the attempt that is the run
#' of record.
#'
#' Preference: has a trigger, then completed, then most trials, then the latest
#' launch. That matches the scanner side -- an aborted run is restarted, and the
#' BIDS converter keeps the completed acquisition and drops short fragments.
#' Every attempt is retained with `is_selected` so nothing is hidden.
select_runs <- function(attempts) {
  attempts %>%
    filter(.data$is_scanner_run) %>%
    group_by(.data$pid, .data$session_date, .data$log_run) %>%
    arrange(desc(.data$has_trigger), desc(.data$completed),
            desc(.data$n_trials), .data$launch_time,
            .by_group = TRUE) %>%
    mutate(
      attempt_n   = row_number(),
      n_attempts  = n(),
      is_selected = .data$attempt_n == 1L
    ) %>%
    ungroup()
}

#' One row per participant-session: how much of the 5-run protocol survived.
summarise_sessions <- function(selected) {
  selected %>%
    filter(.data$is_selected) %>%
    group_by(.data$pid, .data$subject, .data$session_date,
             .data$protocol_version) %>%
    summarise(
      n_runs          = n(),
      n_with_trigger  = sum(.data$has_trigger),
      n_completed     = sum(.data$completed),
      runs_present    = paste(sort(.data$protocol_slot), collapse = ","),
      resp_rate       = suppressWarnings(as.numeric(first(na.omit(.data$resp_rate)))),
      n_source_files  = n_distinct(.data$file),
      session_start   = suppressWarnings(min(.data$trigger_time, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      session_start = as.POSIXct(ifelse(is.finite(.data$session_start),
                                        .data$session_start, NA), tz = "UTC",
                                 origin = "1970-01-01"),
      is_pilot = is_pilot(.data$pid)
    ) %>%
    arrange(.data$pid, .data$session_date)
}
