# intero_physio.R
# Biopac AcqKnowledge (.acq) ingestion, run segmentation and coverage QC.
#
# ---------------------------------------------------------------------------
# WHAT IS IN AN .acq  (verified on P130/P100/P85, 2026-08-06)
#
# Three channels at 2000 Hz, ~45-60 min per participant:
#   CO2100C          %CO2   capnograph -- the PetCO2 source
#   Custom, DA100C   Volts  respiration belt
#   Scanner Pulse    Volts  0/5 V TTL, ONE PULSE PER VOLUME
#
# The trigger channel segments the session into 10 blocks with an inter-pulse
# interval of exactly 2.000 s (= TR):
#
#   196, 5, 196, 5, 242, 5, 242, 5, 172, 5
#
# THE 5-PULSE BLOCKS ARE THE dir-PA FIELDMAP BLIPS, NOT ABORTED RUNS. The BIDS
# converter's ABORT_MAX_VOLS = 100 rule is about BOLD series and must not be
# applied here -- it would discard every fieldmap. Classification is three-way.
#
# TRIGGERS ARE DETECTED AT THE NATIVE RATE, BEFORE ANY DECIMATION. A TTL pulse
# is only ~10 ms wide (about 20 samples at 2000 Hz). Decimating to 100 Hz first
# would leave roughly one sample per pulse and drop an unpredictable fraction of
# them depending on phase, silently corrupting the volume counts that the whole
# alignment rests on.
#
# COVERAGE IS PER RUN, NOT PER PARTICIPANT. Participants took the mouthpiece out
# between scans, so CO2 and belt drop out for entire runs -- most often the
# authenticity run, which needs no paced breathing. Every run is gated on its
# own.
#
# BELT AMPLITUDE IS NOT COMPARABLE BETWEEN PARTICIPANTS (gain differed by ~20x
# across sessions), so no fixed amplitude threshold is applied here. This module
# only *measures*; the cut is chosen later from the sample-wide distribution, as
# respkit's BeforeYouRun.txt instructs.
#
# Author: RADlab

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(purrr)
})

ACQ_FS            <- 2000     # native sampling rate, Hz
ACQ_TRIG_THRESH   <- 2.5      # V, midpoint of the 0/5 V TTL
ACQ_RUN_GAP_S     <- 5        # a gap this long between pulses starts a new block
ACQ_FIELDMAP_N    <- 5L       # pulses in a dir-PA blip
ACQ_MIN_BOLD_N    <- 100L     # below this (and not a fieldmap) = aborted fragment
PHYSIO_TARGET_FS  <- 100      # decimation target; 2000/100 = 20, exact

#' Locate a participant's .acq. Filenames are unpadded (Intero_P52.acq),
#' unlike everything else in the project.
acq_path <- function(pid, dir_physio) {
  n <- pid_num(pid)
  cand <- file.path(dir_physio, sprintf("Intero_P%d.acq", n))
  if (file.exists(cand)) return(cand)
  # tolerate zero-padded variants if they ever appear
  alt <- file.path(dir_physio, sprintf("Intero_P%03d.acq", n))
  if (file.exists(alt)) return(alt)
  NA_character_
}

#' Rising edges of the TTL trigger, in samples. Native rate only.
trigger_edges_ttl <- function(trig, threshold = ACQ_TRIG_THRESH) {
  hi <- trig > threshold
  which(hi & !c(FALSE, hi[-length(hi)]))
}

#' Segment trigger edges into acquisition blocks.
#'
#' @return tibble: block, onset_s, end_s, dur_s, n_pulses, median_iti_s, kind
segment_trigger_blocks <- function(edges, fs = ACQ_FS, gap_s = ACQ_RUN_GAP_S) {
  if (!length(edges)) {
    return(tibble(block = integer(0), onset_s = numeric(0), end_s = numeric(0),
                  dur_s = numeric(0), n_pulses = integer(0),
                  median_iti_s = numeric(0), kind = character(0)))
  }
  gaps <- diff(edges) / fs
  brk  <- c(0L, which(gaps > gap_s), length(edges))
  map_dfr(seq_len(length(brk) - 1L), function(i) {
    ix <- edges[(brk[i] + 1L):brk[i + 1L]]
    n  <- length(ix)
    tibble(
      block        = i,
      onset_s      = ix[1] / fs,
      end_s        = ix[n] / fs,
      dur_s        = (ix[n] - ix[1]) / fs,
      n_pulses     = n,
      median_iti_s = if (n > 1) median(diff(ix) / fs) else NA_real_,
      kind = if (n == ACQ_FIELDMAP_N) "fieldmap"
             else if (n >= ACQ_MIN_BOLD_N) "bold"
             else "fragment"
    )
  })
}

#' Attach BIDS run identities to the BOLD blocks.
#'
#' Blocks are matched to the participant's expected volume counts in acquisition
#' order. An exact pulse/volume match is required to call a block canonical --
#' order alone would mislabel sessions that skipped a run (sub-082 and sub-104
#' acquired slots 1, 3, 5 only).
label_bold_blocks <- function(blocks, expected) {
  bold <- blocks %>% filter(kind == "bold") %>% arrange(onset_s)
  if (!nrow(bold) || !nrow(expected)) {
    return(blocks %>% mutate(protocol_slot = NA_integer_,
                             bids_entity = NA_character_,
                             n_volumes = NA_integer_, pulse_match = NA))
  }
  exp <- expected %>% arrange(protocol_slot)
  slot <- rep(NA_integer_, nrow(bold)); vols <- rep(NA_integer_, nrow(bold))
  used <- rep(FALSE, nrow(exp))
  for (i in seq_len(nrow(bold))) {
    hit <- which(!used & exp$n_volumes == bold$n_pulses[i])
    if (length(hit)) {
      j <- hit[1]; used[j] <- TRUE
      slot[i] <- exp$protocol_slot[j]; vols[i] <- exp$n_volumes[j]
    }
  }
  bold$protocol_slot <- slot
  bold$n_volumes     <- vols
  bold$pulse_match   <- !is.na(slot)
  bold$bids_entity   <- bids_entities(
    INTERO_RUNS$log_run[match(slot, INTERO_RUNS$protocol_slot)])

  blocks %>%
    left_join(bold %>% select(block, protocol_slot, n_volumes, pulse_match,
                              bids_entity), by = "block")
}

#' Per-block signal quality. Measures only; applies no amplitude cut.
#'
#' `co2_p99` is the practical PetCO2 proxy: a mask that is on sits around
#' 4-6 %CO2, one that has been removed collapses to ~0.5.
block_signal_stats <- function(blocks, co2, belt, fs = ACQ_FS, pad_s = 2) {
  n <- length(co2)
  # An explicitly-typed empty frame, not map_dfr's zero-COLUMN result. A
  # participant whose session was abandoned before any scanning (sub-110 has an
  # .acq but no imaging) yields no trigger blocks, and a zero-column tibble
  # makes the downstream join fail with "Join columns in `y` must be present".
  empty <- tibble(block = integer(0), co2_p99 = numeric(0), co2_med = numeric(0),
                  co2_sd = numeric(0), belt_sd = numeric(0),
                  belt_iqr = numeric(0), belt_rng = numeric(0))
  if (!nrow(blocks)) return(empty)
  map_dfr(seq_len(nrow(blocks)), function(i) {
    s <- max(1L, round(blocks$onset_s[i] * fs))
    e <- min(n, round(blocks$end_s[i] * fs) + pad_s * fs)
    if (e <= s) return(tibble(block = blocks$block[i]))
    c2 <- co2[s:e]; bl <- belt[s:e]
    tibble(
      block    = blocks$block[i],
      co2_p99  = unname(quantile(c2, .99, na.rm = TRUE)),
      co2_med  = median(c2, na.rm = TRUE),
      co2_sd   = sd(c2, na.rm = TRUE),
      belt_sd  = sd(bl, na.rm = TRUE),
      belt_iqr = unname(IQR(bl, na.rm = TRUE)),
      belt_rng = diff(range(bl, na.rm = TRUE))
    )
  })
}

#' Read one .acq: triggers at native rate, analog channels decimated.
#'
#' @return list(meta, blocks, co2, belt) where co2/belt are resp_recording
#'   objects at `target_fs`.
read_acq_session <- function(path, id = NULL, target_fs = PHYSIO_TARGET_FS) {
  br <- reticulate::import("bioread")
  d  <- br$read_file(path)
  nm <- vapply(d$channels, function(ch) ch$name, character(1))

  pick <- function(pattern) {
    i <- grep(pattern, nm, ignore.case = TRUE)
    if (length(i)) i[1] else NA_integer_
  }
  i_co2 <- pick("CO2"); i_belt <- pick("DA100C"); i_trig <- pick("Scanner|Pulse|Trig")
  if (any(is.na(c(i_co2, i_belt, i_trig))))
    stop("unexpected channel set in ", basename(path), ": ",
         paste(nm, collapse = " | "))

  fs   <- d$channels[[i_co2]]$samples_per_second
  co2  <- as.numeric(d$channels[[i_co2]]$data)
  belt <- as.numeric(d$channels[[i_belt]]$data)
  trig <- as.numeric(d$channels[[i_trig]]$data)

  # --- triggers first, at the native rate (see header note) ---
  edges  <- trigger_edges_ttl(trig)
  blocks <- segment_trigger_blocks(edges, fs = fs)
  stats  <- block_signal_stats(blocks, co2, belt, fs = fs)
  blocks <- left_join(blocks, stats, by = "block")

  # --- then decimate the analog channels ---
  # Some .acq files carry channels with no samples at all (sub-083, whose
  # session was stopped after run 1). respkit rightly refuses to decimate
  # those, so skip decimation rather than let the whole participant error out
  # and vanish from the coverage table.
  q <- round(fs / target_fs)
  too_short <- length(co2) < 30 * q
  if (too_short) q <- 1L
  co2_d  <- if (q > 1) resp_decimate(co2,  q = q, fs = fs) else co2
  belt_d <- if (q > 1) resp_decimate(belt, q = q, fs = fs) else belt
  if (is.list(co2_d))  co2_d  <- co2_d$signal %||% co2_d[[1]]
  if (is.list(belt_d)) belt_d <- belt_d$signal %||% belt_d[[1]]
  fs_out <- fs / q

  list(
    meta = list(id = id, path = path, fs_native = fs, fs = fs_out, empty = too_short,
                decim_q = q, n_native = length(co2),
                duration_s = length(co2) / fs, channels = nm),
    blocks = blocks,
    # An empty .acq (sub-083's session was stopped after run 1) has no signal to
    # wrap; respkit refuses a zero-length recording, and rightly so. Return NULL
    # channels so the participant still appears in the coverage table as an
    # explicit empty rather than disappearing behind an error.
    co2  = if (length(co2_d))  resp_recording(as.numeric(co2_d),  fs = fs_out,
                                 id = id, units = "%CO2", fs_native = fs) else NULL,
    belt = if (length(belt_d)) resp_recording(as.numeric(belt_d), fs = fs_out,
                                 id = id, units = "Volts", fs_native = fs) else NULL
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# MULTI-PART SESSIONS
#
# Nine participants have their session split across several .acq files. This is
# not an accident of storage: the operators saved as "a", cleared the buffer and
# re-saved as "b" (and "c") to stop the file growing without bound -- the RA's
# session notes describe it for P056 and P058, and P070's restart produced
# "P70_A". Suffixes seen: a/b/c, _1, _A, _PartA.
#
# The parts are CONTINUATIONS, not duplicates: each holds different runs. P058's
# five runs are spread across three files, with part b also capturing two
# aborted authenticity attempts (26 and 47 pulses) before part c got a clean one.
#
# Parts are ordered by `earliest_marker_created_at`, an absolute UTC timestamp
# from the file header -- not by filename suffix, which is a convention rather
# than a guarantee. That timestamp also places every block on one session
# timeline, which is what makes the parts joinable to the logs at all: P058c
# starts 14:02:15.759 local and its authenticity log file is
# ...14h02.15.636.csv, agreeing to within a second.
# ---------------------------------------------------------------------------

ACQ_FILE_RE <- "^Intero_P([0-9]+)((?:[a-c]|_1|_A|_PartA)?)(?: [(][0-9]+[)])?[.]acq$"

#' All .acq parts for every participant, deduplicated and time-ordered.
acq_file_index <- function(dir_physio) {
  fs <- list.files(dir_physio, pattern = "[.]acq$")
  m  <- stringr::str_match(fs, ACQ_FILE_RE)
  ok <- !is.na(m[, 1])
  idx <- tibble(
    file  = fs[ok],
    path  = file.path(dir_physio, fs[ok]),
    pid   = as.integer(m[ok, 2]),
    part  = ifelse(m[ok, 3] == "", "-", m[ok, 3]),
    is_dup_suffix = grepl("[(][0-9]+[)]", fs[ok])
  )
  # Drive's " (1)" copies are byte-identical; drop them exactly as for the logs.
  idx$md5 <- unname(tools::md5sum(idx$path))
  idx %>%
    group_by(.data$md5) %>%
    arrange(.data$is_dup_suffix, .data$file, .by_group = TRUE) %>%
    mutate(is_duplicate = dplyr::row_number() > 1L) %>%
    ungroup() %>%
    filter(!.data$is_duplicate) %>%
    arrange(.data$pid, .data$part)
}

#' Absolute start time of an .acq, on the SAME CLOCK as the PsychoPy logs.
#'
#' The .acq header stores a genuine UTC instant; the logs store bare local
#' wall-clock strings which read_log_file() parses as if they were UTC. Compared
#' directly, every match came out ~18000 s (5 h) adrift.
#'
#' The offset is not a constant: sessions run Sept 2023 - June 2024 and so span
#' both EST (-5) and EDT (-4). Rather than subtract a fixed number, format the
#' UTC instant into Toronto local wall time and re-read it as UTC. That puts it
#' on the logs' naive clock and lets the tz database handle DST.
ACQ_SITE_TZ <- "America/Toronto"

acq_start_time <- function(path, tz = ACQ_SITE_TZ) {
  br <- reticulate::import("bioread")
  d  <- br$read_file(path)
  t  <- try(reticulate::py_get_attr(d, "earliest_marker_created_at"), silent = TRUE)
  if (inherits(t, "try-error") || is.null(t)) return(as.POSIXct(NA))
  t_utc <- as.POSIXct(reticulate::py_str(t), format = "%Y-%m-%d %H:%M:%OS",
                      tz = "UTC")
  if (is.na(t_utc)) return(t_utc)
  as.POSIXct(format(t_utc, tz = tz, usetz = FALSE), tz = "UTC")
}

#' Match BOLD blocks to expected runs with tolerance and monotone slot order.
#'
#' Exact pulse==volume matching fails on the v1 cohort, whose acquisitions were
#' clipped: P058c has 173 pulses against 172 expected volumes, P054 shows 173
#' and 176 against 175. Matching therefore allows a small count difference but
#' still requires slots to advance in acquisition order, so two runs of equal
#' length (loc1 and loc2 are both 190) cannot swap.
match_blocks_tolerant <- function(bold, expected, tol = ACQ_MATCH_TOL) {
  n <- nrow(bold); slot <- rep(NA_integer_, n); vols <- rep(NA_integer_, n)
  d_used <- rep(NA_integer_, n)
  exp <- expected %>% arrange(protocol_slot)
  last <- 0L
  for (i in seq_len(n)) {
    avail <- which(exp$protocol_slot > last)
    if (!length(avail)) break
    dif <- abs(exp$n_volumes[avail] - bold$n_pulses[i])
    j   <- avail[which.min(dif)]
    if (min(dif) <= tol) {
      slot[i] <- exp$protocol_slot[j]; vols[i] <- exp$n_volumes[j]
      d_used[i] <- min(dif); last <- exp$protocol_slot[j]
    }
  }
  bold$protocol_slot <- slot; bold$n_volumes <- vols
  bold$pulse_delta   <- d_used
  bold$pulse_match   <- !is.na(slot)
  bold$bids_entity   <- bids_entities(
    INTERO_RUNS$log_run[match(slot, INTERO_RUNS$protocol_slot)])
  bold
}

ACQ_MATCH_TOL <- 5L   # pulses; v1 acquisitions were clipped by a few volumes

#' Match acq BOLD blocks to log runs by WALLCLOCK, with pulse count as
#' corroboration only.
#'
#' Time is the primary key because it is the one thing that is reliable: the
#' .acq header timestamp and the log's first-trigger wallclock agree to the
#' second wherever both are sound (verified 5/5 on P058 across three files).
#' Pulse counts cannot be primary -- v1 acquisitions were clipped by up to ~20
#' volumes (P054's localizer 2 recorded 173 of 190), which is far beyond any
#' usable count tolerance, and identical-length runs (loc1/loc2 are both 190)
#' are not separable by count at all.
#'
#' @param bold blocks with an absolute `abs_time`
#' @param log_runs the participant's runs of record, with `trigger_time` and
#'   `protocol_slot`
#' @return bold with protocol_slot, bids_entity, dt_s, margin_s, pulse_delta
match_blocks_by_time <- function(bold, log_runs, expected = NULL) {
  n <- nrow(bold)
  out <- bold
  # Every output column is created up front, including on the paths that skip
  # matching. A participant with no imaging at all (sub-076 is absent from
  # inventory.tsv) otherwise reaches the caller missing n_volumes entirely, and
  # the whole run aborts on a participant we already knew had no scans.
  out$protocol_slot <- NA_integer_; out$dt_s <- NA_real_
  out$margin_s <- NA_real_; out$pulse_delta <- NA_integer_
  out$n_volumes <- NA_integer_; out$dup_slot <- FALSE
  if (!n || is.null(log_runs) || !nrow(log_runs)) {
    out$bids_entity <- NA_character_
    return(out)
  }
  lr <- log_runs %>% filter(!is.na(.data$trigger_time))
  if (!nrow(lr)) { out$bids_entity <- NA_character_; return(out) }

  for (i in seq_len(n)) {
    d <- abs(as.numeric(difftime(bold$abs_time[i], lr$trigger_time, units = "secs")))
    o <- order(d)
    out$dt_s[i] <- d[o[1]]
    # margin is what makes the match identifiable, not just close
    out$margin_s[i] <- if (length(d) > 1) d[o[2]] - d[o[1]] else Inf
    out$protocol_slot[i] <- lr$protocol_slot[o[1]]
  }
  out$bids_entity <- bids_entities(
    INTERO_RUNS$log_run[match(out$protocol_slot, INTERO_RUNS$protocol_slot)])

  if (!is.null(expected) && nrow(expected)) {
    ev <- expected$n_volumes[match(out$protocol_slot, expected$protocol_slot)]
    out$n_volumes   <- ev
    out$pulse_delta <- abs(out$n_pulses - ev)
  }
  # A slot claimed by two blocks means one of them is an aborted retry; keep the
  # closer match and mark the other, rather than emitting a duplicate run.
  out$dup_slot <- FALSE
  for (s in unique(na.omit(out$protocol_slot))) {
    ix <- which(out$protocol_slot == s)
    if (length(ix) > 1) {
      keep <- ix[which.min(out$dt_s[ix])]
      out$dup_slot[setdiff(ix, keep)] <- TRUE
    }
  }
  out
}
