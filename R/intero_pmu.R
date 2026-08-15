# intero_pmu.R
# Siemens PMU physiology: the scanner's own respiration belt.
#
# WHY THIS EXISTS ALONGSIDE THE BIOPAC
#
# 21 of 90 Biopac breathing runs grade "fail" -- a flat or octave-erroring belt
# (see intero_resp.R). The scanner logged its own respiratory cushion for the
# same sessions, on independent hardware, and its coverage is better: 258 .resp
# files across 54 sessions, 48 of them with all five runs. It is also easier to
# attribute, because the run number is IN THE FILENAME
# (INTERO23_P052_run3.resp -> protocol slot 3), so none of the wallclock
# matching the Biopac needed applies.
#
# FILE FORMAT (Siemens ASCII PMU, XA30)
#
#   <4 header ints> then a stream of sample groups, terminated by 5003, then a
#   plain-text footer with LogStart/StopMDHTime, LogStart/StopMPCUTime and
#   per-channel Min/Max/Avg.
#
#   5002 ... 6002 delimit free-text comment blocks that appear both in the
#   preamble AND part-way through the data.
#
#   Each sample is followed by four constant flag words
#   (33554432, 98304, 67108864, 163840).
#
# THE FLAGS ARE REMOVED BY VALUE, NOT BY POSITION. Taking every fifth token
# looks correct on the first few hundred samples and then desynchronises,
# because the mid-stream comment blocks do not occupy a whole group: parsing
# P052 run1 positionally returned a "signal" with sd 2.6e7 -- the flag words
# themselves. The signal is 12-bit (0-4095) and the flag words are all >= 98304,
# so the two ranges cannot collide and a value filter is exact.
#
# SAMPLING RATE is 400 Hz, confirmed rather than assumed: 212,890 samples over
# the 532.3 s MDH window is 400.0 Hz, and the same holds for every run checked.
#
# ---------------------------------------------------------------------------
# STATUS: THE DECODE IS VALIDATED. THE ALIGNMENT IS NOT SOLVED.
#
# Validated against the Biopac belt on 50 runs where the Biopac graded "good"
# (intero2023_pmu_validate.R). Searching the whole PMU file for the offset that
# best matches the Biopac trace gives median |r| = 0.617, with 74% of runs above
# 0.5, and the peaks are SHARP -- the best offset beats the next peak more than
# 5 s away by 0.363 on average. Two independent transducers cannot agree that
# well by accident, so the parser is right.
#
# What is NOT solved is where the acquisition starts inside the PMU file:
#
#   * Logging starts BEFORE the scan. Assuming otherwise (take the first
#     n_vol * TR seconds) puts the cross-correlation peak at the -20 s edge of
#     the search for most runs and costs 0.425 of the correlation.
#   * The file's own metadata does not give it. LogStartMDHTime -
#     LogStartMPCUTime ranges from -117 s to +38 s and correlates 0.070 with
#     the empirically recovered offset; at the declared offset, r is 0.024.
#   * There is no trigger to fall back on. The .ext/.ext2 channels are all
#     zeros and the footer reports NrTrig 0. (The footer is a session-wide
#     summary copied verbatim into every file -- .ext2 carries the same "RESP
#     Min Max 1840 8600" as .resp -- which is also why those footer figures
#     never matched the 12-bit sample range.)
#
# The recovered offset is tightly clustered -- median 19.5 s into the file,
# 18.5 s for slot 3 and 20.0 s for slot 4, IQR 6 s, and within a session the
# runs agree to a median of 4 s. So a fixed ~19 s default is roughly right, but
# an IQR of 6 s is 3 TRs, which is too coarse for a per-volume regressor.
#
# CONSEQUENCE: PMU is usable for run-level summaries (rate, variability), and
# for confirming the Biopac. It is NOT yet usable for volume-wise regressors,
# and it therefore cannot currently rescue the 21 dead-Biopac runs -- precisely
# the runs where no Biopac trace exists to align it against. The open option is
# to align via the pacer during breathing runs, which is external to both belts;
# note that would make any subsequent pacer-adherence measure circular.
#
# A SIDE OBSERVATION, UNDERPOWERED: the sign of the PMU-Biopac correlation is
# consistent within a participant for 9 of 13 (6 inverted, 3 matching), which
# is only p = 0.13 on 2 runs each. It hints that the BIOPAC belt polarity
# varies between sessions -- which would explain why resp_polarity() could
# never settle it -- but it does not establish that, and nothing depends on it.
# ---------------------------------------------------------------------------
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(stringr) })

PMU_FS         <- 400          # Hz, verified from MDH/MPCU windows
PMU_FLAGS      <- c(33554432, 98304, 67108864, 163840)
PMU_MAX        <- 4095         # 12-bit full scale
PMU_TARGET_FS  <- 100          # decimation target, matches the Biopac cache
PMU_FILE_RE    <- "^INTERO23_P([0-9]+)_run([0-9]+)[.]resp$"

#' Parse one Siemens PMU .resp file.
#'
#' @return list(signal, fs, n, mdh_start_ms, mdh_stop_ms, mpcu_start_ms,
#'   mpcu_stop_ms, foot_min, foot_max, sat_lo, sat_hi) or NULL if unparseable.
pmu_read <- function(path) {
  if (!file.exists(path) || file.size(path) < 1000) return(NULL)
  raw <- readChar(path, file.size(path), useBytes = TRUE)
  i <- regexpr("5003", raw, fixed = TRUE)
  if (i < 1) return(NULL)
  body <- substr(raw, 1, i - 1)
  foot <- substr(raw, i, nchar(raw))

  tok <- scan(text = body, what = "", quiet = TRUE)
  num <- suppressWarnings(as.numeric(tok))

  # Drop 5002 ... 6002 comment blocks wherever they occur. Vectorised: a
  # running depth from the cumulative counts of the open and close markers.
  # At an opening token the depth is already 1 so it is excluded; at a closing
  # token the depth is back to 0, so that one is excluded explicitly.
  #
  # The blocks cannot simply be ignored: they contain small integers
  # (MSGTYPE 103, currentCushionGain 47, MSGTYPE 400) that fall inside the
  # 12-bit signal range and would be indistinguishable from samples.
  is5002 <- !is.na(num) & num == 5002
  is6002 <- !is.na(num) & num == 6002
  depth <- cumsum(is5002) - cumsum(is6002)
  v <- num[depth == 0 & !is6002]
  v <- v[!is.na(v)]
  if (length(v) < 100) return(NULL)
  v <- v[-(1:4)]                       # header
  sig <- v[!(v %in% PMU_FLAGS)]        # by value; see header
  sig <- sig[sig >= 0 & sig <= PMU_MAX]
  if (length(sig) < PMU_FS * 30) return(NULL)

  gnum <- function(pat) {
    m <- regmatches(foot, regexpr(paste0(pat, "[^\n]*"), foot))
    if (!length(m)) return(NA_real_)
    x <- as.numeric(regmatches(m, gregexpr("-?[0-9]+", m))[[1]])
    if (!length(x)) NA_real_ else x[length(x)]
  }
  fmm <- regmatches(foot, regexpr("RESP Min Max Avg StdDiff:[^\n]*", foot))
  fv <- if (length(fmm)) as.numeric(regmatches(fmm, gregexpr("-?[0-9]+", fmm))[[1]]) else rep(NA_real_, 4)

  list(signal = sig, fs = PMU_FS, n = length(sig),
       mdh_start_ms = gnum("LogStartMDHTime"), mdh_stop_ms = gnum("LogStopMDHTime"),
       mpcu_start_ms = gnum("LogStartMPCUTime"), mpcu_stop_ms = gnum("LogStopMPCUTime"),
       foot_min = fv[1], foot_max = fv[2],
       sat_lo = mean(sig == 0), sat_hi = mean(sig == PMU_MAX))
}

#' Index every PMU .resp file under Data/MRI, with participant and run slot.
pmu_file_index <- function(dir_mri) {
  f <- list.files(dir_mri, pattern = "[.]resp$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(tibble())
  bn <- basename(f)
  m <- str_match(bn, PMU_FILE_RE)
  # dplyr:: qualified: the signal package masks filter() and is loaded by the
  # scripts that call this, so a bare filter() resolves to signal::filter.
  tibble(path = f, file = bn,
         pid = suppressWarnings(as.integer(m[, 2])),
         protocol_slot = suppressWarnings(as.integer(m[, 3]))) %>%
    dplyr::filter(!is.na(.data$pid), !is.na(.data$protocol_slot))
}

#' Take the acquisition window, assuming logging starts with the scan.
#'
#' Returns the first n_vol * tr seconds. When the log is SHORTER than the
#' acquisition the run is returned truncated and `covered` says so, rather than
#' being padded to length -- a short log is a fact about the session, and
#' padding it would hide that.
pmu_window <- function(sig, fs, n_vol, tr = 2) {
  need <- round(n_vol * tr * fs)
  n <- min(need, length(sig))
  list(signal = sig[seq_len(n)], covered = n / need, n_needed = need)
}
