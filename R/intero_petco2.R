# intero_petco2.R
# End-tidal CO2 (PetCO2) extraction from the Biopac capnograph channel.
#
# ---------------------------------------------------------------------------
# WHY THIS DOES NOT USE respkit's BREATH TABLE
#
# respkit's detection and breath_table() are built for a strain-gauge belt: a
# roughly sinusoidal signal where peaks and troughs are symmetric and
# inspiratory/expiratory times are read from the zero crossings. A capnogram is
# not that shape. Each breath is a sawtooth -- near-zero inspiratory baseline,
# a steep expiratory upstroke, an alveolar plateau, then a near-vertical drop as
# the next inhalation flushes the sensor. Peak *detection* transfers; the
# morphology model does not.
#
# So peaks are found with respkit's prominence detector (which is shape-agnostic)
# and PetCO2 is taken as the plateau value directly, rather than routed through
# breath_table's inspiratory/expiratory decomposition, which would report
# meaningless duty cycles for this signal.
#
# PetCO2 = the CO2 concentration at the end of the alveolar plateau, i.e. the
# per-breath maximum. Baseline (the inter-breath trough) should sit near zero if
# the mask is sealed; a rising baseline means rebreathing or a leak, so it is
# retained as a quality measure rather than subtracted away.
#
# UNITS. The channel is already %CO2. To convert to the mmHg conventionally
# reported, multiply by (barometric pressure - 47) / 100; at Toronto sea level
# roughly x7.1. Not applied here -- the raw % is kept and the conversion left
# explicit downstream.
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr) })

PETCO2_MIN_BREATH_S <- 2.0   # fastest plausible breath: 16 brpm x0.7 pacing
PETCO2_MAX_BREATH_S <- 20.0  # beyond this a "breath" is a detection failure
PETCO2_PLAUSIBLE    <- c(1.0, 9.0)   # %CO2; outside this is not an alveolar plateau

#' Detect end-tidal peaks in a capnogram.
#'
#' @param x  CO2 signal, %CO2
#' @param fs sampling rate
#' @param min_prominence absolute %CO2; NULL derives it from the signal so a
#'   session recorded at low gain is not silently rejected (capnogram amplitude
#'   varies several-fold across this sample).
#' @return tibble(breath, t_s, petco2, baseline, period_s)
petco2_detect <- function(x, fs, min_prominence = NULL,
                          min_breath_s = PETCO2_MIN_BREATH_S) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < fs * 10) return(tibble())

  # A FIXED threshold, not one scaled to the signal. The adaptive version
  # (0.4 x the p10-p90 swing) silently merged adjacent breaths on
  # large-amplitude recordings: sub-130's breathing run reported 5.5 brpm
  # against a paced 12, because a healthy capnogram's own swing pushed the
  # threshold above the breath-to-breath variation. Verified against
  # hand-counted traces: 0.2 %CO2 recovers 12-13 peaks per 60 s at 12 brpm
  # pacing on both clean (sub-100) and noisy (sub-130) recordings.
  if (is.null(min_prominence)) min_prominence <- PETCO2_MIN_PROMINENCE

  pk <- try(find_extrema_prominence(x, fs = fs,
                                    min_prominence = min_prominence,
                                    min_distance_s = min_breath_s),
            silent = TRUE)
  idx <- if (inherits(pk, "try-error")) integer(0)
         else if (is.list(pk)) (pk$peaks %||% pk[[1]]) else pk
  idx <- as.integer(idx[is.finite(idx)])
  idx <- idx[idx >= 1 & idx <= n]
  if (length(idx) < 3) return(tibble())

  # Baseline for each breath: the minimum between this peak and the previous
  # one. A sealed mask returns to ~0; a drifting floor is a leak or rebreathing.
  base <- vapply(seq_along(idx), function(i) {
    lo <- if (i == 1) 1L else idx[i - 1]
    min(x[lo:idx[i]], na.rm = TRUE)
  }, numeric(1))

  tibble(
    breath   = seq_along(idx),
    t_s      = (idx - 1) / fs,
    petco2   = x[idx],
    baseline = base,
    period_s = c(NA_real_, diff(idx) / fs)
  )
}

#' Flag implausible breaths without deleting them.
#'
#' Kept rather than dropped so the rejection rate is itself a quality measure:
#' a run where 40% of detections are implausible is a bad recording, and that
#' is more useful to know than a silently shortened series.
petco2_flag <- function(pk) {
  if (!nrow(pk)) return(pk)
  pk %>% mutate(
    ok_value  = .data$petco2 >= PETCO2_PLAUSIBLE[1] & .data$petco2 <= PETCO2_PLAUSIBLE[2],
    ok_period = is.na(.data$period_s) |
                (.data$period_s >= PETCO2_MIN_BREATH_S & .data$period_s <= PETCO2_MAX_BREATH_S),
    ok = .data$ok_value & .data$ok_period
  )
}

#' Resample per-breath PetCO2 onto the fMRI volume grid.
#'
#' PetCO2 is a per-breath measurement, so it is a sparse irregular series that
#' must be interpolated to be a regressor. Linear interpolation between
#' consecutive end-tidal values is the standard choice and is what the
#' respiratory-response-function literature assumes as input.
#'
#' Returns one value per volume. Volumes outside the detected range are NA
#' rather than extrapolated -- a run whose capnogram dies partway should show
#' the gap, not a flat line carried to the end.
petco2_to_volumes <- function(pk, n_vol, tr = 2, run_onset_s = 0) {
  t_vol <- run_onset_s + (seq_len(n_vol) - 1) * tr
  # petco2_flag() returns its (empty) input untouched when nothing was
  # detected, so `ok` may not exist at all -- a run with no capnogram must
  # yield an all-NA regressor, not an error that stops the whole cohort.
  if (!nrow(pk) || !"ok" %in% names(pk)) return(rep(NA_real_, n_vol))
  good <- pk %>% filter(.data$ok)
  if (nrow(good) < 3) return(rep(NA_real_, n_vol))
  approx(x = good$t_s, y = good$petco2, xout = t_vol,
         method = "linear", rule = 1)$y
}

#' Per-run PetCO2 quality summary.
petco2_quality <- function(pk, dur_s) {
  if (!nrow(pk) || !"ok" %in% names(pk)) {
    return(tibble(n_breaths = 0L, rate_brpm = NA_real_, petco2_med = NA_real_,
                  petco2_iqr = NA_real_, baseline_med = NA_real_,
                  frac_ok = NA_real_))
  }
  g <- pk %>% filter(.data$ok)
  tibble(
    n_breaths    = nrow(pk),
    rate_brpm    = if (dur_s > 0) 60 * nrow(g) / dur_s else NA_real_,
    petco2_med   = median(g$petco2, na.rm = TRUE),
    petco2_iqr   = IQR(g$petco2, na.rm = TRUE),
    baseline_med = median(g$baseline, na.rm = TRUE),
    frac_ok      = mean(pk$ok)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

PETCO2_MIN_PROMINENCE <- 0.2   # %CO2; see petco2_detect()

# ---------------------------------------------------------------------------
# CAPNOGRAM VALIDITY
#
# co2_present (a high p99) is necessary but NOT sufficient. Six participants
# sit at 3.6-4.6 %CO2 that never returns toward zero -- the sensor never sees
# fresh air -- yet all six passed co2_present with healthy p99 values, and three
# of them even produced a plausible breath COUNT. Only the baseline separates
# them. Scattered across Feb-May 2024 with valid sessions on the same days
# (sub-114 baseline 4.57 vs sub-116 0.50, both 2024-05-21), so this is per-session
# setup and cannot be excluded by date the way the pre-2024-01-25 capnograph
# outage can.
#
# Computed over the WHOLE run. A single 60 s window misrepresented sub-132,
# which returned 1 peak from a window that landed in a gap while its same-day
# neighbour sub-130 was fine.
# ---------------------------------------------------------------------------

CAPNO_MAX_BASELINE <- 2.0   # %CO2 the trough must fall below
CAPNO_MIN_SWING    <- 0.8   # %CO2 peak-to-trough
CAPNO_RATE_RANGE   <- c(6, 25)   # brpm

capnogram_validity <- function(pk, dur_s) {
  if (!nrow(pk) || dur_s <= 0)
    return(tibble(baseline_med = NA_real_, swing_med = NA_real_,
                  rate_brpm = NA_real_, capnogram_valid = FALSE,
                  fail_reason = "no peaks detected"))
  base <- median(pk$baseline, na.rm = TRUE)
  sw   <- median(pk$petco2 - pk$baseline, na.rm = TRUE)
  rate <- 60 * nrow(pk) / dur_s
  reasons <- c(
    if (!is.finite(base) || base >= CAPNO_MAX_BASELINE) "baseline never flushes",
    if (!is.finite(sw)   || sw   <= CAPNO_MIN_SWING)    "swing too small",
    if (!is.finite(rate) || rate < CAPNO_RATE_RANGE[1] ||
        rate > CAPNO_RATE_RANGE[2])                     "implausible rate"
  )
  tibble(baseline_med = base, swing_med = sw, rate_brpm = rate,
         capnogram_valid = length(reasons) == 0,
         fail_reason = if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_)
}
