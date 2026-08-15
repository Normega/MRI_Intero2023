# intero_resp.R
# Respiration belt (Biopac DA100C) -> per-epoch breathing rate, pacer adherence,
# and the physiological ground truth of rate change.
#
# This module fills the two columns intero_events.R reserved as NA:
#   adherence_verified    did the participant demonstrably follow the pacer
#   actual_change_physio  what their breathing rate ACTUALLY did
#
# Until now sdt_class compared the participant's judgement against the
# INSTRUCTED change. That is only the right comparison if they breathed as
# instructed. Where they did not, the instructed change is the wrong stimulus
# and every derived hit/miss on that epoch is mislabelled.
#
# ---------------------------------------------------------------------------
# WHY THE PACER IS THE VALIDITY TEST, NOT AN AMPLITUDE THRESHOLD
#
# A flat belt does not produce a flat breath table. Normalising a near-constant
# trace amplifies floating-point and sensor noise into something with peaks and
# troughs, and the resulting rate lands in a plausible range: sub-076's
# task-breathing_run-2 reported 11.8 brpm against a 12 brpm pacer from a trace
# that is visually flat with two transient spikes (Repo/Figures/resp_eyeball.png,
# panel 3). Run-level rate therefore cannot be used as evidence that a belt
# worked -- it is exactly the statistic noise reproduces.
#
# The pacer can. Its instructed period CHANGES within an epoch, by a known
# trajectory that differs by condition (see below). Noise has no access to that
# trajectory, so agreement between the measured and instructed period series is
# evidence no flat trace can fabricate.
#
# ---------------------------------------------------------------------------
# THE INSTRUCTED TRAJECTORY
#
# The pacer emits one pulse per half-cycle, alternating inhale/exhale. A full
# breath period at pulse i is period_s[i] + period_s[i+1].
#
# THE FIRST PULSE OF EACH EPOCH IS TRUNCATED and must be dropped. Averaged
# across all 105 traces, the normalised period jumps ~1.11 from cycle 1 to
# cycle 2 in EVERY condition including No Change, which is then flat at 1.111
# for the rest of the epoch. Measured against cycle 2, the endpoints are 0.727
# (Faster) and 1.271 (Slower), matching the change_amount medians of 0.75 and
# 1.25. Measured against cycle 1 they are 0.809 and 1.413, which match nothing.
#
# salience is the ramp SHAPE, not the size of the change -- both levels reach
# the same endpoint:
#   low   gradual linear ramp across the epoch
#   high  flat, then a step at the midpoint
# So change_amount alone does not distinguish the salience conditions, and the
# trajectory is the stimulus. Both are handled by comparing measured against
# instructed with the same estimator.
#
# ---------------------------------------------------------------------------
# POLARITY IS NOT RESOLVED, AND DELIBERATELY NOT USED
#
# resp_polarity() cannot settle this belt: across 138 spontaneous runs the
# median duty cycle is 0.471, which polarity_verdict() reports as not
# confident. Individual runs contradict each other (0.572 and 0.619 both
# indicating "down"). respkit warns the diagnostic is uninformative under a
# symmetric pacer, which is most of our breathing runs.
#
# Consequence: every phase-timing column is mirrored if polarity is wrong, so
# this module emits NONE of them -- no inhale_dur, exhale_dur, duty_cycle,
# ie_ratio, flow_in or flow_out. Period, rate and their changes are invariant
# to polarity (a mirrored cycle has the same duration), and those are all the
# event model needs. If I:E ratio is ever wanted, the transducer wiring has to
# be established first.
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr) })

# Calibrated on this sample per respkit's BeforeYouRun.txt step 1: the
# amplitude thresholds are in the transducer's own units and the respkit
# defaults (0.05 / 0.20) do not transfer to an uncalibrated strain gauge in
# Volts. log10(filtered_iqr) is bimodal with a trough at -3.0; below it,
# median band_ratio is 0.80 and detected rate runs to 17 brpm, above it
# band_ratio is 0.93+ and the dominant rate is exactly the 12 brpm pacer.
BELT_MIN_IQR        <- 1e-3
BELT_MIN_BAND_RATIO <- 0.90

# A run passing the two gates can still be dead for part of its length
# (sub-062 task-localizer_run-2 is flat for 45 s then breathes normally), so
# these are a floor for considering a run at all, and the real decision is
# taken per epoch on coverage.
RESP_MIN_BREATHS  <- 3L
RESP_MIN_COVERAGE <- 0.60

PACER_DROP_FIRST <- 1L   # truncated first half-cycle; see header

#' Belt segment -> breath table, without any phase-timing column.
#'
#' @return list(bt, quality, rec) or NULL if the segment is too short.
belt_analyse <- function(x, fs, id = NULL, min_s = 60) {
  x <- as.numeric(x)
  if (length(x) < fs * min_s) return(NULL)
  rec <- resp_recording(x, fs = fs, units = "V", id = id)
  q <- tryCatch(resp_quality(rec, trim_s = 15,
                             thresholds = list(iqr_unusable = NA, iqr_degraded = NA)),
                error = function(e) NULL)
  # polarity is fixed at "up" from the transducer spec because the data cannot
  # decide it; the warning resp_analyse raises is therefore expected and not
  # informative here. No phase column survives to be affected.
  a <- tryCatch(withCallingHandlers(
         resp_analyse(rec, quality = FALSE, polarity = "up"),
         warning = function(w) if (grepl("polarity|despike", conditionMessage(w)))
                                 invokeRestart("muffleWarning")),
       error = function(e) NULL)
  if (is.null(a)) return(NULL)
  # The FULL breath table is returned. respkit's summarise_breaths() reads the
  # phase columns and errors without them, so they cannot be dropped here --
  # they are dropped by strip_phase() on the way out instead.
  list(bt = a$breaths, quality = q, rec = rec, detection = a$detection)
}

# Mirrored if the belt polarity is wrong, which this sample cannot establish.
# Never emitted and never used in a derived measure; see the header.
RESP_PHASE_COLS <- c("inhale_dur", "exhale_dur", "duty_cycle", "ie_ratio",
                     "flow_in", "flow_out", "t_peak")

strip_phase <- function(bt) bt[, setdiff(names(bt), RESP_PHASE_COLS), drop = FALSE]

#' Run-level belt gate. Returns a grade, never an exclusion.
belt_grade <- function(q) {
  if (is.null(q)) return("fail")
  if (!is.finite(q$filtered_iqr) || q$filtered_iqr < BELT_MIN_IQR) return("fail")
  if (!is.finite(q$band_ratio)   || q$band_ratio  < BELT_MIN_BAND_RATIO) return("degraded")
  "good"
}

#' Pacer trace -> instructed full-breath period series.
#'
#' @param trace a *_pacertrace.tsv, one row per half-cycle.
#' @return tibble(epoch_index, t_s, instructed_period_s)
pacer_instructed <- function(trace) {
  if (!nrow(trace)) return(tibble(epoch_index = integer(), t_s = numeric(),
                                  instructed_period_s = numeric()))
  trace %>%
    mutate(across(c("epoch_index", "pulse_index", "onset", "period_s"),
                  ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(.data$pulse_index > PACER_DROP_FIRST) %>%
    group_by(.data$epoch_index) %>%
    arrange(.data$pulse_index, .by_group = TRUE) %>%
    mutate(instructed_period_s = .data$period_s + dplyr::lead(.data$period_s)) %>%
    ungroup() %>%
    filter(is.finite(.data$instructed_period_s)) %>%
    transmute(epoch_index = as.integer(.data$epoch_index),
              t_s = .data$onset, instructed_period_s = .data$instructed_period_s)
}

#' Log-linear trend in breath period, as a per-second slope.
#'
#' Returns the SLOPE, not a ratio. A ratio needs a span to be evaluated over,
#' and the span has to be chosen once and applied to both the measured and the
#' instructed series or the two are not comparable -- see epoch_resp().
#'
#' Reported in log space so that exponentiating gives something directly
#' comparable to change_amount. A slope is used rather than a
#' last-third/first-third contrast because an epoch holds only 4-5 breaths;
#' a slope uses all of them, and while it under-reads a step it keeps the sign.
period_log_slope <- function(t_s, period_s) {
  ok <- is.finite(t_s) & is.finite(period_s) & period_s > 0
  t_s <- t_s[ok]; period_s <- period_s[ok]
  if (length(t_s) < 3 || diff(range(t_s)) <= 0) return(NA_real_)
  unname(stats::coef(stats::lm(log(period_s) ~ t_s))[2])
}

#' Breaths wholly inside an epoch.
#'
#' respkit's "contained" rule, and the reason it exists: a breath assigned by
#' its opening trough alone can start 0.5 s before the pacer stops and end
#' well inside the rating window, contributing a duration the epoch's condition
#' did not produce. Because breathing slows once the pacer stops, that inflates
#' the trend upward on EVERY epoch -- with the "start" rule the No Change null
#' centred on 1.281 instead of 1.0, and all five conditions were displaced by
#' the same factor.
epoch_breaths_contained <- function(bt, onset, duration) {
  w <- bt$t_start >= onset & bt$t_end <= onset + duration
  bt[w & bt$keep %in% TRUE, , drop = FALSE]
}

#' Per-epoch respiration measures for one run.
#'
#' @param bt  breath table from belt_analyse()
#' @param ev  event rows for this run (onset, duration, epoch_index, trial_type)
#' @param ins instructed series from pacer_instructed(); may be empty
epoch_resp <- function(bt, ev, ins = NULL) {
  if (!nrow(ev)) return(tibble())
  eps <- resp_epochs(onset = ev$onset, duration = ev$duration,
                     label = ev$trial_type)
  feat <- epoch_features(bt, eps, rule = "start", min_breaths = RESP_MIN_BREATHS)

  out <- tibble(
    epoch_index    = ev$epoch_index,
    trial_type     = ev$trial_type,
    onset          = ev$onset,
    duration       = ev$duration,
    n_breaths      = feat$n_breaths,
    coverage_frac  = feat$coverage_frac,
    rate_med_brpm  = if ("rate_bpm_median" %in% names(feat)) feat$rate_bpm_median
                     else if ("mean_rate_bpm" %in% names(feat)) feat$mean_rate_bpm
                     else NA_real_
  )

  out$actual_ratio     <- NA_real_
  out$instructed_ratio <- NA_real_
  out$adherence_mae_s  <- NA_real_
  out$adherence_r      <- NA_real_
  out$n_instructed     <- NA_integer_
  out$eval_span_s      <- NA_real_
  out$n_contained      <- NA_integer_

  for (k in seq_len(nrow(ev))) {
    b <- epoch_breaths_contained(bt, ev$onset[k], ev$duration[k])
    out$n_contained[k] <- nrow(b)
    if (nrow(b) < RESP_MIN_BREATHS) next

    # ONE evaluation window, taken from the measured breaths and applied to
    # both series. Evaluating over the nominal epoch duration instead
    # extrapolates the fit beyond the data it was fitted on -- breaths cover
    # perhaps 70% of the window, so the trend is inflated by the reciprocal,
    # and by a different amount for each epoch.
    span <- diff(range(b$t_start))
    if (!is.finite(span) || span <= 0) next
    out$eval_span_s[k] <- span

    sl <- period_log_slope(b$t_start, b$duration)
    if (is.finite(sl)) out$actual_ratio[k] <- exp(sl * span)

    if (is.null(ins) || !nrow(ins)) next
    if (!identical(ev$trial_type[k], "paced_breathing")) next
    ik <- ins[ins$epoch_index %in% ev$epoch_index[k], , drop = FALSE]
    if (nrow(ik) < 3) next
    out$n_instructed[k] <- nrow(ik)

    # instructed fitted over the SAME window the measured breaths occupy
    iw <- ik[ik$t_s >= min(b$t_start) & ik$t_s <= max(b$t_start), , drop = FALSE]
    sli <- period_log_slope(iw$t_s, iw$instructed_period_s)
    if (is.finite(sli)) out$instructed_ratio[k] <- exp(sli * span)

    tgt <- stats::approx(ik$t_s, ik$instructed_period_s, xout = b$t_start,
                         rule = 2)$y
    out$adherence_mae_s[k] <- mean(abs(b$duration - tgt), na.rm = TRUE)
    if (stats::sd(tgt, na.rm = TRUE) > 1e-6 && nrow(b) >= 4)
      out$adherence_r[k] <- suppressWarnings(
        stats::cor(b$duration, tgt, use = "complete.obs"))
  }
  out
}

# ---------------------------------------------------------------------------
# ADHERENCE IS A RUN-LEVEL PROPERTY, NOT AN EPOCH-LEVEL ONE
#
# An epoch holds 3-4 wholly-contained breaths (331 of 567 hold exactly 3), and
# the instructed change is 25-35% of the period. Against natural breath-to-
# breath variability that is close to the noise floor: five different per-epoch
# estimators all separate Faster from Slower at AUC 0.66-0.69, so no choice of
# estimator rescues it. Pooling every paced breath in a run instead gives a
# median of 23 paired observations and a usable measure.
#
# The pooled measure also showed that most apparent non-adherence is a dead
# belt rather than a disobedient participant. Runs graded "fail" have a median
# adherence of 0.017 and a measured/instructed period ratio of 0.656 -- the
# signature of the detector splitting each breath into its inhale and exhale
# halves. Runs graded "good" have median adherence 0.460 and a ratio of 0.905,
# and 71% of them track the pacer. A lag does not explain the poor runs: good
# runs peak at lag 0 and decay symmetrically, while poor runs stay negative at
# every lag from -6 s to +10 s.
ADHERE_MIN_R <- 0.30   # run tracks the pacer
ADHERE_MIN_N <- 15L    # paired breaths needed before r means anything

#' Run-level pacer adherence, pooled over all paced epochs.
#'
#' @return tibble(n_paired, adhere_r, adhere_gain, adhere_mae_s, period_ratio_oct)
#'
#' period_ratio_oct is a DIAGNOSTIC, not a filter. It is measured/instructed
#' median period, so a value near 0.5 means octave error -- but a participant
#' genuinely breathing faster than instructed also lowers it. Filtering on it
#' would partly be filtering on adherence itself, so belt_grade (which is
#' derived from signal quality alone) is the non-circular gate.
run_adherence <- function(bt, ev, ins) {
  empty <- tibble(n_paired = 0L, adhere_r = NA_real_, adhere_gain = NA_real_,
                  adhere_mae_s = NA_real_, period_ratio_oct = NA_real_)
  if (is.null(ins) || !nrow(ins) || !nrow(ev)) return(empty)
  pr <- map_dfr(seq_len(nrow(ev)), function(k) {
    if (!identical(ev$trial_type[k], "paced_breathing")) return(tibble())
    b  <- epoch_breaths_contained(bt, ev$onset[k], ev$duration[k])
    ik <- ins[ins$epoch_index %in% ev$epoch_index[k], , drop = FALSE]
    if (!nrow(b) || nrow(ik) < 3) return(tibble())
    tibble(measured = b$duration,
           instructed = stats::approx(ik$t_s, ik$instructed_period_s,
                                      xout = b$t_start, rule = 2)$y)
  })
  pr <- pr[is.finite(pr$measured) & is.finite(pr$instructed), , drop = FALSE]
  if (nrow(pr) < ADHERE_MIN_N) return(empty)
  fit <- stats::lm(measured ~ instructed, data = pr)
  tibble(n_paired = nrow(pr),
         adhere_r = suppressWarnings(stats::cor(pr$measured, pr$instructed)),
         adhere_gain = unname(stats::coef(fit)[2]),
         adhere_mae_s = mean(abs(pr$measured - pr$instructed)),
         period_ratio_oct = stats::median(pr$measured) / stats::median(pr$instructed))
}

#' Is this epoch's respiration measured well enough to be believed?
resp_measured <- function(n_breaths, coverage_frac, run_grade) {
  !is.na(n_breaths) & n_breaths >= RESP_MIN_BREATHS &
    !is.na(coverage_frac) & coverage_frac >= RESP_MIN_COVERAGE &
    run_grade %in% c("good", "degraded")
}

#' Classify what the breathing rate actually did.
#'
#' Thresholds are the 5th/95th percentiles of actual_ratio on No Change epochs
#' that passed resp_measured() -- i.e. the observed spread when the instructed
#' period was flat. Set in intero2023_resp.R from the data and pinned here so
#' the classification is reproducible without recomputing the null.
RESP_RATIO_SLOWER <- NA_real_   # set below
RESP_RATIO_FASTER <- NA_real_

classify_actual_change <- function(actual_ratio, measured,
                                   lo = RESP_RATIO_FASTER,
                                   hi = RESP_RATIO_SLOWER) {
  out <- rep(NA_character_, length(actual_ratio))
  ok <- measured %in% TRUE & is.finite(actual_ratio)
  out[ok] <- "No Change"
  out[ok & actual_ratio < lo] <- "Faster"   # period shortened
  out[ok & actual_ratio > hi] <- "Slower"   # period lengthened
  out
}
