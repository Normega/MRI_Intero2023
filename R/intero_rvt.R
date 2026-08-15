# intero_rvt.R
# Respiration volume per time (RVT) and breathing rate, sampled per functional
# volume.
#
# RVT is Birn et al. (2006): for each breath, the peak-to-trough excursion
# divided by the breath period, then interpolated onto the volume grid. It is
# the standard respiratory nuisance regressor and it captures something PetCO2
# does not -- a change in depth at constant rate moves RVT while barely moving
# end-tidal CO2, and vice versa.
#
# THE RRF IS THE CORRECT KERNEL HERE, unlike for PetCO2. Birn's respiration
# response function was derived as the impulse response FROM RVT TO BOLD, so
# RVT convolved with it is a prediction of the BOLD signal and the expected
# correlation is POSITIVE. The same kernel applied to PetCO2 predicts the wrong
# sign, because ventilating more raises RVT while lowering end-tidal CO2 (see
# intero_rrf.R). This module therefore uses rrf_birn() deliberately, and the
# sign of the result is a check on that reasoning rather than an assumption.
#
# UNITS. The belt is an uncalibrated strain gauge whose gain was set by hand
# each session, so RVT amplitude is NOT comparable across participants. Every
# series is z-scored within run before it is written. That is standard for RVT
# regressors and it is not optional here: raw RVT would make a group analysis a
# comparison of amplifier settings. Rate, in breaths per minute, IS comparable
# and is left in its own units.
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

RVT_MIN_BREATHS <- 5L    # below this a run cannot support an interpolated series

#' Per-volume RVT from a breath table.
#'
#' respkit's breath_table() already computes rvt = amplitude / duration per
#' cycle. Only kept breaths are used: a cycle flagged bad_amplitude has a
#' "peak" below one of its own troughs, and its RVT is meaningless.
#'
#' The value is placed at the MIDPOINT of each breath rather than its onset.
#' RVT summarises a whole cycle, so anchoring it to the cycle's start shifts
#' the series systematically early by half a breath -- about 2.3 s here, which
#' is more than one TR.
rvt_series <- function(bt, n_vol, tr = 2) {
  if (!nrow(bt) || !"rvt" %in% names(bt)) return(rep(NA_real_, n_vol))
  b <- bt[bt$keep %in% TRUE & is.finite(bt$rvt), , drop = FALSE]
  if (nrow(b) < RVT_MIN_BREATHS) return(rep(NA_real_, n_vol))
  t_mid <- (b$t_start + b$t_end) / 2
  t_vol <- (seq_len(n_vol) - 1) * tr
  stats::approx(t_mid, b$rvt, xout = t_vol, method = "linear", rule = 1)$y
}

#' Per-volume breathing rate (breaths per minute), same anchoring as RVT.
rate_series <- function(bt, n_vol, tr = 2) {
  if (!nrow(bt) || !"duration" %in% names(bt)) return(rep(NA_real_, n_vol))
  b <- bt[bt$keep %in% TRUE & is.finite(bt$duration) & bt$duration > 0, , drop = FALSE]
  if (nrow(b) < RVT_MIN_BREATHS) return(rep(NA_real_, n_vol))
  t_mid <- (b$t_start + b$t_end) / 2
  t_vol <- (seq_len(n_vol) - 1) * tr
  stats::approx(t_mid, 60 / b$duration, xout = t_vol, method = "linear", rule = 1)$y
}

#' z-score within run, ignoring NA, so amplitude is comparable across sessions.
z_run <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' RVT convolved with the RRF, on the fine grid then resampled to TR.
#'
#' Convolving the TR-sampled series directly would alias the kernel's 3.1 s
#' rise against a 2 s sample interval.
rvt_convolved <- function(bt, n_vol, tr = 2, dt = RRF_DT, lag_s = 0) {
  if (!nrow(bt) || !"rvt" %in% names(bt)) return(rep(NA_real_, n_vol))
  b <- bt[bt$keep %in% TRUE & is.finite(bt$rvt), , drop = FALSE]
  if (nrow(b) < RVT_MIN_BREATHS) return(rep(NA_real_, n_vol))
  t_mid <- (b$t_start + b$t_end) / 2
  n_f <- ceiling(n_vol * tr / dt) + 1L
  t_f <- (seq_len(n_f) - 1) * dt - lag_s
  xf <- stats::approx(t_mid, b$rvt, xout = t_f, method = "linear", rule = 1)$y
  yf <- convolve_causal(xf, kernel_grid(rrf_birn, dt))
  stats::approx((seq_len(n_f) - 1) * dt, yf, xout = (seq_len(n_vol) - 1) * tr, rule = 1)$y
}
