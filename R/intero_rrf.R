# intero_rrf.R
# Response-function convolution for the PetCO2 regressor.
#
# WHY CONVOLUTION IS NOT OPTIONAL
#
# PetCO2 is measured at the mouth; the BOLD effect it drives happens seconds
# later and is temporally smeared. Entering raw PetCO2 as a regressor asks the
# model to explain a delayed, dispersed response with an undelayed, sharp
# predictor. It underfits, and the vascular variance it fails to absorb stays
# in the residual where it can masquerade as neural -- which is exactly the
# variance the regressor was added to remove.
#
# TWO KERNELS ARE PROVIDED, AND THE DATA CHOOSE
#
# rrf_birn() is the respiration response function of Birn et al. (2008,
# NeuroImage 40:644-654). IT WAS DERIVED FOR RESPIRATION VOLUME PER TIME, not
# for PetCO2, and is used here because it is the standard kernel for
# respiratory nuisance regression -- not because it is known to be the right
# shape for end-tidal CO2.
#
# ITS SIGN IS WRONG FOR PetCO2, AND NOT BY ACCIDENT. RVT and end-tidal CO2 move
# in opposite directions: ventilating more blows off CO2, so high RVT means low
# PetCO2. The RRF's deep negative lobe encodes that inverse relationship, and
# its signed area is about -16. Convolving PetCO2 with it therefore predicts
# global BOLD to FALL when arterial CO2 rises, which is the opposite of
# established cerebrovascular physiology. Treat any positive correlation
# obtained with this kernel as a red flag, not a result.
#
# hrf_spm() is the canonical double-gamma. Some of the PetCO2 literature uses
# it instead, on the argument that a CO2-driven vasodilatory response is closer
# to a normal haemodynamic response than to the biphasic RRF.
#
# Rather than assert which is right, intero2023_rrf.R fits both against the
# global signal at a grid of lags and reports which explains more variance. The
# kernel is an empirical question and this dataset can answer it.
#
# LAG IS SEPARATE FROM THE KERNEL. Both kernels have their own peak latency,
# but there is also a bulk transit delay from lung to brain (of order 5-10 s)
# that varies with cardiac output and therefore with participant. Convolution
# handles the response shape; the fitted lag handles the transit.
#
# Author: RADlab

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

RRF_DT      <- 0.5    # s, grid the convolution is done on before resampling to TR
RRF_LEN_S   <- 45     # s, kernel support; RRF is essentially zero past ~40 s

# The range must be wide enough for the profile to TURN OVER inside it. At
# [0, 20] every kernel's mean profile rose monotonically to the boundary, which
# says only that the optimum is elsewhere. Negative lags are included as a
# control: a profile peaking below zero means BOLD leads CO2, which is
# impossible, and indicates misalignment rather than physiology.
RRF_LAG_MIN <- -10
RRF_LAG_MAX <- 40

#' Birn et al. (2008) respiration response function.
#'
#' RRF(t) = 0.6 t^2.1 exp(-t/1.6) - 0.0023 t^3.54 exp(-t/4.25)
rrf_birn <- function(t) {
  ifelse(t < 0, 0,
         0.6 * t^2.1 * exp(-t / 1.6) - 0.0023 * t^3.54 * exp(-t / 4.25))
}

#' SPM canonical haemodynamic response: two gamma densities, peak 6 s,
#' undershoot 16 s, undershoot 1/6 the amplitude.
hrf_spm <- function(t) {
  ifelse(t < 0, 0,
         stats::dgamma(t, shape = 6, scale = 1) -
           stats::dgamma(t, shape = 16, scale = 1) / 6)
}

#' Sample a kernel on the convolution grid, normalised by ABSOLUTE area.
#'
#' NOT by signed area. The Birn RRF is biphasic and its lobes do not cancel:
#' the positive lobe integrates to about +5.7 and the negative undershoot to
#' about -21.9, so its signed area is roughly -16. Dividing by that inverts the
#' kernel, which silently flips the sign of every regressor built from it and of
#' every correlation used to validate them. Absolute area is stable, cannot be
#' near zero for a non-trivial kernel, and preserves sign.
#'
#' Scale is therefore arbitrary, which costs nothing here: correlation is
#' invariant to positive rescaling, and a GLM beta absorbs it. Sign is what had
#' to be protected.
kernel_grid <- function(fn, dt = RRF_DT, len_s = RRF_LEN_S) {
  t <- seq(0, len_s, by = dt)
  k <- fn(t)
  a <- sum(abs(k)) * dt
  if (!is.finite(a) || a == 0) return(k)
  k / a * dt
}

#' Causal convolution with zero padding, NAs preserved.
#'
#' The input is DEMEANED first, which is what makes zero padding correct: the
#' run did not begin with no CO2 in the mask, it began at the participant's
#' baseline, so padding with zeros after demeaning pads with that baseline.
#' Padding the raw series with zeros would inject a spurious step at run onset
#' the size of the whole baseline.
#'
#' Interior NAs are linearly interpolated before convolving, because a single
#' missing breath would otherwise blank the following 45 s of regressor. The
#' output is re-masked to NA wherever the input was NA, so gaps stay visible.
convolve_causal <- function(x, k) {
  n <- length(x)
  na <- is.na(x)
  if (all(na)) return(rep(NA_real_, n))
  xi <- x
  if (any(na)) {
    idx <- which(!na)
    xi <- stats::approx(idx, x[idx], xout = seq_len(n), rule = 2)$y
  }
  xi <- xi - mean(xi, na.rm = TRUE)
  m <- length(k)
  y <- stats::filter(c(rep(0, m - 1), xi), k, method = "convolution", sides = 1)
  y <- as.numeric(y)[m:(m + n - 1)]
  y[na] <- NA_real_
  y
}

#' Build a TR-sampled PetCO2 regressor, optionally convolved and lagged.
#'
#' @param pk    breath table from petco2_detect()/petco2_flag()
#' @param n_vol volumes in the run
#' @param tr    repetition time
#' @param kernel one of "none", "rrf", "hrf"
#' @param lag_s transit delay applied to the PetCO2 series before convolution
#'
#' The series is built on a fine grid (RRF_DT) rather than at TR, because
#' convolving a 2 s-sampled signal with a kernel whose rise is ~2 s aliases the
#' rise. It is resampled to TR only at the end.
petco2_regressor <- function(pk, n_vol, tr = 2, kernel = c("none", "rrf", "hrf"),
                             lag_s = 0, dt = RRF_DT) {
  kernel <- match.arg(kernel)
  if (!nrow(pk) || !"ok" %in% names(pk)) return(rep(NA_real_, n_vol))
  dur <- n_vol * tr
  n_f <- ceiling(dur / dt) + 1L
  # lag by sampling the fine grid EARLIER in the physiological series: volume
  # at time t is driven by CO2 at t - lag.
  t_f <- (seq_len(n_f) - 1) * dt - lag_s
  good <- pk[pk$ok %in% TRUE, , drop = FALSE]
  if (nrow(good) < 3) return(rep(NA_real_, n_vol))
  xf <- stats::approx(good$t_s, good$petco2, xout = t_f, method = "linear",
                      rule = 1)$y
  yf <- if (kernel == "none") xf - mean(xf, na.rm = TRUE) else
    convolve_causal(xf, kernel_grid(if (kernel == "rrf") rrf_birn else hrf_spm, dt))
  t_vol <- (seq_len(n_vol) - 1) * tr
  stats::approx((seq_len(n_f) - 1) * dt, yf, xout = t_vol, rule = 1)$y
}

#' Band-pass both series to the frequency range a physiological response lives in.
#'
#' Detrending alone leaves the correlation dominated by very slow common
#' variance, which is what made the lag unidentifiable: at 0.01 Hz one cycle is
#' 100 s, so shifting by 40 s moves less than half a cycle and the profile stays
#' flat. Restricting to 0.01-0.10 Hz raised the median per-run correlation from
#' 0.240 to 0.288 and halved the proportion of runs whose best lag pinned to a
#' search boundary (25% -> 16%).
PETCO2_BAND <- c(0.01, 0.10)

bandpass_ts <- function(x, fs, band = PETCO2_BAND) {
  ok <- is.finite(x)
  if (sum(ok) < 50) return(rep(NA_real_, length(x)))
  xi <- stats::approx(which(ok), x[ok], xout = seq_along(x), rule = 2)$y
  xi <- xi - mean(xi)
  hi <- min(band[2], 0.95 * fs / 2)
  b <- signal::butter(2, c(band[1], hi) / (fs / 2), type = "pass")
  y <- as.numeric(signal::filtfilt(b, xi))
  y[!ok] <- NA_real_
  y
}

cor_band <- function(x, y, fs, band = PETCO2_BAND) {
  a <- bandpass_ts(x, fs, band); b <- bandpass_ts(y, fs, band)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 50) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok]))
}

#' Remove polynomial drift, so a correlation is not carried by shared trend.
#'
#' Both the global signal and PetCO2 drift over a run -- the scanner thermally,
#' the participant physiologically -- and a raw correlation between two
#' drifting series mostly measures that both drifted. Degree 3 over a 4-8 min
#' run is roughly a 1/128 Hz high pass, the conventional choice.
detrend_poly <- function(x, degree = 3) {
  n <- length(x); i <- seq_len(n)
  ok <- is.finite(x)
  if (sum(ok) < degree + 2) return(rep(NA_real_, n))
  fit <- stats::lm(x[ok] ~ stats::poly(i[ok], degree))
  out <- rep(NA_real_, n)
  out[ok] <- stats::residuals(fit)
  out
}

#' Correlation of a regressor with an observed series, both detrended.
cor_detrended <- function(x, y, degree = 3) {
  a <- detrend_poly(x, degree); b <- detrend_poly(y, degree)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 10) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok]))
}

#' Search transit lag, returning the whole profile rather than just the peak.
#'
#' The profile matters: a real vascular response gives a smooth curve with one
#' broad maximum, while noise gives a ragged one whose argmax is not meaningful.
lag_profile <- function(pk, n_vol, y, tr = 2, kernel = "rrf",
                        lags = seq(RRF_LAG_MIN, RRF_LAG_MAX, by = 2),
                        band = PETCO2_BAND) {
  r <- vapply(lags, function(L) {
    x <- petco2_regressor(pk, n_vol, tr = tr, kernel = kernel, lag_s = L)
    if (is.null(band)) cor_detrended(x, y) else cor_band(x, y, fs = 1 / tr, band)
  }, numeric(1))
  tibble(lag_s = lags, r = r)
}

# ---------------------------------------------------------------------------
# THE LAG IS NOT IDENTIFIABLE FROM THESE DATA, SO IT IS NOT FITTED
#
# The band-passed group profile is smooth and correctly signed -- strongly
# negative at negative lags, crossing zero at about 9 s, then rising -- but it
# reaches a plateau near 18 s and stays flat out to 40 s instead of turning
# over. The physiology is therefore bounded BELOW (the response cannot be
# faster than ~9 s) and unbounded above within the range the data can resolve.
#
# Reporting a fitted median lag from that profile would be inventing precision:
# per-run argmaxes scatter across the whole range and 16% still pin to a
# boundary. Instead the regressors are written as a BASIS of fixed lags, and
# the GLM is left to weight them. That is the standard treatment for a
# physiological regressor whose delay is uncertain, and it degrades gracefully
# -- if the true lag is 15 s, the 10 s and 20 s columns straddle it.
PETCO2_LAG_BASIS <- c(0, 10, 20, 30)
