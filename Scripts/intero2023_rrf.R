# intero2023_rrf.R
# Fits the PetCO2 -> BOLD response empirically against fMRIPrep's global signal,
# then writes convolved, lag-corrected PetCO2 regressors.
#
# THE SANITY CHECK THIS PERFORMS
#
# PetCO2 must predict GLOBAL BOLD positively. That is settled physiology:
# rising arterial CO2 dilates cerebral vessels and raises the BOLD signal
# almost everywhere. If the regressor does not do that in these data, the
# regressor is wrong -- and any later finding about it removing or sparing
# effects of interest would be uninterpretable. This runs before any of that.
#
# Global signal comes from Data/Ready (fMRIPrep 24.1.1 confounds), so no images
# are read.
#
# [2026-11-24: NO LONGER EXCLUDED. The shift below is a property of the old
#  Data/Ready snapshot only; the blocks and the BIDS tree of record map these
#  runs by scanner slot. Original note kept for the record:]
# sub-082 and sub-104 ARE EXCLUDED. Their run labels are shifted: the file
# named task-breath1 holds 172 volumes where the breathing slot expects 242,
# and 172 is the authenticity length. Every other run in the tree agrees with
# inventory.tsv exactly (242 of 244), so the mapping is sound apart from these
# two, and guessing which physiology belongs to which image is precisely the
# error that reaches results unnoticed.
#
# Produces:
#   Repo/Data/petco2_rrf_fit.csv               per run: kernel, best lag, r
#   Repo/Data/regressors/<sub>_<entity>_petco2.tsv   GLM-ready regressors
#   Repo/Figures/petco2_rrf.png                kernels and lag profiles
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
  library(signal)   # masks filter(); dplyr:: is qualified below where it matters
})

# Paths default to the Windows shared drive; INTERO_ROOT / RESPKIT_R / INTERO_FMRIPREP
# override them (2026-11-24: lets the same script run against the Linux tree).
dir_root  <- Sys.getenv("INTERO_ROOT", "I:/Shared drives/MRI_Intero2023")
dir_out   <- file.path(dir_root, "Repo", "Data")
dir_R     <- file.path(dir_root, "Repo", "R")
dir_fig   <- file.path(dir_root, "Repo", "Figures")
dir_reg   <- file.path(dir_out, "regressors")
dir_cache <- file.path(dir_root, "Data", "Physio", "cache")
dir_ready <- Sys.getenv("INTERO_FMRIPREP", file.path(dir_root, "Data", "Ready"))
TR <- 2

rk <- Sys.getenv("RESPKIT_R", "I:/Shared drives/TechnicalReference/respkit/R")
for (f in list.files(rk, pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in c("intero_ids.R", "intero_qc.R", "intero_petco2.R", "intero_rrf.R"))
  source(file.path(dir_R, f))
dir.create(dir_reg, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))
cat(sprintf("RRF peak at %.1f s, HRF peak at %.1f s\n",
            which.max(rrf_birn(seq(0, 45, 0.1))) * 0.1 - 0.1,
            which.max(hrf_spm(seq(0, 45, 0.1))) * 0.1 - 0.1))

# EXCLUDE of sub-082/sub-104 removed 2026-11-24 -- see the note added to the header.
EXCLUDE <- character(0)
TASK_MAP <- c(loc1 = "task-localizer_run-1", loc2 = "task-localizer_run-2",
              breath1 = "task-breathing_run-1", breath2 = "task-breathing_run-2",
              auth1 = "task-authenticity")

conf <- tibble(path = list.files(dir_ready, pattern = "desc-confounds_timeseries[.]tsv$",
                                 recursive = TRUE, full.names = TRUE)) %>%
  mutate(bn = basename(path),
         subject = str_match(bn, "^(sub-[0-9]+)_")[, 2],
         bids_entity = confounds_entity(bn, TASK_MAP)) %>%   # old or BIDS naming
  dplyr::filter(!is.na(bids_entity), !subject %in% EXCLUDE)

pcr <- read_csv(file.path(dir_out, "petco2_runs.csv"), show_col_types = FALSE,
                progress = FALSE) %>% dplyr::filter(capnogram_valid)
cat(sprintf("valid capnogram runs: %d ; confound files usable: %d\n",
            nrow(pcr), nrow(conf)))

fits <- list(); profiles <- list(); nreg <- 0L
cache <- new.env()

for (i in seq_len(nrow(pcr))) {
  sb <- pcr$subject[i]; ent <- pcr$bids_entity[i]
  cp <- conf %>% dplyr::filter(subject == sb, bids_entity == ent)
  if (!nrow(cp)) next

  key <- paste0(sb, "_petco2")
  if (!exists(key, envir = cache)) {
    f <- file.path(dir_cache, paste0(sb, "_petco2.rds"))
    if (!file.exists(f)) next
    assign(key, readRDS(f), envir = cache)
  }
  r <- get(key, envir = cache)
  if (is.null(r$series[[ent]])) next
  pk <- r$series[[ent]]$breaths
  if (!nrow(pk) || !"ok" %in% names(pk)) next

  d <- suppressWarnings(read_tsv(cp$path[1], show_col_types = FALSE, progress = FALSE,
                                 col_types = cols(.default = col_character())))
  g <- suppressWarnings(as.numeric(d$global_signal))
  n_vol <- length(g)
  fd <- suppressWarnings(as.numeric(d$framewise_displacement))
  if (n_vol < 60) next

  row <- tibble(subject = sb, bids_entity = ent, n_vol = n_vol,
                task = str_match(ent, "task-([a-z]+)")[, 2],
                fd_med = median(fd, na.rm = TRUE))

  for (kn in c("none", "rrf", "hrf")) {
    lp <- lag_profile(pk, n_vol, g, tr = TR, kernel = kn)
    profiles[[length(profiles) + 1]] <- lp %>%
      mutate(subject = sb, bids_entity = ent, kernel = kn)
    if (all(is.na(lp$r))) { row[[paste0("r_", kn)]] <- NA_real_
                            row[[paste0("lag_", kn)]] <- NA_real_; next }
    b <- which.max(lp$r)      # max, not max|.|: the sign is a prediction
    row[[paste0("r_", kn)]]   <- lp$r[b]
    row[[paste0("lag_", kn)]] <- lp$lag_s[b]
    row[[paste0("r0_", kn)]]  <- lp$r[lp$lag_s == 0]
  }
  fits[[length(fits) + 1]] <- row
}

ft <- bind_rows(fits); pf <- bind_rows(profiles)
write_csv(ft, file.path(dir_out, "petco2_rrf_fit.csv"))

cat(sprintf("\nruns fitted: %d  participants: %d\n", nrow(ft), n_distinct(ft$subject)))

cat("\n== SANITY CHECK: does PetCO2 predict global BOLD positively? ==\n")
for (kn in c("none", "rrf", "hrf")) {
  rr <- ft[[paste0("r_", kn)]]
  cat(sprintf("  %-4s  median r = %+.3f   IQR [%+.3f, %+.3f]   positive in %.0f%% of runs\n",
              kn, median(rr, na.rm = TRUE),
              quantile(rr, .25, na.rm = TRUE), quantile(rr, .75, na.rm = TRUE),
              100 * mean(rr > 0, na.rm = TRUE)))
}
cat("\n-- variance explained at the fitted lag (r^2, median) --\n")
for (kn in c("none", "rrf", "hrf"))
  cat(sprintf("  %-4s  r2 = %.3f\n", kn, median(ft[[paste0("r_", kn)]]^2, na.rm = TRUE)))

cat("\n-- does convolution beat raw PetCO2, within run? --\n")
cmp <- ft %>% dplyr::filter(is.finite(r_none), is.finite(r_rrf), is.finite(r_hrf))
cat(sprintf("  rrf > none in %.0f%% of %d runs (median gain in r: %+.3f)\n",
            100 * mean(cmp$r_rrf > cmp$r_none), nrow(cmp),
            median(cmp$r_rrf - cmp$r_none)))
cat(sprintf("  hrf > none in %.0f%% of %d runs (median gain in r: %+.3f)\n",
            100 * mean(cmp$r_hrf > cmp$r_none), nrow(cmp),
            median(cmp$r_hrf - cmp$r_none)))
cat(sprintf("  hrf > rrf  in %.0f%% of runs (median difference: %+.3f)\n",
            100 * mean(cmp$r_hrf > cmp$r_rrf), median(cmp$r_hrf - cmp$r_rrf)))

cat("\n-- fitted transit lag (s) --\n")
for (kn in c("rrf", "hrf"))
  print(c(kernel = kn, round(quantile(ft[[paste0("lag_", kn)]], c(.25,.5,.75), na.rm = TRUE), 1)))

# A profile whose argmax sits on either boundary is not an estimate. Report how
# often that happens rather than quoting the median as if it were one.
cat("\n-- fraction of runs whose best lag is AT a search boundary --\n")
for (kn in c("none", "rrf", "hrf")) {
  L <- ft[[paste0("lag_", kn)]]
  cat(sprintf("  %-4s  at %+ds: %.0f%%   at %+ds: %.0f%%\n", kn,
              RRF_LAG_MIN, 100 * mean(L == RRF_LAG_MIN, na.rm = TRUE),
              RRF_LAG_MAX, 100 * mean(L == RRF_LAG_MAX, na.rm = TRUE)))
}
cat("\n-- group-mean profile: does it turn over inside the range? --\n")
gp0 <- pf %>% group_by(kernel, lag_s) %>%
  summarise(r = mean(r, na.rm = TRUE), .groups = "drop")
for (kn in c("none", "rrf", "hrf")) {
  g <- gp0 %>% dplyr::filter(kernel == kn) %>% arrange(lag_s)
  b <- which.max(g$r)
  cat(sprintf("  %-4s  peak mean r = %+.3f at lag %+.0fs %s\n", kn, g$r[b], g$lag_s[b],
              if (b == 1 || b == nrow(g)) "<-- AT BOUNDARY, not an estimate" else "(interior)"))
}

# ---------------------------------------------------------------------------
# THE NUMBERS ABOVE ARE MAXIMA OVER 26 LAGS AND ARE THEREFORE OPTIMISTIC.
# Taking the best of 26 correlated candidates inflates r substantially: the
# lag-selected median of +0.288 becomes +0.082 at a fixed 20 s lag. Quoting the
# selected value as an effect size would overstate how much of the global
# signal PetCO2 actually accounts for, so the unselected values are computed
# here and are what should be reported.
# ---------------------------------------------------------------------------
cat("\n== FIXED-LAG correlations (no selection) and the lag basis ==\n")
fx <- list()
for (i in seq_len(nrow(ft))) {
  sb <- ft$subject[i]; ent <- ft$bids_entity[i]; n_vol <- ft$n_vol[i]
  cp <- conf %>% dplyr::filter(subject == sb, bids_entity == ent)
  d <- suppressWarnings(read_tsv(cp$path[1], show_col_types = FALSE, progress = FALSE,
                                 col_types = cols(.default = col_character())))
  g <- suppressWarnings(as.numeric(d$global_signal))
  pk <- get(paste0(sb, "_petco2"), envir = cache)$series[[ent]]$breaths
  rr <- vapply(PETCO2_LAG_BASIS, function(L)
    cor_band(petco2_regressor(pk, n_vol, TR, "hrf", L), g, 1 / TR), numeric(1))
  X <- do.call(cbind, lapply(PETCO2_LAG_BASIS, function(L)
    bandpass_ts(petco2_regressor(pk, n_vol, TR, "hrf", L), 1 / TR)))
  y <- bandpass_ts(g, 1 / TR)
  ok <- is.finite(y) & apply(is.finite(X), 1, all)
  R2 <- if (sum(ok) > 30) summary(stats::lm(y[ok] ~ X[ok, ]))$r.squared else NA_real_
  fx[[length(fx) + 1]] <- tibble(subject = sb, bids_entity = ent, task = ft$task[i],
    lag = PETCO2_LAG_BASIS, r = rr, R2_basis = R2)
}
fxd <- bind_rows(fx)
for (L in PETCO2_LAG_BASIS) {
  v <- fxd$r[fxd$lag == L]
  cat(sprintf("  lag %2ds  median r = %+.3f   positive in %3.0f%% of runs\n",
              L, median(v, na.rm = TRUE), 100 * mean(v > 0, na.rm = TRUE)))
}
b1 <- fxd %>% distinct(subject, bids_entity, task, R2_basis)
cat(sprintf("\n  %d-lag HRF basis: median R2 = %.3f  IQR [%.3f, %.3f]\n",
            length(PETCO2_LAG_BASIS), median(b1$R2_basis, na.rm = TRUE),
            quantile(b1$R2_basis, .25, na.rm = TRUE),
            quantile(b1$R2_basis, .75, na.rm = TRUE)))
cat("\n  basis R2 and fixed 20s r by task:\n")
print(as.data.frame(fxd %>% dplyr::filter(lag == 20) %>% group_by(task) %>%
  summarise(n = n(), r_lag20 = round(median(r, na.rm = TRUE), 3),
            R2_basis = round(median(R2_basis, na.rm = TRUE), 3), .groups = "drop")))
write_csv(fxd, file.path(dir_out, "petco2_fixedlag_fit.csv"))

cat("\n-- by task (lag-selected, optimistic) --\n")
print(as.data.frame(ft %>% group_by(task) %>% summarise(n = n(),
  r_none = round(median(r_none, na.rm = TRUE), 3),
  r_rrf  = round(median(r_rrf,  na.rm = TRUE), 3),
  r_hrf  = round(median(r_hrf,  na.rm = TRUE), 3),
  lag_rrf = round(median(lag_rrf, na.rm = TRUE), 1), .groups = "drop")))

# ---- group lag profile, the shape test -------------------------------------
gp <- pf %>% group_by(kernel, lag_s) %>%
  summarise(r = mean(r, na.rm = TRUE), .groups = "drop")

png(file.path(dir_fig, "petco2_rrf.png"), width = 1400, height = 560, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
tt <- seq(0, 45, 0.1)
plot(tt, rrf_birn(tt) / max(rrf_birn(tt)), type = "l", lwd = 2, col = "steelblue",
     xlab = "time (s)", ylab = "normalised", main = "kernels")
lines(tt, hrf_spm(tt) / max(hrf_spm(tt)), lwd = 2, col = "firebrick")
abline(h = 0, col = "grey70"); legend("topright", c("RRF (Birn 2008)", "HRF (SPM)"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
with(gp %>% dplyr::filter(kernel == "rrf"), plot(lag_s, r, type = "b", pch = 16,
     col = "steelblue", ylim = range(gp$r, na.rm = TRUE), xlab = "transit lag (s)",
     ylab = "mean r with global signal", main = "lag profile"))
with(gp %>% dplyr::filter(kernel == "hrf"), lines(lag_s, r, type = "b", pch = 16, col = "firebrick"))
with(gp %>% dplyr::filter(kernel == "none"), lines(lag_s, r, type = "b", pch = 16, col = "grey40"))
abline(h = 0, col = "grey70")
legend("topright", c("rrf", "hrf", "none"), col = c("steelblue", "firebrick", "grey40"),
       lwd = 2, bty = "n")
hist(ft$r_rrf, breaks = 30, col = "grey80", border = "white",
     main = "r at fitted lag (RRF)", xlab = "r"); abline(v = 0, col = "red", lwd = 2)
par(op); dev.off()

# ---- write the regressors --------------------------------------------------
# A BASIS of fixed lags, not a fitted one. See PETCO2_LAG_BASIS in intero_rrf.R:
# the lag is bounded below at ~9 s but not identifiable above, so any single
# fitted value would be false precision. Fitting per run against the global
# signal and then using that same signal downstream would additionally let the
# regressor absorb whatever it happened to fit, including neural variance.
#
# HRF, not RRF. The RRF is the wrong-signed kernel for PetCO2 (it was derived
# for RVT, which moves opposite to end-tidal CO2); it is written only as
# petco2_rrf_diag for diagnostics and should not enter a model.
cat(sprintf("\nwriting regressors: HRF at lags %s s\n",
            paste(PETCO2_LAG_BASIS, collapse = ", ")))

for (i in seq_len(nrow(ft))) {
  sb <- ft$subject[i]; ent <- ft$bids_entity[i]; n_vol <- ft$n_vol[i]
  r <- get(paste0(sb, "_petco2"), envir = cache)
  pk <- r$series[[ent]]$breaths
  out <- tibble(petco2_raw = petco2_regressor(pk, n_vol, TR, "none", 0))
  for (L in PETCO2_LAG_BASIS)
    out[[sprintf("petco2_hrf_lag%02d", L)]] <-
      petco2_regressor(pk, n_vol, TR, "hrf", L)
  out$petco2_rrf_diag <- petco2_regressor(pk, n_vol, TR, "rrf", 20)
  write_tsv(out %>% mutate(across(everything(), ~round(.x, 6))),
            file.path(dir_reg, sprintf("%s_%s_petco2.tsv", sb, ent)), na = "n/a")
  nreg <- nreg + 1L
}
cat(sprintf("regressor files written: %d (columns: raw, %d HRF lags, rrf diagnostic)\n",
            nreg, length(PETCO2_LAG_BASIS)))

qc <- pmap_dfr(list(ft$subject, ft$bids_entity, ft$r_rrf), function(s, e, rr)
  qc_entry(s, "physio", "petco2_predicts_global_bold",
           if (!is.finite(rr)) "fail" else if (rr > 0.3) "good"
           else if (rr > 0.1) "info" else "degraded",
           entity = e, value = rr,
           note = sprintf("r with global signal, RRF-convolved at fitted lag"),
           source = "fMRIPrep confounds global_signal"))
qc_write(qc_merge(qc_read(file.path(dir_out, "qc_ledger.csv")), qc),
         file.path(dir_out, "qc_ledger.csv"))

cat(sprintf("\nRun ended: %s\n", Sys.time()))
