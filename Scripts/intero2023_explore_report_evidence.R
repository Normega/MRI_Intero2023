# intero2023_explore_report_evidence.R -- EXPLORATORY E14 (deviations.md, "After the update", 2026-11-28).
# Behaviour only, no imaging. Does the post-epoch report (Faster / No Change / Slower) follow the
# instructed pacer change (which is VISIBLE: the circle itself speeds up or slows down), the
# participant's own measured breathing change, or both?
#
#   Rscript Scripts/intero2023_explore_report_evidence.R            # the logged analysis
#   Rscript Scripts/intero2023_explore_report_evidence.R --dryrun   # label-free code check: responses
#       permuted within participant (seed 1), outputs and figure to $E14_DRYRUN_DIR (default tempdir())
#
# Inputs : Data/breathing_trials.csv, Data/resp_epochs.csv (merged on subject x bids_entity x
#          trial_type x epoch_index), and the per-breath belt tables in
#          ../sourcedata/physio/cache/<sub>_resp.rds (not in the repo; used only for the free-epoch
#          transition ratio and the within-epoch trend standard errors).
# Outputs: Analysis/results/EXPLORATORY_E14_{counts,primary,measurement_error,nochange,
#          dissociation,free,confidence,sensitivity}.csv ; Figures/EXPLORATORY_E14_report_evidence.png
#
# Conventions follow Scripts/intero2023_behaviour.R: the report is a 3-category identification,
# ordered Faster < No Change < Slower, ordinal probit (ordinal::clmm) with a participant intercept.
# Ratios are breath PERIOD ratios, so log ratio > 0 = period lengthened = slower.
#   x = log(actual_ratio): measured within-epoch log-linear trend of contained breaths, end vs start
#   z = log(instructed_ratio): the pacer's trend fitted over the SAME window (R/intero_resp.R)
#   sal = +1/2 step ("high"), -1/2 ramp ("low"), 0 No Change

# Set Up ---------
## Load libraries ---------
packages <- c("tidyverse", "ordinal", "lme4", "lmerTest")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE, width = 160)
for (thispack in packages) {
  suppressPackageStartupMessages(library(thispack, character.only = TRUE, quietly = TRUE,
                                         verbose = FALSE))
}

args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)])))
repo <- dirname(here)
dryrun <- "--dryrun" %in% commandArgs(trailingOnly = TRUE)
cache <- file.path(dirname(repo), "sourcedata", "physio", "cache")
rd <- if (dryrun) Sys.getenv("E14_DRYRUN_DIR", tempdir()) else file.path(repo, "Analysis", "results")
dir.create(rd, showWarnings = FALSE, recursive = TRUE)
fig <- if (dryrun) file.path(rd, "EXPLORATORY_E14_report_evidence_dryrun.png") else
  file.path(repo, "Figures", "EXPLORATORY_E14_report_evidence.png")

SEED <- 20261128
NBOOT_ME <- 1000      # participant bootstrap for the measurement-error correction
NBOOT_P  <- 2000      # participant bootstrap for proportions
if (dryrun) { NBOOT_ME <- 50; NBOOT_P <- 200 }   # code check only
PLAN_SIGMA <- 0.40    # plan noise floor, SD of the per-epoch log period ratio
STRICT_X <- 0.20      # stricter dissociation rule: |x| >= .20
STALE_CACHE <- "sub-070"   # breath cache is stale for this participant (Data/STALE_TABLES.md)
LV <- c("Faster", "No Change", "Slower")

wald <- function(fit, terms) {
  cs <- coef(summary(fit))
  tibble::tibble(term = terms, est = cs[terms, 1], se = cs[terms, 2],
                 lo = cs[terms, 1] - 1.96 * cs[terms, 2], hi = cs[terms, 1] + 1.96 * cs[terms, 2],
                 z = cs[terms, 1] / cs[terms, 2], p = 2 * pnorm(-abs(cs[terms, 1] / cs[terms, 2])))
}
fit_clmm <- function(f, d) ordinal::clmm(f, data = d, link = "probit", Hess = TRUE)

# Data ---------
key <- c("subject", "bids_entity", "trial_type", "epoch_index")
b <- readr::read_csv(file.path(repo, "Data", "breathing_trials.csv"), na = c("NA", "n/a", ""))
r <- readr::read_csv(file.path(repo, "Data", "resp_epochs.csv"), na = c("NA", "n/a", ""))
stopifnot(!anyDuplicated(b[, key]), !anyDuplicated(r[, key]))

if (dryrun) {
  set.seed(1)
  b <- b |> dplyr::group_by(subject, trial_type) |>
    dplyr::mutate(response = sample(response), confidence = sample(confidence)) |> dplyr::ungroup()
  cat("DRY RUN: responses and confidence permuted within participant x epoch type\n")
}

# previous paced report in the same run (signed -1/0/+1; 0 and first = 1 for the first paced epoch)
prev <- b |> dplyr::filter(trial_type == "paced_breathing") |>
  dplyr::arrange(subject, bids_entity, epoch_index) |>
  dplyr::group_by(subject, bids_entity) |>
  dplyr::mutate(prev_resp = dplyr::lag(c(Faster = -1, `No Change` = 0, Slower = 1)[response]),
                first = as.integer(is.na(prev_resp)),
                prev_resp = dplyr::coalesce(unname(prev_resp), 0)) |>
  dplyr::ungroup() |> dplyr::select(subject, bids_entity, epoch_index, prev_resp, first)

m <- dplyr::inner_join(
  r |> dplyr::select(-response, -sdt_class),
  b |> dplyr::select(dplyr::all_of(key), change_amount, confidence, response),
  by = key)
stopifnot(nrow(m) == nrow(r))

# Primary sample: paced, belt usable (measured = >= 3 breaths, coverage >= .60, belt good/degraded),
# finite measured and instructed ratio, a response. No further exclusion.
p <- m |>
  dplyr::filter(trial_type == "paced_breathing", measured %in% TRUE, is.finite(actual_ratio),
                is.finite(instructed_ratio), !is.na(response)) |>
  dplyr::left_join(prev, by = c("subject", "bids_entity", "epoch_index")) |>
  dplyr::mutate(x = log(actual_ratio), z = log(instructed_ratio),
                sal = dplyr::case_when(change_direction == "No Change" ~ 0,
                                       salience == "high" ~ 0.5, salience == "low" ~ -0.5),
                resp = factor(response, levels = LV, ordered = TRUE),
                subject = factor(subject))
stopifnot(!anyNA(p$sal), !anyNA(p$prev_resp))

# Breath cache: within-epoch trend SE (No Change epochs) and free-epoch transition ratio ---------
rows <- list()
for (s in setdiff(unique(m$subject), STALE_CACHE)) {
  f <- file.path(cache, paste0(s, "_resp.rds"))
  if (!file.exists(f)) next
  bx <- readRDS(f)$breaths
  ev <- m |> dplyr::filter(subject == s)
  for (k in seq_len(nrow(ev))) {
    bt <- bx[[ev$bids_entity[k]]]
    if (is.null(bt)) next
    w <- bt$t_start >= ev$onset[k] & bt$t_end <= ev$onset[k] + ev$duration[k] & bt$keep %in% TRUE
    bb <- bt[w, , drop = FALSE]
    se <- NA_real_; xr <- NA_real_
    if (nrow(bb) >= 3) {
      span <- diff(range(bb$t_start))
      ft <- lm(log(bb$duration) ~ bb$t_start)
      xr <- unname(coef(ft)[2]) * span
      se <- suppressWarnings(summary(ft)$coefficients[2, 2]) * span   # 3 breaths: 1 df
    }
    rows[[length(rows) + 1]] <- tibble::tibble(
      subject = s, bids_entity = ev$bids_entity[k], trial_type = ev$trial_type[k],
      epoch_index = ev$epoch_index[k], n_cont = nrow(bb), x_rebuilt = xr, se_x = se,
      log_med = if (nrow(bb)) log(median(bb$duration)) else NA_real_,
      log_last2 = if (nrow(bb) >= 2) mean(log(tail(bb$duration, 2))) else NA_real_)
  }
}
bc <- dplyr::bind_rows(rows) |>
  dplyr::left_join(m |> dplyr::select(dplyr::all_of(key), actual_ratio), by = key)
chk <- with(bc, abs(x_rebuilt - log(actual_ratio)))
cat(sprintf("cache rebuild of actual_ratio: %d of %d match to 1e-6\n",
            sum(chk < 1e-6, na.rm = TRUE), sum(is.finite(chk))))

# Counts ---------
nrun <- function(d) dplyr::n_distinct(paste(d$subject, d$bids_entity))
n_sub_p <- dplyr::n_distinct(p$subject); n_ep_p <- nrow(p)
counts <- tibble::tibble(set = "paced_primary", epochs = nrow(p),
                         participants = dplyr::n_distinct(p$subject), runs = nrun(p))
cat("\nprimary sample:\n"); print(table(p$change_direction, p$salience, useNA = "ifany"))

# 1. Primary model ---------
f1 <- fit_clmm(resp ~ z + x + sal + (1 | subject), p)
sd_zx <- c(sd(p$z), sd(p$x), NA)
prim <- wald(f1, c("z", "x", "sal")) |>
  dplyr::mutate(sd_pred = sd_zx,
                est_per_sd = est * sd_pred, lo_per_sd = lo * sd_pred, hi_per_sd = hi * sd_pred,
                p_holm = c(p.adjust(p[1:2], "holm"), NA))
r_xz <- cor(p$x, p$z)
pc <- p |> dplyr::group_by(subject) |>
  dplyr::mutate(xc = x - mean(x), zc = z - mean(z)) |> dplyr::ungroup()
vif <- function(v, others) 1 / (1 - summary(lm(reformulate(others, v), p))$r.squared)
diag_tab <- tibble::tibble(
  term = c("r(x,z)", "r(x,z) within participant", "VIF z", "VIF x", "VIF sal",
           "participant intercept SD", "threshold Faster|NoChange", "threshold NoChange|Slower"),
  est = c(r_xz, cor(pc$xc, pc$zc), vif("z", c("x", "sal")), vif("x", c("z", "sal")),
          vif("sal", c("x", "z")), sqrt(as.numeric(ordinal::VarCorr(f1)$subject)), f1$Theta))
cat("\nPRIMARY: resp ~ z + x + sal + (1 | subject), probit\n")
print(as.data.frame(prim), digits = 3); print(as.data.frame(diag_tab), digits = 3)
readr::write_csv(dplyr::bind_rows(prim |> dplyr::mutate(model = "primary", .before = 1),
                                  diag_tab |> dplyr::mutate(model = "diagnostic", .before = 1)),
                 file.path(rd, "EXPLORATORY_E14_primary.csv"))

# 2. Measurement error ---------
# Regression calibration (linear approximation on the latent probit scale). With X = X* + u and
# lambda = 1 - sigma_u^2 / Var(X | z, sal):
#   b_x = beta_x * lambda ;  b_z = beta_z + beta_x * gamma * (1 - lambda)
# gamma = coefficient of z in lm(x ~ z + sal). lambda_star = reliability at which beta_z = 0.
me_core <- function(d) {
  ff <- fit_clmm(resp ~ z + x + sal + (1 | subject), d)
  lx <- lm(x ~ z + sal, d)
  list(bz = unname(coef(ff)["z"]), bx = unname(coef(ff)["x"]),
       gam = unname(coef(lx)["z"]), vres = summary(lx)$sigma^2)
}
correct <- function(cc, sig) {
  lam <- 1 - sig^2 / cc$vres
  if (!is.finite(lam) || lam <= 0) return(c(lambda = lam, bx_c = NA, bz_c = NA))
  c(lambda = lam, bx_c = cc$bx / lam, bz_c = cc$bz - cc$bx * cc$gam * (1 - lam) / lam)
}
lam_star <- function(cc) {
  v <- cc$bx * cc$gam / (cc$bz + cc$bx * cc$gam)
  if (cc$bx > 0 && cc$bz > 0) v else NA_real_
}
cc0 <- me_core(p)
sig_nc <- sd(p$x[p$change_direction == "No Change"])
se_nc <- bc |> dplyr::filter(trial_type == "paced_breathing", is.finite(se_x)) |>
  dplyr::semi_join(p |> dplyr::filter(change_direction == "No Change") |>
                     dplyr::mutate(subject = as.character(subject)),
                   by = c("subject", "bids_entity", "epoch_index"))
sig_se <- sqrt(mean(se_nc$se_x^2))
sigmas <- c(illustr_0.25 = 0.25, illustr_0.30 = 0.30, illustr_0.35 = 0.35,
            nochange_sd = sig_nc, plan_0.40 = PLAN_SIGMA, trend_se_nochange = sig_se)

set.seed(SEED)
subs <- levels(p$subject)
boot <- purrr::map_dfr(seq_len(NBOOT_ME), function(i) {
  ss <- sample(subs, length(subs), replace = TRUE)
  d <- purrr::imap_dfr(ss, function(s, j) p |> dplyr::filter(subject == s) |>
                         dplyr::mutate(subject = paste0("b", j)))
  d$subject <- factor(d$subject)
  cc <- tryCatch(me_core(d), error = function(e) NULL)
  if (is.null(cc)) return(tibble::tibble())
  tibble::as_tibble(c(list(i = i, lam_star = lam_star(cc)),
    unlist(lapply(names(sigmas), function(nm) {
      v <- correct(cc, sigmas[[nm]]); setNames(as.list(v), paste(nm, names(v), sep = ":")) }),
      recursive = FALSE)))
})
qci <- function(v) if (sum(is.finite(v)) < 0.9 * NBOOT_ME) c(NA, NA) else
  unname(quantile(v, c(.025, .975), na.rm = TRUE))
me <- purrr::map_dfr(names(sigmas), function(nm) {
  v <- correct(cc0, sigmas[[nm]])
  tibble::tibble(sigma_label = nm, sigma_u = sigmas[[nm]], var_x_given_z = cc0$vres,
                 gamma = cc0$gam, lambda = v[["lambda"]],
                 bx_obs = cc0$bx, bz_obs = cc0$bz, bx_corrected = v[["bx_c"]],
                 bx_c_lo = qci(boot[[paste0(nm, ":bx_c")]])[1],
                 bx_c_hi = qci(boot[[paste0(nm, ":bx_c")]])[2],
                 bz_corrected = v[["bz_c"]],
                 bz_c_lo = qci(boot[[paste0(nm, ":bz_c")]])[1],
                 bz_c_hi = qci(boot[[paste0(nm, ":bz_c")]])[2])
})
ls0 <- lam_star(cc0)
me <- dplyr::bind_rows(me, tibble::tibble(
  sigma_label = "tipping_point", lambda = ls0,
  sigma_u = if (is.finite(ls0)) sqrt(cc0$vres * (1 - ls0)) else NA_real_,
  var_x_given_z = cc0$vres, gamma = cc0$gam, bx_obs = cc0$bx, bz_obs = cc0$bz,
  lambda_lo = qci(boot$lam_star)[1], lambda_hi = qci(boot$lam_star)[2],
  n_boot_finite = sum(is.finite(boot$lam_star))))
cat(sprintf("\nMEASUREMENT ERROR (bootstrap fits %d of %d; lambda_star finite in %d)\n",
            nrow(boot), NBOOT_ME, sum(is.finite(boot$lam_star))))
print(as.data.frame(me), digits = 3)
readr::write_csv(me, file.path(rd, "EXPLORATORY_E14_measurement_error.csv"))

# Participant-cluster bootstrap for a proportion and a difference ---------
pboot <- function(d, fn) {
  set.seed(SEED)
  ss <- unique(as.character(d$subject))
  v <- replicate(NBOOT_P, {
    k <- sample(ss, length(ss), replace = TRUE)
    fn(dplyr::bind_rows(lapply(k, function(s) d[as.character(d$subject) == s, ])))
  })
  unname(quantile(v, c(.025, .975), na.rm = TRUE))
}

# 3a. No Change paced epochs ---------
pnc <- p |> dplyr::filter(change_direction == "No Change") |> dplyr::mutate(subject = droplevels(subject))
# Pre-set rule: clmm is the estimate unless it fails or is ill-conditioned (Hessian condition
# number > 1e6 or SE(x) < 1e-3), in which case clm (no participant term) is the estimate.
fnc <- tryCatch(fit_clmm(resp ~ x + (1 | subject), pnc), error = function(e) NULL)
fnc0 <- ordinal::clm(resp ~ x, data = pnc, link = "probit")
nc_ok <- !is.null(fnc) && is.finite(coef(summary(fnc))["x", 2]) &&
  coef(summary(fnc))["x", 2] >= 1e-3 && summary(fnc)$info$cond.H < 1e6
ncw <- dplyr::bind_rows(
  if (!is.null(fnc)) wald(fnc, "x") |>
    dplyr::mutate(model = "clmm", cond_H = summary(fnc)$info$cond.H, used = nc_ok),
  wald(fnc0, "x") |> dplyr::mutate(model = "clm", cond_H = summary(fnc0)$info$cond.H,
                                   used = !nc_ok)) |>
  dplyr::mutate(epochs = nrow(pnc), participants = dplyr::n_distinct(pnc$subject))
nc_bins <- pnc |> dplyr::mutate(bin = cut(x, c(-Inf, -STRICT_X, STRICT_X, Inf),
                                          labels = c("measured faster (x < -.20)", "|x| <= .20",
                                                     "measured slower (x > .20)"))) |>
  dplyr::count(bin, response) |> tidyr::pivot_wider(names_from = response, values_from = n,
                                                     values_fill = 0)
cat("\nNO CHANGE EPOCHS:\n"); print(as.data.frame(ncw), digits = 3); print(as.data.frame(nc_bins))
readr::write_csv(dplyr::bind_rows(ncw |> dplyr::mutate(kind = "model", .before = 1),
                                  nc_bins |> dplyr::mutate(kind = "counts", .before = 1)),
                 file.path(rd, "EXPLORATORY_E14_nochange.csv"))

# 3b. Dissociation epochs: breathing on the wrong side of 1 ---------
pch <- p |> dplyr::filter(change_direction != "No Change") |>
  dplyr::mutate(opp = ifelse(change_direction == "Faster", "Slower", "Faster"),
                follow = dplyr::case_when(response == change_direction ~ "pacer",
                                          response == "No Change" ~ "No Change",
                                          TRUE ~ "breath_or_wrong"),
                wrong_side = sign(x) == -sign(z) & x != 0,
                set = dplyr::case_when(wrong_side & abs(x) >= STRICT_X ~ "dissociated_strict",
                                       wrong_side ~ "dissociated_mild",
                                       TRUE ~ "concordant"))
ppacer <- function(d) { dd <- d[d$follow != "No Change", ]; mean(dd$follow == "pacer") }
sets <- list(concordant = pch |> dplyr::filter(set == "concordant"),
             dissociated = pch |> dplyr::filter(wrong_side),
             dissociated_strict = pch |> dplyr::filter(set == "dissociated_strict"))
dis <- purrr::imap_dfr(sets, function(d, nm) {
  ci <- pboot(d, ppacer)
  tibble::tibble(set = nm, epochs = nrow(d), participants = dplyr::n_distinct(d$subject),
                 n_follow_pacer = sum(d$follow == "pacer"), n_no_change = sum(d$follow == "No Change"),
                 n_follow_breath = sum(d$follow == "breath_or_wrong"),
                 p_pacer_given_directional = ppacer(d), lo = ci[1], hi = ci[2],
                 p_no_change = mean(d$follow == "No Change"))
})
dd_fun <- function(lab) function(d) ppacer(d[d$wrong_side & (lab == "any" | abs(d$x) >= STRICT_X), ]) -
  ppacer(d[d$set == "concordant", ])
for (lab in c("any", "strict")) {
  ci <- pboot(pch, dd_fun(lab))
  dis <- dplyr::bind_rows(dis, tibble::tibble(
    set = paste0("difference_", ifelse(lab == "any", "dissociated", "dissociated_strict"),
                 "_minus_concordant"),
    p_pacer_given_directional = dd_fun(lab)(pch), lo = ci[1], hi = ci[2]))
}
cat("\nDISSOCIATION EPOCHS:\n"); print(as.data.frame(dis), digits = 3)
readr::write_csv(dis, file.path(rd, "EXPLORATORY_E14_dissociation.csv"))

# 4. Free epochs ---------
pz <- m |> dplyr::filter(trial_type == "paced_breathing") |>
  dplyr::select(subject, bids_entity, epoch_index, p_measured = measured,
                cue = change_amount) |>
  dplyr::mutate(cue = log(cue))
fr <- m |> dplyr::filter(trial_type == "free_breathing", measured %in% TRUE, !is.na(response),
                         subject != STALE_CACHE) |>
  dplyr::left_join(pz, by = c("subject", "bids_entity", "epoch_index")) |>
  dplyr::left_join(bc |> dplyr::filter(trial_type == "free_breathing") |>
                     dplyr::select(subject, bids_entity, epoch_index, f_n = n_cont, log_med),
                   by = c("subject", "bids_entity", "epoch_index")) |>
  dplyr::left_join(bc |> dplyr::filter(trial_type == "paced_breathing") |>
                     dplyr::select(subject, bids_entity, epoch_index, p_n = n_cont, log_last2),
                   by = c("subject", "bids_entity", "epoch_index")) |>
  dplyr::filter(p_measured %in% TRUE, f_n >= 2, p_n >= 2) |>
  dplyr::mutate(trans = log_med - log_last2, within = log(actual_ratio),
                resp = factor(response, levels = LV, ordered = TRUE), subject = factor(subject))
frw <- fr |> dplyr::filter(is.finite(within)) |> dplyr::mutate(subject = droplevels(subject))
ff1 <- fit_clmm(resp ~ trans + (1 | subject), fr)
ff2 <- fit_clmm(resp ~ trans + within + (1 | subject), frw)
ff3 <- fit_clmm(resp ~ trans + cue + (1 | subject), fr)
free <- dplyr::bind_rows(
  wald(ff1, "trans") |> dplyr::mutate(model = "F1 trans", epochs = nrow(fr),
                                      participants = dplyr::n_distinct(fr$subject)),
  wald(ff2, c("trans", "within")) |> dplyr::mutate(model = "F2 trans + within", epochs = nrow(frw),
                                                    participants = dplyr::n_distinct(frw$subject)),
  wald(ff3, c("trans", "cue")) |> dplyr::mutate(model = "F3 trans + preceding pacer endpoint",
                                                 epochs = nrow(fr),
                                                 participants = dplyr::n_distinct(fr$subject)))
free_desc <- fr |> dplyr::group_by(response) |>
  dplyr::summarise(n = dplyr::n(), trans_median = median(trans), .groups = "drop")
cat("\nFREE EPOCHS:\n"); print(as.data.frame(free), digits = 3); print(as.data.frame(free_desc))
cat(sprintf("  cor(trans, cue) = %.3f ; cor(trans, within) = %.3f\n", cor(fr$trans, fr$cue),
            cor(frw$trans, frw$within)))
readr::write_csv(dplyr::bind_rows(free, free_desc |> dplyr::mutate(model = "descriptive")),
                 file.path(rd, "EXPLORATORY_E14_free.csv"))

# 5. Confidence ---------
pcf <- p |> dplyr::filter(is.finite(confidence)) |> dplyr::mutate(az = abs(z), ax = abs(x))
fc <- lmerTest::lmer(confidence ~ az + ax + sal + (1 | subject), data = pcf, REML = TRUE)
cs <- coef(summary(fc))
conf <- tibble::tibble(term = rownames(cs), est = cs[, 1], se = cs[, 2], df = cs[, 3],
                       lo = cs[, 1] - qt(.975, cs[, 3]) * cs[, 2],
                       hi = cs[, 1] + qt(.975, cs[, 3]) * cs[, 2], t = cs[, 4], p = cs[, 5],
                       epochs = nrow(pcf), participants = dplyr::n_distinct(pcf$subject))
cat("\nCONFIDENCE:\n"); print(as.data.frame(conf), digits = 3)
readr::write_csv(conf, file.path(rd, "EXPLORATORY_E14_confidence.csv"))

# 6. Sensitivity ---------
two <- p |> dplyr::group_by(subject) |> dplyr::filter(dplyr::n_distinct(bids_entity) == 2) |>
  dplyr::ungroup() |> dplyr::mutate(subject = droplevels(subject))
fs1 <- fit_clmm(resp ~ z + x + sal + prev_resp + first + (1 | subject), p)
fs2 <- fit_clmm(resp ~ z + x + sal + (1 | subject), two)
sens <- dplyr::bind_rows(
  wald(fs1, c("z", "x", "sal", "prev_resp", "first")) |>
    dplyr::mutate(model = "S1 + previous paced report", epochs = n_ep_p,
                  participants = n_sub_p),
  wald(fs2, c("z", "x", "sal")) |>
    dplyr::mutate(model = "S2 two-run participants", epochs = nrow(two),
                  participants = dplyr::n_distinct(two$subject)))
cat("\nSENSITIVITY:\n"); print(as.data.frame(sens), digits = 3)
readr::write_csv(sens, file.path(rd, "EXPLORATORY_E14_sensitivity.csv"))

counts <- dplyr::bind_rows(counts,
  tibble::tibble(set = "paced_nochange", epochs = nrow(pnc),
                 participants = dplyr::n_distinct(pnc$subject), runs = nrun(pnc)),
  tibble::tibble(set = "paced_change", epochs = nrow(pch),
                 participants = dplyr::n_distinct(pch$subject), runs = nrun(pch)),
  tibble::tibble(set = "free_transition", epochs = nrow(fr),
                 participants = dplyr::n_distinct(fr$subject), runs = nrun(fr)),
  tibble::tibble(set = "free_transition_and_within", epochs = nrow(frw),
                 participants = dplyr::n_distinct(frw$subject), runs = nrun(frw)),
  tibble::tibble(set = "two_run_participants", epochs = nrow(two),
                 participants = dplyr::n_distinct(two$subject), runs = nrun(two)),
  tibble::tibble(set = "trend_se_nochange_epochs", epochs = nrow(se_nc)))
readr::write_csv(counts, file.path(rd, "EXPLORATORY_E14_counts.csv"))
print(as.data.frame(counts))

# Figure ---------
{
  col_p <- "#2a78d6"; col_b <- "#eb6834"; col_n <- "#b4b2ad"; ink <- "#52514e"
  png(fig, width = 1650, height = 560, res = 130)
  op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 3, 1), col.axis = ink, col.lab = ink, fg = ink)
  # A: per-SD coefficients
  est <- prim$est_per_sd[1:2]; lo <- prim$lo_per_sd[1:2]; hi <- prim$hi_per_sd[1:2]
  plot(c(1, 2), est, xlim = c(0.5, 2.5), ylim = range(0, lo, hi) * 1.15, pch = 19, cex = 1.6,
       col = c(col_p, col_b), xaxt = "n", xlab = "", ylab = "probit coefficient per SD",
       main = "A. Report ~ pacer + breath (paced)", bty = "n")
  segments(1:2, lo, 1:2, hi, lwd = 2, col = c(col_p, col_b)); abline(h = 0, col = col_n)
  axis(1, 1:2, c("instructed\n(visible pacer)", "measured\n(own breathing)"), padj = 0.5)
  # B: dissociation
  dz <- dis |> dplyr::filter(!grepl("difference", set))
  mat <- t(as.matrix(dz[, c("n_follow_pacer", "n_no_change", "n_follow_breath")]) / dz$epochs)
  barplot(mat, names.arg = c("breath with\npacer", "breath against\npacer", "against,\n|x| >= .20"),
          col = c(col_p, col_n, col_b), border = "white", ylab = "proportion of epochs",
          main = "B. Change epochs: what was reported", ylim = c(0, 1.15), yaxt = "n")
  axis(2, seq(0, 1, .2))
  legend("top", c("pacer direction", "No Change", "opposite to pacer"), fill = c(col_p, col_n, col_b),
         border = "white", bty = "n", horiz = TRUE, cex = 0.85)
  # C: free epochs
  fr$bin <- cut(fr$trans, quantile(fr$trans, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                labels = c("faster", "middle", "slower"))
  col_f <- "#1baf7a"; col_s <- "#6250d6"
  tb <- prop.table(table(fr$response, fr$bin), 2)[LV, ]
  barplot(tb, col = c(col_f, col_n, col_s), border = "white", ylab = "proportion of epochs",
          xlab = "free vs end of paced: measured period (tertile)",
          main = "C. Free epochs: report by transition", ylim = c(0, 1.15), yaxt = "n")
  axis(2, seq(0, 1, .2))
  legend("top", paste("report", LV), fill = c(col_f, col_n, col_s), border = "white", bty = "n",
         horiz = TRUE, cex = 0.85)
  par(op); dev.off()
  cat("figure written\n")
}
