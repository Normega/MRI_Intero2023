# intero2023_loc_t4_mixed.R -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §8 T4):
# mixed models of the per-epoch C1-pattern expression index.
# Choices: Analysis/localizer_stageC3_implementation_notes.md §1 (committed before any value).
#
#   Rscript Scripts/intero2023_loc_t4_mixed.R
#
# Input : Analysis/results/EXPLORATORY_LOC_C3_T4_epochs.csv (Scripts/intero2023_loc_t4.py index)
# Output: Analysis/results/EXPLORATORY_LOC_C3_T4_mixed.csv
# Registered form (plan §8): index ~ confidence + outcome + (1 | participant), paced epochs,
# outcome = aware / unaware / nochange (wrong-direction epochs are not one of the three outcomes).
# lme4 REML, lmerTest Satterthwaite df; outcome contrasts with emmeans (Satterthwaite).

# Set Up ---------
## Load libraries ---------
packages <- c("tidyverse", "lme4", "lmerTest", "emmeans")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}
emm_options(lmer.df = "satterthwaite")

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
rd <- file.path(dirname(here), "Analysis", "results")
d0 <- readr::read_csv(file.path(rd, "EXPLORATORY_LOC_C3_T4_epochs.csv")) |>
  dplyr::mutate(subject = factor(subject), run_id = factor(paste(subject, run)))

INDICES <- c("dot_own", "corr_own", "dot_loso", "corr_loso",
             "dot_own_SAL", "dot_own_SMN", "corr_own_SAL", "corr_own_SMN")

fixed_rows <- function(m, idx, model, data) {
  s <- summary(m)$coefficients
  ci <- confint(m, method = "Wald", parm = "beta_")
  tibble::tibble(index = idx, model = model, term = rownames(s), estimate = s[, "Estimate"],
                 se = s[, "Std. Error"], df = s[, "df"], t = s[, "t value"], p = s[, "Pr(>|t|)"],
                 ci_lo = ci[rownames(s), 1], ci_hi = ci[rownames(s), 2],
                 n_participants = dplyr::n_distinct(data$subject), n_epochs = nrow(data))
}

contrast_rows <- function(m, idx, model, data) {
  em <- emmeans::emmeans(m, ~ outcome)
  pr <- summary(pairs(em, reverse = FALSE, adjust = "none"), infer = c(TRUE, TRUE))
  tibble::tibble(index = idx, model = model, term = paste0("contrast: ", pr$contrast),
                 estimate = pr$estimate, se = pr$SE, df = pr$df, t = pr$t.ratio, p = pr$p.value,
                 ci_lo = pr$lower.CL, ci_hi = pr$upper.CL,
                 n_participants = dplyr::n_distinct(data$subject), n_epochs = nrow(data))
}

omnibus_row <- function(m, idx, model, data, term = "outcome") {
  a <- anova(m, type = 3, ddf = "Satterthwaite")
  tibble::tibble(index = idx, model = model, term = paste0("F: ", term),
                 estimate = a[term, "F value"], df = a[term, "DenDF"], p = a[term, "Pr(>F)"],
                 n_participants = dplyr::n_distinct(data$subject), n_epochs = nrow(data))
}

icc_row <- function(data, idx, model) {
  m <- lmer(y ~ 1 + (1 | subject), data = data, REML = TRUE)
  vc <- as.data.frame(VarCorr(m))
  tibble::tibble(index = idx, model = model, term = "ICC (participant intercept share)",
                 estimate = vc$vcov[1] / sum(vc$vcov),
                 se = NA_real_, sd_between = sqrt(vc$vcov[1]), sd_within = sqrt(vc$vcov[2]),
                 n_participants = dplyr::n_distinct(data$subject), n_epochs = nrow(data))
}

out <- list()
for (idx in INDICES) {
  paced <- d0 |>
    dplyr::filter(epoch_type == "paced", outcome %in% c("aware", "unaware", "nochange"),
                  !is.na(confidence), !is.na(.data[[idx]])) |>
    dplyr::mutate(y = .data[[idx]], outcome = factor(outcome, levels = c("aware", "unaware", "nochange")),
                  conf_b = ave(confidence, subject), conf_w = confidence - conf_b,
                  conf_b = conf_b - mean(conf_b))
  # registered form
  m1 <- lmer(y ~ confidence + outcome + (1 | subject), data = paced, REML = TRUE)
  out[[length(out) + 1]] <- fixed_rows(m1, idx, "registered", paced)
  out[[length(out) + 1]] <- contrast_rows(m1, idx, "registered", paced)
  out[[length(out) + 1]] <- omnibus_row(m1, idx, "registered", paced)
  out[[length(out) + 1]] <- icc_row(paced, idx, "registered")
  if (idx %in% c("dot_own", "corr_own")) {
    # sensitivity: run intercepts nested in participant (the breathing runs have no clean baseline)
    m2 <- lmer(y ~ confidence + outcome + (1 | subject) + (1 | run_id), data = paced, REML = TRUE)
    out[[length(out) + 1]] <- fixed_rows(m2, idx, "run_intercepts", paced)
    out[[length(out) + 1]] <- contrast_rows(m2, idx, "run_intercepts", paced)
    # sensitivity: confidence split into within- and between-participant parts
    m3 <- lmer(y ~ conf_w + conf_b + outcome + (1 | subject), data = paced, REML = TRUE)
    out[[length(out) + 1]] <- fixed_rows(m3, idx, "confidence_within_between", paced)
    # sensitivity: wrong-direction epochs as a fourth outcome
    p4 <- d0 |>
      dplyr::filter(epoch_type == "paced", !is.na(confidence), !is.na(.data[[idx]])) |>
      dplyr::mutate(y = .data[[idx]],
                    outcome = factor(outcome, levels = c("aware", "unaware", "nochange", "wrongdir")))
    m4 <- lmer(y ~ confidence + outcome + (1 | subject), data = p4, REML = TRUE)
    out[[length(out) + 1]] <- contrast_rows(m4, idx, "with_wrongdir", p4)
    # added descriptive check: paced vs free epochs
    pf <- d0 |>
      dplyr::filter(!is.na(.data[[idx]]), outcome != "wrongdir") |>
      dplyr::mutate(y = .data[[idx]], epoch_type = factor(epoch_type, levels = c("free", "paced")))
    m5 <- lmer(y ~ epoch_type + (1 | subject) + (1 | run_id), data = pf, REML = TRUE)
    out[[length(out) + 1]] <- fixed_rows(m5, idx, "paced_vs_free", pf)
  }
}
res <- dplyr::bind_rows(out)
readr::write_csv(res |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 6))),
                 file.path(rd, "EXPLORATORY_LOC_C3_T4_mixed.csv"))
print(res |> dplyr::filter(index %in% c("dot_own", "corr_own")) |>
        dplyr::select(index, model, term, estimate, ci_lo, ci_hi, p) |> as.data.frame(), digits = 3)
