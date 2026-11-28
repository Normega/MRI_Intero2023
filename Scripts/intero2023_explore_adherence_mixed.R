# intero2023_explore_adherence_mixed.R -- EXPLORATORY E12 Step 3 (deviations.md, "After the update", 2026-11-28;
# run after the PI's post-hoc override of the failed Step 0 gate, as the override row specifies):
# epoch-level mixed model of paced PLV on the T4 C1-pattern expression index.
#
#   Rscript Scripts/intero2023_explore_adherence_mixed.R
#
# Inputs : Analysis/results/EXPLORATORY_E12_epochs.csv (paced PLV, change type, direction)
#          Analysis/results/EXPLORATORY_LOC_C3_T4_epochs.csv (dot_own index per epoch)
# Output : Analysis/results/EXPLORATORY_E12_mixed.csv
# Model  : plv ~ dot_own_c + change_type + dir + (1 | subject), REML, Satterthwaite df, two-sided.
#          dot_own_c = dot_own centred within participant x run; change_type = step / ramp / none (factor);
#          dir = +1/2 Faster, -1/2 Slower, 0 none.

# Set Up ---------
## Load libraries ---------
packages <- c("tidyverse", "lme4", "lmerTest")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
rd <- file.path(dirname(here), "Analysis", "results")

adh <- readr::read_csv(file.path(rd, "EXPLORATORY_E12_epochs.csv")) |>
  dplyr::filter(!is.na(plv)) |>
  dplyr::select(subject, entity = bids_entity, epoch_index, plv, change_type, direction)
t4 <- readr::read_csv(file.path(rd, "EXPLORATORY_LOC_C3_T4_epochs.csv")) |>
  dplyr::filter(epoch_type == "paced") |>
  dplyr::mutate(epoch_index = as.integer(epoch_index)) |>
  dplyr::select(subject, entity, run, epoch_index, dot_own)

d <- dplyr::inner_join(adh, t4, by = c("subject", "entity", "epoch_index")) |>
  dplyr::group_by(subject, run) |>
  dplyr::mutate(dot_own_c = dot_own - mean(dot_own)) |>
  dplyr::ungroup() |>
  dplyr::mutate(change_type = factor(change_type, levels = c("none", "ramp", "step")),
                dir = dplyr::case_when(direction == "Faster" ~ 0.5, direction == "Slower" ~ -0.5, TRUE ~ 0),
                subject = factor(subject))

m <- lmer(plv ~ dot_own_c + change_type + dir + (1 | subject), data = d, REML = TRUE)
s <- summary(m)$coefficients
ci <- confint(m, method = "Wald", parm = "beta_")
out <- tibble::tibble(analysis = "EXPLORATORY E12 Step 3 (post-hoc override of the gate)",
                      term = rownames(s), estimate = s[, "Estimate"], se = s[, "Std. Error"], df = s[, "df"],
                      t = s[, "t value"], p = s[, "Pr(>|t|)"], ci_lo = ci[rownames(s), 1], ci_hi = ci[rownames(s), 2],
                      n_participants = dplyr::n_distinct(d$subject), n_epochs = nrow(d),
                      sd_index_within = sd(d$dot_own_c), sd_plv = sd(d$plv))
readr::write_csv(out, file.path(rd, "EXPLORATORY_E12_mixed.csv"))
print(out, width = 200)
