# intero2023_physio_sheets.R
# Rewrites ONLY the physiological columns of the committed first-level
# confound sheets, from Data/regressors, and updates has_rvt / has_petco2 /
# rvt_source in the two manifests. Added 2026-11-24.
#
# WHY A SEPARATE WRITER. Regenerating the sheets with
# intero2023_covariates_linux.R would also re-derive csf / white_matter from the
# Linux fMRIPrep pass, which differ from the committed sheets (VERIFICATION_REPORT
# §2: median r 0.985 / 0.996). The physio fixes must not touch those columns, so
# this script carries every non-physio column through as the exact text that is
# committed, and replaces only:
#
#   rvt_rrf_z, rvt_z, resp_rate_c                  (belt; from <sub>_<ent>_resp.tsv,
#                                                   or the PMU fallback file)
#   petco2_hrf_lag00, _lag10, _lag20, _lag30        (from <sub>_<ent>_petco2.tsv)
#
# using the sheets' own convention (FirstLevel/README.md; identical code to
# intero2023_covariates_linux.R): each column is padded/trimmed to n_vol,
# demeaned within run, NA -> 0, rounded to 6 dp; a belt column whose coverage is
# < 50% is omitted; the PetCO2 basis enters whole or not at all.
#
# RVT SOURCE (see intero2023_pmu_rvt.R): Biopac if <sub>_<ent>_resp.tsv exists
# and its rvt_rrf_z covers >= 50% of volumes; else <sub>_<ent>_desc-pmu_resp.tsv
# if present; else none. Recorded per run as rvt_source = biopac | pmu | none.
#
# A sheet whose physio columns come out identical to the committed text is not
# rewritten. The script prints every sheet it rewrites.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble); library(purrr) })

dir_repo <- Sys.getenv("INTERO_REPO", "/media/normega/FB_Cluster_Data/intero2024/MRI_Intero2023")
dir_fl   <- file.path(dir_repo, "FirstLevel")
dir_reg  <- file.path(dir_repo, "Data", "regressors")
PHYSIO_MIN_COV <- 0.5
BELT <- c(rvt_rrf_z = "rvt_rrf_z", rvt_z = "rvt_z", resp_rate_c = "resp_rate_brpm")
CO2  <- sprintf("petco2_hrf_lag%02d", c(0, 10, 20, 30))

demean0 <- function(x) {                 # verbatim from intero2023_covariates_linux.R
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0, length(x)))
  x[ok] <- x[ok] - mean(x[ok]); x[!ok] <- 0; x
}
fit_len <- function(x, n) { length(x) <- n; x }
read_reg <- function(f) read_tsv(f, na = "n/a", show_col_types = FALSE, progress = FALSE)
usable_biopac <- function(f, n_vol)
  file.exists(f) && mean(is.finite(fit_len(read_reg(f)$rvt_rrf_z, n_vol))) >= PHYSIO_MIN_COV

physio_block <- function(sb, ent, n_vol) {
  out <- list(); src <- "none"
  bf <- file.path(dir_reg, sprintf("%s_%s_resp.tsv", sb, ent))
  pf <- file.path(dir_reg, sprintf("%s_%s_desc-pmu_resp.tsv", sb, ent))
  rf <- if (usable_biopac(bf, n_vol)) bf else if (file.exists(pf)) pf else if (file.exists(bf)) bf else NA
  if (!is.na(rf)) {
    r <- read_reg(rf)
    for (nm in names(BELT)) {
      v <- fit_len(r[[BELT[[nm]]]], n_vol)
      if (mean(is.finite(v)) >= PHYSIO_MIN_COV) out[[nm]] <- demean0(v)
    }
    if (length(out)) src <- if (rf == pf) "pmu" else "biopac"
  }
  cf <- file.path(dir_reg, sprintf("%s_%s_petco2.tsv", sb, ent))
  has_co2 <- FALSE
  if (file.exists(cf)) {
    p <- read_reg(cf)
    if (all(CO2 %in% names(p)) &&
        all(vapply(CO2, function(cc) mean(is.finite(fit_len(p[[cc]], n_vol))) >= PHYSIO_MIN_COV, TRUE))) {
      for (cc in CO2) out[[cc]] <- demean0(fit_len(p[[cc]], n_vol)); has_co2 <- TRUE
    }
  }
  list(cols = as_tibble(out), rvt_source = src, has_rvt = src != "none", has_co2 = has_co2)
}

# physio values -> the exact text readr writes in the sheets (round 6, write_tsv)
as_sheet_text <- function(tb) {
  if (!ncol(tb)) return(NULL)            # no physio for this run: bind_cols() drops NULL
  txt <- format_tsv(tb %>% mutate(across(everything(), ~ round(.x, 6))), na = "n/a")
  read_tsv(I(txt), col_types = cols(.default = col_character()), na = character(),
           show_col_types = FALSE, progress = FALSE)
}

update_manifest <- function(mfile) {
  m <- read_csv(mfile, col_types = cols(.default = col_character()), na = character(),
                show_col_types = FALSE, progress = FALSE)
  if (!"rvt_source" %in% names(m)) m <- m %>% mutate(rvt_source = NA_character_, .after = has_petco2)
  for (i in seq_len(nrow(m))) {
    sb <- m$subject[i]; ent <- m$bids_entity[i]; n_vol <- as.integer(m$n_vol[i])
    sh <- file.path(dir_fl, sprintf("%s_%s_desc-model_confounds.tsv", sb, ent))
    d  <- read_tsv(sh, col_types = cols(.default = col_character()), na = character(),
                   show_col_types = FALSE, progress = FALSE)
    stopifnot(nrow(d) == n_vol)
    pb <- physio_block(sb, ent, n_vol)
    new <- bind_cols(d %>% select(-any_of(c(names(BELT), CO2))), as_sheet_text(pb$cols))
    if (!identical(names(new), names(d)) || !identical(as.list(new), as.list(d[names(new)]))) {
      write_tsv(new, sh, na = "n/a")
      cat(sprintf("  rewrote %s_%s  [%s | rvt_source=%s, petco2=%s]\n", sb, ent,
                  paste(intersect(c(names(BELT), CO2), names(new)), collapse = ","), pb$rvt_source, pb$has_co2))
    }
    m$has_rvt[i] <- toupper(as.character(pb$has_rvt)); m$has_petco2[i] <- toupper(as.character(pb$has_co2))
    m$rvt_source[i] <- pb$rvt_source
  }
  write_csv(m, mfile, na = "")
  cat(sprintf("%s: %d runs; rvt_source %s\n", basename(mfile), nrow(m),
              paste(names(table(m$rvt_source)), table(m$rvt_source), sep = "=", collapse = " ")))
}

for (mf in c("firstlevel_manifest.csv", "firstlevel_manifest_candidates.csv"))
  update_manifest(file.path(dir_fl, mf))
