# intero2023_tissue_confounds.R
# Rewrites ONLY the csf and white_matter columns of the committed first-level
# confound sheets, from the fMRIPrep derivatives of record. Added 2026-11-25,
# before the OSF freeze and before any model was fitted.
#
# WHY. VERIFICATION_REPORT.md §2a: in the 95 registered sheets these two columns
# do not reproduce from the derivatives of record (median r 0.985 / 0.996); every
# other column does, bit-exactly. They came from an earlier anatomical pass that
# cannot be inspected. This script makes the tissue regressors match the stated
# inputs. The nuisance set itself is unchanged (csf and white_matter were always
# in it).
#
# HOW. Every other column is carried through as the exact committed text. The
# two columns are taken from <sub>_<ent>_desc-confounds_timeseries.tsv under the
# sheets' own convention (intero2023_covariates_linux.R): numeric, NA -> 0,
# rounded to 6 dp, written by readr. No BOLD voxel is read.
#
# IDENTITY CHECK. A sheet is rewritten only if its row count equals the fMRIPrep
# file's and its 12 motion columns and FD match the fMRIPrep values to 6 dp;
# otherwise the script stops. The manifest is not touched (no gate depends on
# csf or white_matter).
#
# Output: prints one line per sheet (r old vs new for each column) and writes
# Analysis/tissue_confounds_regen.csv.
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble); library(purrr) })

dir_repo <- Sys.getenv("INTERO_REPO", "/media/normega/FB_Cluster_Data/intero2024/MRI_Intero2023")
dir_fmri <- Sys.getenv("INTERO_FMRIPREP", "/media/normega/FB_Cluster_Data/intero2024/derivatives/fmriprep")
dir_fl   <- file.path(dir_repo, "FirstLevel")
TISSUE   <- c("csf", "white_matter")
MOTION   <- c(paste0(rep(c("trans", "rot"), each = 3), "_", c("x", "y", "z")))
CHECK    <- c(MOTION, paste0(MOTION, "_derivative1"), "framewise_displacement")

num0 <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.finite(x), x, 0) }
as_sheet_text <- function(tb) {          # identical to intero2023_physio_sheets.R
  txt <- format_tsv(tb %>% mutate(across(everything(), ~ round(.x, 6))), na = "n/a")
  read_tsv(I(txt), col_types = cols(.default = col_character()), na = character(),
           show_col_types = FALSE, progress = FALSE)
}

runs <- bind_rows(lapply(c("firstlevel_manifest.csv", "firstlevel_manifest_candidates.csv"), function(f)
  read_csv(file.path(dir_fl, f), col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
    select(subject, bids_entity))) %>% distinct()

res <- pmap_dfr(runs, function(subject, bids_entity) {
  sh <- file.path(dir_fl, sprintf("%s_%s_desc-model_confounds.tsv", subject, bids_entity))
  if (!file.exists(sh)) return(NULL)
  fp <- file.path(dir_fmri, subject, "func",
                  sprintf("%s_%s_desc-confounds_timeseries.tsv", subject, bids_entity))
  if (!file.exists(fp)) stop("no fMRIPrep confounds for ", subject, " ", bids_entity)
  d <- read_tsv(sh, col_types = cols(.default = col_character()), na = character(),
                show_col_types = FALSE, progress = FALSE)
  f <- read_tsv(fp, col_types = cols(.default = col_character()), show_col_types = FALSE, progress = FALSE)
  stopifnot(nrow(d) == nrow(f), all(TISSUE %in% names(d)))
  for (cc in CHECK) if (max(abs(num0(d[[cc]]) - round(num0(f[[cc]]), 6))) > 1e-6)
    stop("identity check failed: ", subject, " ", bids_entity, " column ", cc)
  newcols <- as_sheet_text(tibble(csf = num0(f$csf), white_matter = num0(f$white_matter)))
  r_csf <- cor(num0(d$csf), num0(newcols$csf)); r_wm <- cor(num0(d$white_matter), num0(newcols$white_matter))
  changed <- !identical(d$csf, newcols$csf) || !identical(d$white_matter, newcols$white_matter)
  if (changed) {
    d$csf <- newcols$csf; d$white_matter <- newcols$white_matter
    write_tsv(d, sh, na = "n/a")
  }
  tibble(subject, bids_entity, n_vol = nrow(d), rewritten = changed,
         r_csf = round(r_csf, 4), r_white_matter = round(r_wm, 4))
})

write_csv(res, file.path(dir_repo, "Analysis", "tissue_confounds_regen.csv"))
cat(sprintf("%d sheets checked, %d rewritten. r(old,new) csf median %.3f min %.3f; white_matter median %.3f min %.3f\n",
            nrow(res), sum(res$rewritten), median(res$r_csf), min(res$r_csf),
            median(res$r_white_matter), min(res$r_white_matter)))
