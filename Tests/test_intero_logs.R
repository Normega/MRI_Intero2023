# test_intero_logs.R
# Validation for the Intero2023 log parser. Plain base-R assertions so this
# runs with no extra packages; port to testthat if it gets installed.
#
#   Rscript Repo/Tests/test_intero_logs.R
#
# The response-recovery tests are the important ones: the off-by-one and
# carry-forward behaviour of these logs is the failure mode most likely to
# silently corrupt the behavioural analysis.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(purrr); library(tibble)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_R    <- file.path(dir_root, "Repo", "R")
dir_logs <- file.path(dir_root, "Data", "Logs")
dir_vals <- file.path(dir_root, "Data", "Values")

source(file.path(dir_R, "intero_ids.R"))
source(file.path(dir_R, "intero_logs.R"))

.pass <- 0; .fail <- 0
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) {
    cat(sprintf("   error: %s\n", conditionMessage(e))); FALSE
  }))
  if (ok) { .pass <<- .pass + 1; cat(sprintf("  PASS  %s\n", label)) }
  else    { .fail <<- .fail + 1; cat(sprintf("  FAIL  %s\n", label)) }
}

# ---- 1. ID normalisation --------------------------------------------------
cat("\n-- ID normalisation --\n")
check("all spellings of participant 52 agree",
      all(pid_num(c("P52", "P052", "sub-052", "Intero_P52.acq",
                    "P052_20231025", 52, "52")) == 52L))
check("dates in filenames are not mistaken for IDs",
      pid_num("P052_20231025") == 52L)
check("unpadded P90 and padded P090 collapse",
      pid_num("P90") == pid_num("P090"))
check("BIDS label is zero-padded", pid_bids("P52") == "sub-052")
check("unparseable input yields NA", is.na(pid_num("no-id-here")))
check("pilot flags are exactly P00/P01/P50/P999",
      all(is_pilot(c("P50", "P999", "P01", "P00"))) &&
      !any(is_pilot(c("P51", "P60", "P130"))))
check("protocol version splits at P60",
      protocol_version("P59") == "v1" && protocol_version("P60") == "v2")
check("run labels map to the Linux-copy BIDS naming",
      bids_entities("ChangeDetection1") == "task-breathing_run-1" &&
      bids_entities("Authenticity")     == "task-authenticity" &&
      bids_entities("Localizer2")       == "task-localizer_run-2")

# ---- 2. duration parsing --------------------------------------------------
cat("\n-- duration parsing --\n")
check("H:MM:SS.ffffff parses to seconds",
      abs(parse_hms("0:25:55.433102") - 1555.433102) < 1e-6)
check("H:MM:SS with no fraction parses", parse_hms("0:00:00") == 0)
check("blank and malformed yield NA, not 0",
      all(is.na(parse_hms(c("", NA, "garbage", "12")))))

# ---- 3. the off-by-one, against known ground truth ------------------------
# P130's Authenticity run: the 16 value words come from Values/P130.txt, and
# each word's rating is written on the FOLLOWING row -- the last one landing on
# the trailing Completion row. If the shift is dropped or doubled, the count of
# recovered ratings changes and this test fails.
cat("\n-- response recovery: P130 Authenticity --\n")

f130 <- file.path(dir_logs,
                  "Log.Intero2023.P130.SessionTest_2024-06-03_11h02.37.251.csv")
d130 <- read_log_file(f130)
auth <- d130 %>% filter(Run == "Authenticity") %>% shift_responses()
stim <- auth %>% filter(Description == "AuthenticRating_Start")
vals <- readLines(file.path(dir_vals, "P130.txt"), warn = FALSE)
vals <- trimws(vals[nzchar(trimws(vals))])

check("16 value words presented", nrow(stim) == 16)
check("presented stimuli are exactly the Values file, as a set",
      setequal(trimws(stim$Stimulus), vals))
check("every stimulus recovers a rating (none lost to the shift)",
      all(!is.na(stim$Authenticity_resp)))
check("recovered ratings are in the 1-7 response range",
      all(as.numeric(stim$Authenticity_resp) %in% 1:7))
check("first stimulus recovers the value logged on the NEXT row",
      stim$Authenticity_resp[1] == auth$Authenticity[which(
        auth$Description == "AuthenticRating_Start")[1] + 1])
check("last stimulus recovers the value from the Completion row",
      stim$Authenticity_resp[16] ==
        auth$Authenticity[auth$Description == "Completion"][1])
check("the raw column would have LOST the last rating",
      is.na(stim$Authenticity[1]) || trimws(stim$Authenticity[1]) == "")

# ---- 4. staleness must not be read as new responses -----------------------
# In the breathing runs RT is populated and carries forward unchanged into
# rows where nothing was answered. Recovery must not credit those rows.
cat("\n-- staleness: P130 breathing run --\n")

cd <- d130 %>% filter(Run == "ChangeDetection1") %>% shift_responses()
rt_raw <- suppressWarnings(as.numeric(cd$RT))
n_nonblank_raw <- sum(is.finite(rt_raw))
n_recovered    <- sum(!is.na(cd$RT_resp))
n_distinct_rt  <- length(unique(rt_raw[is.finite(rt_raw)]))

cat(sprintf("   raw non-blank RT rows=%d, distinct RT values=%d, recovered=%d\n",
            n_nonblank_raw, n_distinct_rt, n_recovered))
check("recovered responses are far fewer than non-blank raw rows",
      n_recovered < n_nonblank_raw)
check("recovered count equals the number of distinct RT values",
      n_recovered == n_distinct_rt)
check("no recovered RT is a duplicate of its predecessor",
      { r <- suppressWarnings(as.numeric(cd$RT_resp)); r <- r[is.finite(r)]
        length(r) == length(unique(r)) })

# ---- 5. duplicate / empty file handling -----------------------------------
cat("\n-- file hygiene --\n")
scan <- scan_log_dir(dir_logs)
check("the two byte-identical content duplicates are detected",
      sum(scan$files$is_duplicate) == 2)
check("duplicates are exactly (1)-suffixed files, and never empty ones",
      all(scan$files$is_dup_suffix[scan$files$is_duplicate]) &&
      !any(scan$files$is_empty[scan$files$is_duplicate]))
check("all 18 zero-byte launches are marked empty, not duplicates",
      sum(scan$files$is_empty) == 18 &&
      !any(scan$files$is_duplicate & scan$files$is_empty))
check("empty files contribute no runs",
      !any(scan$attempts$file %in% scan$files$file[scan$files$is_empty]))
check("every filename parsed", all(scan$files$parsed_ok))
check("the pilot staircase run is not counted as a scanner run",
      !"ChangeThreshold" %in% INTERO_RUNS$log_run &&
      all(!scan$attempts$is_scanner_run[scan$attempts$log_run == "ChangeThreshold"]))

# ---- 6. run selection -----------------------------------------------------
cat("\n-- run selection --\n")
sel <- scan$attempts %>% filter(!is_duplicate) %>% select_runs()
one_per <- sel %>% filter(is_selected) %>%
  count(pid, session_date, log_run) %>% pull(n)
check("exactly one run of record per participant/session/run",
      all(one_per == 1))
check("a triggered attempt always beats an untriggered one",
      sel %>% group_by(pid, session_date, log_run) %>%
        filter(any(has_trigger)) %>%
        summarise(ok = has_trigger[is_selected], .groups = "drop") %>%
        pull(ok) %>% all())

# ---- summary --------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", .pass, .fail))
if (.fail > 0) quit(status = 1)
