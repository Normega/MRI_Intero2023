# test_intero_events.R
# Validation for the breathing-run event extraction.
#   Rscript -e 'testthat::test_file("Repo/Tests/test_intero_events.R")'
#
# Ground truth is P130, whose two breathing runs were read off the raw log by
# hand. If the epoch boundaries or the rating-carrier row ever change, these
# hard-coded values are what will catch it.

library(testthat)
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_R    <- file.path(dir_root, "Repo", "R")
dir_logs <- file.path(dir_root, "Data", "Logs")

source(file.path(dir_R, "intero_ids.R"))
source(file.path(dir_R, "intero_logs.R"))
source(file.path(dir_R, "intero_events.R"))

d130 <- read_log_file(file.path(
  dir_logs, "Log.Intero2023.P130.SessionTest_2024-06-03_11h02.37.251.csv"))
r1 <- extract_breathing_events(dplyr::filter(d130, Run == "ChangeDetection1"))
ev <- r1$events
pc <- r1$pacer

test_that("a run yields the expected epoch inventory", {
  expect_equal(sum(ev$trial_type == "paced_breathing"), 7L)
  expect_equal(sum(ev$trial_type == "free_breathing"),  7L)
  expect_equal(sum(ev$trial_type == "rating"),         14L)
  expect_equal(sum(ev$trial_type == "valence_rating"),  1L)
  expect_equal(nrow(ev), 29L)
})

test_that("epochs are ordered, non-negative and non-overlapping", {
  expect_false(is.unsorted(ev$onset))
  expect_true(all(ev$duration >= 0))
  # every epoch ends no later than the next one begins
  e <- ev$onset + ev$duration
  expect_true(all(e[-nrow(ev)] <= ev$onset[-1] + 1e-6))
})

test_that("epoch durations are in the designed range", {
  paced  <- ev %>% filter(trial_type == "paced_breathing")
  free   <- ev %>% filter(trial_type == "free_breathing")
  rating <- ev %>% filter(trial_type == "rating")
  expect_true(all(paced$duration  > 10 & paced$duration  < 35))
  expect_true(all(free$duration   > 15 & free$duration   < 25))
  expect_true(all(rating$duration > 8  & rating$duration < 16))
})

test_that("the paced design matrix matches the raw log for P130 run 1", {
  paced <- ev %>% filter(trial_type == "paced_breathing") %>% arrange(epoch_index)
  expect_equal(paced$change_direction,
               c("Faster", "No Change", "Faster", "Slower",
                 "Slower", "Faster", "Slower"))
  expect_equal(paced$salience,
               c("low", NA, "high", "low", "low", "low", "high"))
  expect_equal(paced$change_amount,
               c(0.7, 1.0, 0.7, 1.3, 1.3, 0.7, 1.3))
  expect_equal(paced$response,
               c("No Change", "Slower", "Faster", "Slower",
                 "Slower", "Faster", "Slower"))
})

test_that("signal-detection classes are assigned correctly", {
  paced <- ev %>% filter(trial_type == "paced_breathing") %>% arrange(epoch_index)
  expect_equal(paced$sdt_class,
               c("miss", "false_alarm", "hit", "hit", "hit", "hit", "hit"))
  # and the classifier itself, including the case P130 happens not to contain
  expect_equal(sdt_class("Faster", "Slower"), "wrong_direction")
  expect_equal(sdt_class("Slower", "Slower"), "hit")
  expect_equal(sdt_class("No Change", "No Change"), "correct_rejection")
  expect_true(is.na(sdt_class("Faster", NA)))
})

test_that("free periods carry no paced attributes and no fake accuracy", {
  free <- ev %>% filter(trial_type == "free_breathing")
  # the log scores every free period 0 because TrialType 'Free' can never
  # equal a Faster/Slower/NoChange response -- that must not survive
  expect_true(all(is.na(free$accuracy)))
  expect_true(all(is.na(free$change_amount)))
  expect_true(all(is.na(free$salience)))
  expect_true(all(is.na(free$sdt_class)))
  # but the judgement and the modulators are real
  expect_true(all(!is.na(free$response)))
  expect_true(all(free$confidence >= 1 & free$confidence <= 7))
  expect_true(all(free$arousal    >= 1 & free$arousal    <= 7))
})

test_that("confidence and arousal vary within a run, valence does not", {
  cyc <- ev %>% filter(trial_type %in% c("paced_breathing", "free_breathing"))
  expect_gt(length(unique(cyc$confidence)), 1)
  expect_gt(length(unique(cyc$arousal)), 1)
  # valence is rated once per run, so it is a run-level covariate only
  expect_equal(length(unique(ev$valence_run)), 1L)
  expect_equal(unique(ev$valence_run), 6)
})

test_that("ratings come from the *_Complete row without an extra shift", {
  # P130 run 1, cycle 1: a Faster/low trial reported as No Change,
  # confidence 5, arousal 3. If shift_responses were applied on top, these
  # would land on cycle 2 instead.
  p1 <- ev %>% filter(trial_type == "paced_breathing", epoch_index == 1)
  expect_equal(p1$response, "No Change")
  expect_equal(p1$confidence, 5)
  expect_equal(p1$arousal, 3)
  expect_equal(p1$accuracy, 0)
})

test_that("the rating epochs mirror the breathing epoch they belong to", {
  rp <- ev %>% filter(trial_type == "rating", rating_of == "paced_breathing")
  pb <- ev %>% filter(trial_type == "paced_breathing")
  expect_equal(rp$response[order(rp$epoch_index)],
               pb$response[order(pb$epoch_index)])
  expect_equal(sum(ev$rating_of == "free_breathing", na.rm = TRUE), 7L)
})

test_that("the pacer trace excludes boundary rows", {
  expect_equal(length(unique(pc$epoch_index)), 7L)
  # 7 epochs x 9 real pacer events; the 7 tied boundary markers are dropped
  expect_equal(nrow(pc), 63L)
  # no zero-length half-cycles survive -- those were the artefact
  expect_true(all(pc$period_s > 0.5, na.rm = TRUE))
  expect_true(all(pc$direction %in% c(-1, 0, 1) | is.na(pc$direction)))
})

test_that("the pacer trace recovers the ramp, not a step", {
  # A Faster trial must end quicker than it began: for low salience the pacer
  # ramps toward ChangeAmount rather than stepping, which is the whole basis
  # of the awareness manipulation.
  f1 <- pc %>% filter(epoch_index == 1, !is.na(period_s))
  expect_equal(unique(f1$change_direction), "Faster")
  expect_equal(unique(f1$salience), "low")
  h <- ceiling(nrow(f1) / 2)
  expect_lt(mean(tail(f1$period_s, h)), mean(head(f1$period_s, h)))
  # and a Slower trial must do the opposite
  s4 <- pc %>% filter(epoch_index == 4, !is.na(period_s))
  expect_equal(unique(s4$change_direction), "Slower")
  h4 <- ceiling(nrow(s4) / 2)
  expect_gt(mean(tail(s4$period_s, h4)), mean(head(s4$period_s, h4)))
})

test_that("physio join slots exist and are empty", {
  expect_true(all(is.na(ev$adherence_verified)))
  expect_true(all(is.na(ev$actual_change_physio)))
  expect_true(all(ev$stimulus_change_known[ev$trial_type == "paced_breathing"]))
  expect_false(any(ev$stimulus_change_known[ev$trial_type == "free_breathing"]))
})

test_that("a v1 run is refused rather than mis-parsed", {
  d51 <- read_log_file(file.path(
    dir_logs, "Log.Intero2023.P51.SessionTest_2023-09-17_21h37.43.983.csv"))
  v1 <- extract_breathing_events(dplyr::filter(d51, Run == "ChangeDetection1"))
  expect_equal(nrow(v1$events), 0L)
})

# ---------------------------------------------------------------------------
# v1 (P51-P59) harmonisation
# ---------------------------------------------------------------------------

d51 <- read_log_file(file.path(
  dir_logs, "Log.Intero2023.P51.SessionTest_2023-09-17_21h37.43.983.csv"))
v1 <- extract_breathing_run(dplyr::filter(d51, Run == "ChangeDetection1"))

test_that("v1 is dispatched to the v1 reader and yields balanced cycles", {
  expect_gt(nrow(v1$events), 0)
  expect_equal(sum(v1$events$trial_type == "paced_breathing"), 5L)
  expect_equal(sum(v1$events$trial_type == "free_breathing"),  5L)
  expect_equal(sum(v1$events$trial_type == "rating"),         10L)
  # P51 ran 5 cycles; P053-P059 ran 7. Cycle count is not a version marker.
  expect_true(all(v1$events$duration >= 0))
  expect_false(is.unsorted(v1$events$onset))
})

test_that("v1 and v2 emit the SAME schema so the cohorts bind cleanly", {
  expect_setequal(names(v1$events), names(ev))
  expect_setequal(names(v1$pacer),  names(pc))
  combined <- dplyr::bind_rows(v1$events, ev)
  expect_equal(nrow(combined), nrow(v1$events) + nrow(ev))
})

test_that("v1 cycle 1 matches the raw log", {
  p1 <- v1$events %>% filter(trial_type == "paced_breathing", epoch_index == 1)
  expect_equal(p1$change_direction, "Slower")
  expect_equal(p1$salience, "low")
  expect_equal(p1$change_amount, 1.25)
  expect_equal(p1$response, "Slower")
  expect_equal(p1$confidence, 4)
  expect_equal(p1$arousal, 4)
  expect_equal(p1$sdt_class, "hit")
})

test_that("the v1/v2 valence difference is preserved, not averaged away", {
  # v1 rates valence after EVERY period; v2 rates it once per run
  expect_true(all(!is.na(
    v1$events$valence_epoch[v1$events$trial_type == "paced_breathing"])))
  expect_true(all(is.na(ev$valence_epoch)))
  expect_gt(length(unique(na.omit(v1$events$valence_epoch))), 1)
  # run-level value stays available for both, so they pool at run level
  expect_equal(length(unique(v1$events$valence_run)), 1L)
  expect_true(is.finite(unique(v1$events$valence_run)))
  expect_equal(unique(v1$events$n_sub_ratings), 4L)
  expect_equal(unique(ev$n_sub_ratings), 3L)
})

test_that("v1 free periods carry the same traps handling as v2", {
  free1 <- v1$events %>% filter(trial_type == "free_breathing")
  expect_true(all(is.na(free1$accuracy)))
  expect_true(all(is.na(free1$sdt_class)))
  expect_true(all(!is.na(free1$response)))
})

test_that("v1 pacer trace recovers the individualised rate and its ramp", {
  p1 <- v1$pacer %>% filter(epoch_index == 1, !is.na(period_s))
  expect_gt(nrow(p1), 5)
  # P51 was paced at 8 brpm -> 3.75 s half-cycle, ramping x1.25 for Slower
  expect_gt(median(p1$period_s), 3.0)
  h <- ceiling(nrow(p1) / 2)
  expect_gt(mean(tail(p1$period_s, h)), mean(head(p1$period_s, h)))
})
