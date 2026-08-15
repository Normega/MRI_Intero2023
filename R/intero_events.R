# intero_events.R
# Event extraction for the Intero2023 breathing (change-detection) runs.
#
# ---------------------------------------------------------------------------
# DESIGN, as established from Task/InteroPilot.4.51.py and verified against
# the logs. Read this before changing the epoch boundaries.
#
# One breathing run = 7 cycles. Each cycle is:
#
#   Fixation -> [ChangeTrial x ~10] -> BreathRating_Start -> Change Trial Complete
#            -> [Free]              -> BreathRating_Start -> Free Period Complete
#
# so it yields FOUR modelled epochs:
#   paced_breathing   the pacer is driving the breath, rate ramping to ChangeAmount
#   rating            composite rating window (see below)
#   free_breathing    unpaced
#   rating            composite rating window
#
# THE RATING WINDOW IS ITS OWN EPOCH. getratings() runs three scales back to
# back -- speedrating() -> Response+Accuracy, confrating() -> Confidence,
# arousalrating() -> Arousal -- inside a single window bracketed by
# BreathRating_Start and the *_Complete row. It is ~12 s, comparable to the
# 17-20 s breathing epochs, so folding it into the trial would confound
# breathing with decision and motor report. In v2 the three sub-ratings are NOT
# separately timestamped (Task lines 396 and 418 are commented out), so the
# window is modelled as one composite epoch and cannot be decomposed. v1
# (P51-P59) DOES log SpeedRating_Start / ConfidenceRating_Start /
# ArousalRating_Start separately -- richer, and handled elsewhere.
#
# WHY THE *_Complete ROW CARRIES THE RATINGS. Elsewhere in these logs a
# response is written one row late (see shift_responses in intero_logs.R). Here
# no shift is needed: the *_Complete row IS the row written after getratings()
# returns, so its Response/Accuracy/Confidence/Arousal are this cycle's values.
# Applying shift_responses on top would push them onto the NEXT cycle.
#
# TWO TRAPS, both structural rather than behavioural:
#   * Accuracy = (Response == TrialType) is computed by the task. Free periods
#     have TrialType 'Free' and a Faster/Slower/NoChange response, so they can
#     never match and every free period is scored 0. Forced to NA here.
#   * ChangeAmount and Salience carry forward from the preceding paced trial
#     into the free rows. Nulled here.
#
# VALENCE IS RUN-LEVEL. valencerating() is called once per run (Task line 836);
# the column then persists by carry-forward. It has no within-run variance and
# is NOT a within-run parametric modulator. Carried as valence_run for
# convenience. The per-cycle modulators are confidence and arousal.
#
# Author: RADlab

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble); library(purrr)
})

CHANGE_TYPES <- c("Faster", "Slower", "No Change")

# Two ChangeTrial rows closer together than this are the same pacer event: a
# boundary/carry-over marker plus the real row. See the dedup block below for
# why 0.25 s is the right order of magnitude.
PACER_TIE_TOL_S <- 0.25

# Task writes "No Change" with a space; be tolerant of spacing variants.
.norm_resp <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x <- sub("(?i)^no\\s*change$", "No Change", x, perl = TRUE)
  x
}

.blank_to_na <- function(x) {
  x <- trimws(as.character(x)); x[x == ""] <- NA_character_; x
}

.num <- function(x) suppressWarnings(as.numeric(.blank_to_na(x)))

#' Signal-detection class for one paced trial.
#'
#' `wrong_direction` is kept distinct from `hit` and `miss`: the participant
#' detected a change but misidentified its sign, which is neither a clean
#' detection nor a clean failure of awareness. Collapsing it into `miss` would
#' inflate the "change without awareness" cell with trials that were plainly
#' noticed.
sdt_class <- function(trial_type, response) {
  tt <- .norm_resp(trial_type); rp <- .norm_resp(response)
  dplyr::case_when(
    is.na(tt) | is.na(rp)                          ~ NA_character_,
    tt == "No Change" & rp == "No Change"          ~ "correct_rejection",
    tt == "No Change" & rp != "No Change"          ~ "false_alarm",
    tt != "No Change" & rp == tt                   ~ "hit",
    tt != "No Change" & rp == "No Change"          ~ "miss",
    tt != "No Change"                              ~ "wrong_direction",
    TRUE                                           ~ NA_character_
  )
}

# ---------------------------------------------------------------------------
# v1 (P51-P59). Same four epochs per cycle, but the log looks different:
#
#   * The pacer half-cycles are `Fixation` rows carrying TrialType/Salience/
#     ChangeAmount -- there are no ChangeTrial rows at all.
#   * Ratings are FOUR separately-timestamped events, not three collapsed into
#     one: SpeedRating_Start -> ConfidenceRating_Start -> ArousalRating_Start ->
#     ValenceRating_Start -> *_Complete. So v1 has finer rating timing than v2,
#     and it rates VALENCE ONCE PER PERIOD rather than once per run.
#   * Cycle count varies and is NOT a version marker: P51 ran 5 per run,
#     P053-P059 ran 7, same as v2. Never assume a fixed number.
#   * No `Free` marker row: the free period runs from the end of the paced
#     rating window to the first SpeedRating_Start of the free block.
#   * The paced pacer rate is individualised (P51 = 8 brpm) rather than fixed
#     at 12, so the pacer trace is the only record of the actual stimulus.
#
# Harmonisation: emit the SAME epochs and columns as v2, reading all four
# ratings off the *_Complete row, where they have all accumulated by
# carry-forward exactly as in v2. The genuine v1/v2 differences are preserved
# as data rather than erased -- valence_epoch is populated for v1 only, and
# n_sub_ratings records 4 vs 3. valence_run is the run mean for v1 so the two
# cohorts remain poolable at run level.
# ---------------------------------------------------------------------------

extract_breathing_events_v1 <- function(d) {
  empty <- list(events = tibble(), pacer = tibble())
  if (nrow(d) == 0) return(empty)

  desc <- trimws(as.character(d$Description))
  tt   <- d$run_time_s
  ttype <- .norm_resp(d$TrialType)

  i_spd <- which(desc == "SpeedRating_Start")
  i_ctc <- which(desc == "Change Trial Complete")
  i_fpc <- which(desc == "Free Period Complete")
  # pacer half-cycles: Fixation rows that carry a paced TrialType
  i_fix <- which(desc == "Fixation" & !is.na(ttype) & ttype %in% CHANGE_TYPES)

  if (!length(i_spd) || !length(i_ctc)) return(empty)

  ev <- list(); pc <- list()
  n_cycles <- length(i_ctc)
  prev_end_t <- 0; prev_end <- 0L

  for (k in seq_len(n_cycles)) {
    k_ctc <- i_ctc[k]

    spd_p <- i_spd[i_spd < k_ctc & i_spd > prev_end]
    if (!length(spd_p)) { prev_end <- k_ctc; next }
    spd_p <- min(spd_p)

    fix_k <- i_fix[i_fix > prev_end & i_fix < spd_p]
    if (!length(fix_k)) { prev_end <- k_ctc; next }
    i_last <- max(fix_k)

    ev[[length(ev) + 1]] <- tibble(
      onset            = tt[min(fix_k)],
      duration         = max(0, tt[spd_p] - tt[min(fix_k)]),
      trial_type       = "paced_breathing",
      epoch_index      = k,
      change_direction = ttype[i_last],
      change_amount    = .num(d$ChangeAmount[i_last]),
      salience         = .blank_to_na(d$Salience[i_last]),
      n_pacer_cycles   = length(fix_k),
      response   = .norm_resp(d$Response[k_ctc]),
      accuracy   = .num(d$Accuracy[k_ctc]),
      confidence = .num(d$Confidence[k_ctc]),
      arousal    = .num(d$Arousal[k_ctc]),
      valence_epoch = .num(d$Valence[k_ctc])
    )
    ev[[length(ev) + 1]] <- tibble(
      onset = tt[spd_p], duration = max(0, tt[k_ctc] - tt[spd_p]),
      trial_type = "rating", epoch_index = k, rating_of = "paced_breathing",
      response   = .norm_resp(d$Response[k_ctc]),
      accuracy   = .num(d$Accuracy[k_ctc]),
      confidence = .num(d$Confidence[k_ctc]),
      arousal    = .num(d$Arousal[k_ctc]),
      valence_epoch = .num(d$Valence[k_ctc])
    )

    on <- tt[fix_k]
    pc[[length(pc) + 1]] <- tibble(
      epoch_index = k, pulse_index = seq_along(fix_k), onset = on,
      period_s = c(diff(on), NA_real_),
      direction = NA_real_,          # v1 logs no Direction column value
      change_direction = ttype[i_last],
      salience = .blank_to_na(d$Salience[i_last])
    )

    # free period: from the end of the paced rating to the free block's first
    # SpeedRating_Start. There is no Free marker row in v1.
    fpc_k <- i_fpc[i_fpc > k_ctc]
    if (length(fpc_k)) {
      k_fpc <- min(fpc_k)
      spd_f <- i_spd[i_spd > k_ctc & i_spd < k_fpc]
      if (length(spd_f)) {
        spd_f <- min(spd_f)
        ev[[length(ev) + 1]] <- tibble(
          onset = tt[k_ctc], duration = max(0, tt[spd_f] - tt[k_ctc]),
          trial_type = "free_breathing", epoch_index = k,
          response   = .norm_resp(d$Response[k_fpc]),
          accuracy   = NA_real_,     # structurally 0 in v1 too; not real
          confidence = .num(d$Confidence[k_fpc]),
          arousal    = .num(d$Arousal[k_fpc]),
          valence_epoch = .num(d$Valence[k_fpc])
        )
        ev[[length(ev) + 1]] <- tibble(
          onset = tt[spd_f], duration = max(0, tt[k_fpc] - tt[spd_f]),
          trial_type = "rating", epoch_index = k, rating_of = "free_breathing",
          response   = .norm_resp(d$Response[k_fpc]),
          accuracy   = NA_real_,
          confidence = .num(d$Confidence[k_fpc]),
          arousal    = .num(d$Arousal[k_fpc]),
          valence_epoch = .num(d$Valence[k_fpc])
        )
      }
      prev_end <- k_fpc
    } else prev_end <- k_ctc
  }

  if (!length(ev)) return(empty)

  events <- bind_rows(ev) %>%
    mutate(
      # v1 rates valence per period; the run mean keeps it comparable with v2's
      # single run-level rating without pretending the two are the same thing.
      valence_run = mean(.data$valence_epoch[.data$trial_type %in%
                           c("paced_breathing", "free_breathing")], na.rm = TRUE),
      n_sub_ratings = 4L,
      sdt_class = ifelse(.data$trial_type == "paced_breathing",
                         sdt_class(.data$change_direction, .data$response),
                         NA_character_),
      stimulus_change_known = .data$trial_type == "paced_breathing",
      adherence_verified   = NA,
      actual_change_physio = NA_character_
    ) %>%
    arrange(.data$onset)

  list(events = events, pacer = if (length(pc)) bind_rows(pc) else tibble())
}

#' Dispatch on what the log actually contains, not on the participant number.
#'
#' Protocol version is inferred from structure so a mislabelled or transitional
#' session cannot be parsed with the wrong reader.
extract_breathing_run <- function(d) {
  if (nrow(d) == 0) return(list(events = tibble(), pacer = tibble()))
  desc <- trimws(as.character(d$Description))
  if (any(desc == "ChangeTrial"))          extract_breathing_events(d)
  else if (any(desc == "SpeedRating_Start")) extract_breathing_events_v1(d)
  else list(events = tibble(), pacer = tibble())
}

#' Extract modelled epochs + the pacer trace for ONE v2 breathing run.
#'
#' @param d rows of a single ChangeDetection run, in log order, as returned by
#'   read_log_file()
#' @return list(events = tibble, pacer = tibble)
extract_breathing_events <- function(d) {
  empty <- list(events = tibble(), pacer = tibble())
  if (nrow(d) == 0) return(empty)

  desc <- trimws(as.character(d$Description))
  tt   <- d$run_time_s

  i_brs <- which(desc == "BreathRating_Start")
  i_ctc <- which(desc == "Change Trial Complete")
  i_fpc <- which(desc == "Free Period Complete")
  i_fr  <- which(desc == "Free")
  i_ct  <- which(desc == "ChangeTrial")

  # v1 runs have no ChangeTrial rows at all; refuse rather than emit nonsense.
  if (!length(i_ct) || !length(i_ctc)) return(empty)

  ev <- list(); pc <- list()

  # ---- run-level valence ---------------------------------------------------
  i_vs <- which(desc == "ValenceRating_Start")
  i_vr <- which(desc == "ValenceRating")
  valence_run <- if (length(i_vr)) .num(d$Valence[i_vr[1]]) else NA_real_
  if (length(i_vs) && length(i_vr) && !is.na(tt[i_vs[1]])) {
    ev[[length(ev) + 1]] <- tibble(
      onset = tt[i_vs[1]],
      duration = max(0, tt[i_vr[1]] - tt[i_vs[1]]),
      trial_type = "valence_rating", epoch_index = NA_integer_
    )
  }

  # ---- cycles --------------------------------------------------------------
  n_cycles <- length(i_ctc)
  prev_end <- 0L

  for (k in seq_len(n_cycles)) {
    k_ctc <- i_ctc[k]

    # rating window for the paced trial: the BreathRating_Start just before it
    brs_p <- i_brs[i_brs < k_ctc & i_brs > prev_end]
    if (!length(brs_p)) { prev_end <- k_ctc; next }
    brs_p <- max(brs_p)

    # paced epoch: the ChangeTrial rows between the previous cycle and this
    # rating window. The first is a boundary row carrying the PREVIOUS trial's
    # stale TrialType, but it shares a timestamp with the true first row, so
    # taking the earliest index gives the right onset either way. Attributes
    # are read from the LAST ChangeTrial row, which is unambiguously this
    # trial's.
    ct_k <- i_ct[i_ct > prev_end & i_ct < brs_p]
    if (!length(ct_k)) { prev_end <- k_ctc; next }
    i_first <- min(ct_k); i_last <- max(ct_k)

    # Drop boundary rows. A ChangeTrial sharing its timestamp with the next
    # ChangeTrial is a carry-over marker holding the PREVIOUS trial's stale
    # TrialType, not a pacer event. Left in, it contributes a 0-second
    # half-cycle that both understates the mean pacer period and injects a
    # phantom transition into the trace the belt gets regressed against.
    # Keeping the later row of each tied pair also keeps the correct TrialType.
    # Tolerance, not equality: the boundary row usually shares the exact
    # timestamp but is sometimes written ~1 ms earlier (P130 run 1 cycle 4).
    # The threshold is set from the design rather than from jitter -- the
    # shortest possible real half-cycle is the fastest paced rate (16 brpm ->
    # 1.875 s half-period) scaled by the largest speed-up (x0.7), about 1.3 s.
    # 0.25 s therefore separates the two cases by a wide margin in both
    # directions.
    ct_p <- ct_k
    if (length(ct_p) > 1) {
      tied <- c(abs(diff(tt[ct_p])) < PACER_TIE_TOL_S, FALSE)
      tied[is.na(tied)] <- FALSE
      ct_p <- ct_p[!tied]
    }

    ev[[length(ev) + 1]] <- tibble(
      onset            = tt[i_first],
      duration         = max(0, tt[brs_p] - tt[i_first]),
      trial_type       = "paced_breathing",
      epoch_index      = k,
      change_direction = .norm_resp(d$TrialType[i_last]),
      change_amount    = .num(d$ChangeAmount[i_last]),
      salience         = .blank_to_na(d$Salience[i_last]),
      n_pacer_cycles   = length(ct_p)
    )

    # rating window (paced)
    ev[[length(ev) + 1]] <- tibble(
      onset       = tt[brs_p],
      duration    = max(0, tt[k_ctc] - tt[brs_p]),
      trial_type  = "rating",
      epoch_index = k,
      rating_of   = "paced_breathing",
      response    = .norm_resp(d$Response[k_ctc]),
      accuracy    = .num(d$Accuracy[k_ctc]),
      confidence  = .num(d$Confidence[k_ctc]),
      arousal     = .num(d$Arousal[k_ctc])
    )

    # attach the paced ratings back onto the breathing epoch too, so a single
    # events file supports both "model the breath" and "model the report"
    j <- length(ev) - 1L
    ev[[j]]$response   <- .norm_resp(d$Response[k_ctc])
    ev[[j]]$accuracy   <- .num(d$Accuracy[k_ctc])
    ev[[j]]$confidence <- .num(d$Confidence[k_ctc])
    ev[[j]]$arousal    <- .num(d$Arousal[k_ctc])

    # ---- pacer trace for this epoch ----------------------------------------
    on <- tt[ct_p]
    pc[[length(pc) + 1]] <- tibble(
      epoch_index = k,
      pulse_index = seq_along(ct_p),
      onset       = on,
      # Interval to the next pacer event: the instantaneous half-cycle period.
      # This is what encodes the ramp -- the change is gradual for low salience
      # and a step at the midpoint for high, so the trajectory (not just
      # ChangeAmount) is the stimulus ground truth for physio adherence.
      period_s    = c(diff(on), NA_real_),
      direction   = .num(d$Direction[ct_p]),
      change_direction = .norm_resp(d$TrialType[i_last]),
      salience    = .blank_to_na(d$Salience[i_last])
    )

    # ---- free period ---------------------------------------------------------
    fr_k <- i_fr[i_fr > k_ctc]
    fpc_k <- i_fpc[i_fpc > k_ctc]
    if (length(fr_k) && length(fpc_k)) {
      i_free <- min(fr_k); k_fpc <- min(fpc_k)
      brs_f <- i_brs[i_brs > i_free & i_brs < k_fpc]
      if (length(brs_f)) {
        brs_f <- min(brs_f)
        ev[[length(ev) + 1]] <- tibble(
          onset       = tt[i_free],
          duration    = max(0, tt[brs_f] - tt[i_free]),
          trial_type  = "free_breathing",
          epoch_index = k,
          # change_amount / salience deliberately absent: they carry forward
          # from the paced trial and mean nothing here.
          response    = .norm_resp(d$Response[k_fpc]),
          accuracy    = NA_real_,      # structurally 0 in the log; not real
          confidence  = .num(d$Confidence[k_fpc]),
          arousal     = .num(d$Arousal[k_fpc])
        )
        ev[[length(ev) + 1]] <- tibble(
          onset       = tt[brs_f],
          duration    = max(0, tt[k_fpc] - tt[brs_f]),
          trial_type  = "rating",
          epoch_index = k,
          rating_of   = "free_breathing",
          response    = .norm_resp(d$Response[k_fpc]),
          accuracy    = NA_real_,
          confidence  = .num(d$Confidence[k_fpc]),
          arousal     = .num(d$Arousal[k_fpc])
        )
      }
      prev_end <- k_fpc
    } else {
      prev_end <- k_ctc
    }
  }

  if (!length(ev)) return(empty)

  events <- bind_rows(ev) %>%
    mutate(
      valence_run = valence_run,
      # v2 rates valence once per run, so there is no per-epoch value; the
      # column exists only so v1 and v2 tables bind without coercion.
      valence_epoch = NA_real_,
      n_sub_ratings = 3L,
      sdt_class = ifelse(.data$trial_type == "paced_breathing",
                         sdt_class(.data$change_direction, .data$response),
                         NA_character_),
      # For paced trials the stimulus trajectory is fully recoverable from the
      # pacer trace; for free periods there is no referent at all until the
      # belt is analysed.
      stimulus_change_known = .data$trial_type == "paced_breathing",
      # Reserved for the physio pass. Left NA so results can be merged on
      # (subject, run, epoch_index) without regenerating events.
      adherence_verified   = NA,
      actual_change_physio = NA_character_
    ) %>%
    arrange(.data$onset)

  list(events = events, pacer = if (length(pc)) bind_rows(pc) else tibble())
}
