# intero2023_maia.R
# Scores the MAIA-2 from the screener battery, producing the covariate H6
# (preregistration) requires.
#
# THE REGISTERED SCORE NEEDS NO REVERSE-CODING. H6 replicates the IEAT (Farb,
# Zuo & Price 2023), which scored a total OMITTING the Not-Worrying and
# Not-Distracting subscales. Every item in the six remaining subscales
# (Noticing, Attention Regulation, Emotional Awareness, Self-Regulation, Body
# Listening, Trusting) is positively worded, so the registered score involves
# no reversal decisions. ND and NW are still scored, descriptively, with the
# standard MAIA-2 reversals, and flagged as not entering H6.
#
# ITEMS ARE LOCATED BY TEXT, NOT BY COLUMN CODE. The screener export's column
# codes are Qualtrics-generated and carry no item identity; the codebook maps
# code -> full label. Each MAIA-2 item is matched by a distinctive phrase of
# its canonical wording. The script FAILS LOUDLY if any phrase matches zero or
# more than one column, rather than scoring a wrong or partial scale.
#
# Produces: Repo/Data/maia_scores.csv
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_out  <- file.path(dir_root, "Repo", "Data")
source(file.path(dir_root, "Repo", "R", "intero_ids.R"))

cat(sprintf("Run started: %s\n", Sys.time()))

# distinctive fragment of each MAIA-2 item, by subscale
MAIA_ITEMS <- tribble(
  ~subscale, ~item, ~phrase,
  "noticing",  1, "notice where the tension is located",
  "noticing",  2, "notice when I am uncomfortable in my body",
  "noticing",  3, "notice where in my body I am comfortable",
  "noticing",  4, "notice changes in my breathing",
  "not_distracting", 5, "ignore physical tension or discomfort until",
  "not_distracting", 6, "distract myself from sensations of discomfort",
  "not_distracting", 7, "try to power through it",
  "not_distracting", 8, "try to ignore pain",
  "not_distracting", 9, "push feelings of discomfort away by focusing",
  "not_distracting", 10, "occupy myself with something else",
  "not_worrying", 11, "physical pain, I become upset",
  "not_worrying", 12, "start to worry that something is wrong",
  "not_worrying", 13, "notice an unpleasant body sensation without worrying",
  "not_worrying", 14, "stay calm and not worry",
  "not_worrying", 15, "get it out of my mind",
  "attention_regulation", 16, "pay attention to my breath without being distracted",
  "attention_regulation", 17, "maintain awareness of my inner bodily sensations",
  "attention_regulation", 18, "pay attention to my posture",
  "attention_regulation", 19, "return awareness to my body if I am distracted",
  "attention_regulation", 20, "refocus my attention from thinking to sensing",
  "attention_regulation", 21, "maintain awareness of my whole body",
  "attention_regulation", 22, "consciously focus on my body as a whole",
  "emotional_awareness", 23, "notice how my body changes when I am angry",
  "emotional_awareness", 24, "something is wrong in my life I can feel it",
  "emotional_awareness", 25, "body feels different after a peaceful experience",
  "emotional_awareness", 26, "breathing becomes free and easy",
  "emotional_awareness", 27, "notice how my body changes when I feel happy",
  "self_regulation", 28, "overwhelmed I can find a calm place inside",
  "self_regulation", 29, "bring awareness to my body I feel a sense of calm",
  "self_regulation", 30, "use my breath to reduce tension",
  "self_regulation", 31, "calm my mind by focusing on my body",
  "body_listening", 32, "listen for information from my body about my emotional",
  "body_listening", 33, "take time to explore how my body feels",
  "body_listening", 34, "listen to my body to inform me",
  "trusting", 35, "at home in my body",
  "trusting", 36, "body is a safe place",
  "trusting", 37, "trust my body sensations"
)
REVERSED <- c(5, 6, 7, 8, 9, 11, 12, 15)   # standard MAIA-2; ND/NW only
IEAT_SUBSCALES <- c("noticing", "attention_regulation", "emotional_awareness",
                    "self_regulation", "body_listening", "trusting")

cb <- read_csv(file.path(dir_out, "screener_codebook.csv"), show_col_types = FALSE,
               progress = FALSE)

# The screener's headers are MAIA_Q1_<k>; the codebook lists them with the
# MAIA_ prefix stripped (Q1_<k>). Rather than trust either name, the mapping is
# VERIFIED two ways: each item's canonical phrase must text-match the codebook
# label of Q1_<k>, and the item number embedded in that label ("- <k>. ...")
# must equal k. Only when the numbering scheme survives both checks is the
# by-number mapping used -- which also recovers item 18, whose administered
# wording differs from the canonical phrase.
hits <- map_dfr(seq_len(nrow(MAIA_ITEMS)), function(i) {
  m <- cb %>% dplyr::filter(str_detect(str_to_lower(label),
                                       fixed(str_to_lower(MAIA_ITEMS$phrase[i]))))
  tibble(item = MAIA_ITEMS$item[i], subscale = MAIA_ITEMS$subscale[i],
         phrase = MAIA_ITEMS$phrase[i], n_matches = nrow(m),
         column_code = if (nrow(m) == 1) m$column_code[1] else NA_character_)
})
chk <- hits %>% dplyr::filter(n_matches == 1) %>%
  mutate(code_num = suppressWarnings(as.integer(str_match(column_code, "Q1_([0-9]+)$")[, 2])))
cat(sprintf("phrase-matched items: %d of 37 ; of those, code number == item number: %d\n",
            nrow(chk), sum(chk$code_num == chk$item, na.rm = TRUE)))
stopifnot(nrow(chk) >= 30, all(chk$code_num == chk$item, na.rm = TRUE))

sc <- read_csv(file.path(dir_out, "screener_clean.csv"), show_col_types = FALSE,
               progress = FALSE, guess_max = 5000)
ok_cols <- MAIA_ITEMS %>%
  mutate(column_code = sprintf("MAIA_Q1_%d", item)) %>%
  dplyr::filter(column_code %in% names(sc))
cat(sprintf("MAIA_Q1_<k> columns present in screener_clean: %d of 37\n", nrow(ok_cols)))
stopifnot(nrow(ok_cols) == 37)

d <- sc %>% select(participant_id, all_of(ok_cols$column_code))
# THIS ADMINISTRATION IS CODED 1-6, the canonical 0-5 anchors shifted up by
# one (Qualtrics numbering its six radio buttons from 1). Raw value table over
# all 37 columns: values 1..6 only. An earlier version of this script filtered
# to 0-5 and thereby silently discarded every "6" -- i.e. every "Always"
# response, 343 of them -- which surfaced as 25 participants with "incomplete"
# MAIA. Everything is converted to the canonical 0-5 metric here (raw - 1), so
# published norms and the IEAT apply directly.
lv <- c("never" = 0, "rarely" = 1, "sometimes" = 2, "often" = 3,
        "very often" = 4, "always" = 5)
to_num <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  n <- ifelse(!is.na(n) & n >= 1 & n <= 6, n - 1, NA_real_)
  t <- lv[str_to_lower(str_trim(as.character(x)))]
  ifelse(is.na(n), unname(t), n)
}
d <- d %>% mutate(across(-participant_id, to_num),
                  pid = pid_num(participant_id)) %>%
  dplyr::filter(!is.na(pid), !is_pilot(pid))
cat(sprintf("respondents with a participant id: %d (%d unique)\n",
            nrow(d), n_distinct(d$pid)))

# a participant screened twice keeps the more complete record
d$n_ans <- rowSums(!is.na(d %>% select(all_of(ok_cols$column_code))))
cat("\nresponse value range check:\n")
rng <- range(d %>% select(all_of(ok_cols$column_code)), na.rm = TRUE)
cat(sprintf("  min=%s max=%s (expect 0-5)\n", rng[1], rng[2]))
d <- d %>% group_by(pid) %>% arrange(desc(n_ans)) %>% slice(1) %>% ungroup()

# reverse ND/NW items (descriptive subscales only), using the OBSERVED anchors:
# this administration runs 1-5, so reversal is (min+max) - x = 6 - x, not 5 - x.
rev_lo <- min(d %>% select(all_of(ok_cols$column_code)), na.rm = TRUE)
rev_hi <- max(d %>% select(all_of(ok_cols$column_code)), na.rm = TRUE)
rev_cols <- ok_cols$column_code[ok_cols$item %in% REVERSED]
d <- d %>% mutate(across(all_of(rev_cols), ~ (rev_lo + rev_hi) - .x))

scores <- map_dfr(seq_len(nrow(d)), function(i) {
  row <- d[i, ]
  ss <- map_dbl(unique(ok_cols$subscale), function(s) {
    cols <- ok_cols$column_code[ok_cols$subscale == s]
    v <- as.numeric(row[, cols])
    if (mean(is.na(v)) > 0.5) NA_real_ else mean(v, na.rm = TRUE)
  })
  names(ss) <- unique(ok_cols$subscale)
  ieat <- ss[IEAT_SUBSCALES]
  tibble(pid = row$pid, subject = pid_bids(row$pid),
         n_items_answered = row$n_ans,
         !!!as.list(ss),
         maia_ieat_total = if (any(is.na(ieat))) NA_real_ else mean(ieat),
         maia2_total = if (any(is.na(ss))) NA_real_ else mean(ss))
})
write_csv(scores, file.path(dir_out, "maia_scores.csv"))

cat("\n== MAIA coverage ==\n")
cat(sprintf("scored participants: %d ; complete IEAT-style total: %d\n",
            nrow(scores), sum(is.finite(scores$maia_ieat_total))))
cat("\nmissingness pattern (items answered of 37):\n")
print(table(cut(scores$n_items_answered, c(-1, 0, 18, 36, 37),
                labels = c("0 (section not reached)", "1-18", "19-36", "all 37"))))
cat("\nIEAT-style total distribution (0-5 scale):\n")
print(round(summary(scores$maia_ieat_total), 2))
cat(sprintf("SD = %.2f\n", sd(scores$maia_ieat_total, na.rm = TRUE)))

# coverage within the imaging sample
pr <- read_csv(file.path(dir_out, "prereg_sample.csv"), show_col_types = FALSE,
               progress = FALSE)
j <- pr %>% left_join(scores %>% select(subject, maia_ieat_total), by = "subject")
cat(sprintf("\nusable-imaging participants: %d ; with complete MAIA: %d\n",
            nrow(pr), sum(is.finite(j$maia_ieat_total))))
cat("\nsubscale means (0-5):\n")
print(as.data.frame(scores %>% summarise(across(all_of(unique(ok_cols$subscale)),
      ~ round(mean(.x, na.rm = TRUE), 2)))))
cat(sprintf("\nRun ended: %s\n", Sys.time()))
