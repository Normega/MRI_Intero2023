# intero2023_clean.R
# De-identifies and scores Intero2023 questionnaire data
# Produces: screener_clean.csv, postscan_clean.csv, run_notes.csv
# Does NOT filter responses — all rows retained for participant flow reporting
#
# Script version: 1.0
# Author: RADlab

# Set Up ---------
## Load libraries ---------
packages <- c("readxl", "readr", "tidyverse", "lubridate")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Paths ---------
dir_source <- "I:/Shared drives/MRI_Intero2023/OSF/Questionnaires/"
dir_out    <- "I:/Shared drives/MRI_Intero2023/Repo/Data/"

f_screener <- file.path(dir_source, "Screener_Intero2023.xlsx")
f_postscan <- file.path(dir_source, "Post_Scan_Intero2023.xlsx")
f_tracker  <- file.path(dir_source, "Intero2023_Participant_Tracker.xlsx")
f_notes    <- file.path(dir_source, "Interoception_2023_Session_notes.docx")

dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

# ── 1. Participant ID Lookup ────────────────────────────────────────────────

tracker_raw <- readxl::read_excel(f_tracker, sheet = "Sheet1")

# Build email → P### map (keep only rows with both email and P###)
id_lookup <- tracker_raw |>
  dplyr::rename(participant_id = `Participant #`, email = Email) |>
  dplyr::filter(!is.na(participant_id), !is.na(email)) |>
  dplyr::mutate(email = trimws(tolower(email))) |>
  dplyr::select(participant_id, email) |>
  dplyr::distinct(email, .keep_all = TRUE)

cat(sprintf("ID lookup built: %d participants\n", nrow(id_lookup)))

# Email-scrubbing safety net: replaces any email address in any column that
# contains text with "[email removed]" before a data frame is written to CSV.
# Coerces to character first so glue/factor/other classes are still scrubbed.
strip_emails <- function(df) {
  email_re <- "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}"
  scrub_one <- function(x) {
    if (is.numeric(x) || is.logical(x)) return(x)   # leave numeric/logical alone
    gsub(email_re, "[email removed]", as.character(x), perl = TRUE)
  }
  df[] <- lapply(df, scrub_one)
  df
}

# ── 2. Screener / Baseline (Qualtrics double-header) ───────────────────────

## 2a. Read headers robustly --------------------------------------------------
# Read first 2 rows forcing ALL columns to text so text headers in
# numeric-typed columns are not coerced to NA by readxl type guessing.
header_rows <- readxl::read_excel(f_screener, sheet = "Sheet0",
                                  col_names = FALSE, n_max = 2,
                                  col_types = "text")

header_codes  <- as.character(unlist(header_rows[1, ], use.names = FALSE))
header_labels <- as.character(unlist(header_rows[2, ], use.names = FALSE))

# Replace NA/empty codes with placeholders
na_idx <- which(is.na(header_codes) | header_codes == "NA" | header_codes == "")
if (length(na_idx) > 0) header_codes[na_idx] <- paste0("col_", na_idx)

# Save codebook
codebook <- tibble::tibble(column_code = header_codes, label = header_labels)
readr::write_csv(codebook, file.path(dir_out, "screener_codebook.csv"))

## 2b. Compute final column names (block-prefix BEFORE dedup) -----------------
# rename_by_label defined below in 2c; we need it first, so source order matters.
# Define it here so it is available before the data read.
rename_by_label <- function(codes, labels) {
  new_codes <- codes
  block_state <- "meta"
  for (i in seq_along(codes)) {
    lbl <- if (is.na(labels[i])) "" else as.character(labels[i])
    cd  <- if (is.na(codes[i]))  "" else as.character(codes[i])
    
    if (grepl("How often have you been upset", lbl, fixed = TRUE))        block_state <- "PSS"
    if (grepl("Feeling nervous, anxious, or on edge", lbl, fixed = TRUE)) block_state <- "PHQ4"
    if (grepl("When I am tense I notice", lbl, fixed = TRUE))             block_state <- "MAIA"
    if (grepl("I am able to watch myself thinking", lbl, fixed = TRUE))   block_state <- "MPoD"
    if (grepl("is reserved", lbl, fixed = TRUE))                          block_state <- "BFI"
    if (grepl("In the past few days.*emotions", lbl))                     block_state <- "Emotions"
    if (grepl("What does the term.*Wellbeing.*mean", lbl))               block_state <- "WB"
    if (grepl("In most ways my life", lbl, fixed = TRUE))                 block_state <- "SWLS"
    if (grepl("I felt that my choices were based", lbl, fixed = TRUE))    block_state <- "SDT"
    if (grepl("Overall, how satisfied are you with life", lbl, fixed = TRUE)) block_state <- "WBmulti"
    if (grepl("I lead a purposeful", lbl, fixed = TRUE))                  block_state <- "PWB8"
    if (grepl("Please describe what meditation", lbl, fixed = TRUE))      block_state <- "Meditation"
    if (grepl("We are interested in your day-to-day", lbl, fixed = TRUE)) block_state <- "MPoD"
    
    if (grepl("^Q1_", cd) || (cd == "Q1" && block_state %in% c("WB", "PSS", "PHQ4", "MAIA", "MPoD", "BFI"))) {
      new_codes[i] <- paste0(block_state, "_", cd)
    }
    if (block_state == "WBmulti" && grepl("^Q[5-9]_|^Domain 6", cd)) {
      new_codes[i] <- paste0("WBmulti_", cd)
    }
  }
  new_codes
}

# Block-prefix the ORIGINAL codes (makes scale blocks unique), then dedup as
# a final safety net for any genuine remaining collisions.
final_names <- rename_by_label(header_codes, header_labels)
final_names <- make.unique(final_names, sep = "_dup")

## 2b.1 Read data rows with final names directly -----------------------------
screener_raw <- readxl::read_excel(f_screener, sheet = "Sheet0",
                                   col_names = final_names, skip = 2)
screener_raw <- readr::type_convert(screener_raw, guess_integer = FALSE)

cat(sprintf("Screener raw: %d rows x %d cols\n", nrow(screener_raw), ncol(screener_raw)))

## 2b.5 Join P### (Email/Q99 names are unmodified by rename) ------------------
screener_raw <- screener_raw |>
  dplyr::mutate(
    email_lower = dplyr::coalesce(
      trimws(tolower(Email)),
      trimws(tolower(Q99))
    )
  ) |>
  dplyr::left_join(id_lookup, by = c("email_lower" = "email")) |>
  dplyr::select(-email_lower)

cat(sprintf("Screener P### matched: %d of %d rows\n",
            sum(!is.na(screener_raw$participant_id)), nrow(screener_raw)))

## 2d. Drop PII columns -------------------------------------------------------
pii_cols <- c("Name", "Email", "Q99",             # name/email fields
              "RecipientLastName", "RecipientFirstName", "RecipientEmail",
              "IPAddress", "LocationLatitude", "LocationLongitude",
              "ExternalReference", "DistributionChannel", "UserLanguage",
              "Q76_Id", "Q76_Name", "Q76_Size", "Q76_Type",  # consent signature
              "Emails", "drawdata")                # appended cols

screener_clean <- screener_raw |>
  dplyr::select(-dplyr::any_of(pii_cols)) |>
  dplyr::relocate(participant_id, .before = 1)

## 2c. (rename_by_label defined above) ----------------------------------------

## 2e. Flag attention/catch checks --------------------------------------------
# PSS attention check: PSS_Q1_8 should be the "Always 7" endpoint.
# Qualtrics stored it as text "Always7"; strip non-digits and compare to 7.
# CatchQ correct answer: "Kind"
screener_clean <- screener_clean |>
  dplyr::mutate(
    .pss_attn_num = suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(PSS_Q1_8)))),
    flag_pss_attention = dplyr::case_when(
      is.na(.pss_attn_num)   ~ NA,
      .pss_attn_num == 7     ~ FALSE,
      TRUE                   ~ TRUE
    ),
    flag_catch_q = dplyr::case_when(
      is.na(CatchQ)          ~ NA,
      CatchQ == "Kind"       ~ FALSE,
      TRUE                   ~ TRUE
    )
  ) |>
  dplyr::select(-.pss_attn_num)

## 2f. Score validated scales -------------------------------------------------

# Qualtrics exported some scales with the choice LABEL embedded plus its value
# (e.g. "Never1", "Always7", "Agree strongly5", "Neither agree nor disagree3").
# These three blocks are affected: PSS, MPoD, BFI. Extract the trailing numeric
# value so they score correctly. Scales that stored bare numbers are untouched.
recode_anchored <- function(x) {
  x <- as.character(x)
  # Pull the trailing run of digits ("Agree strongly5" -> "5", "7" -> "7")
  num <- sub("^.*?([0-9]+)\\s*$", "\\1", x)
  num[!grepl("[0-9]", x)] <- NA            # no digit present -> NA
  suppressWarnings(as.numeric(num))
}

anchored_cols <- c(paste0("PSS_Q1_", 1:11),
                   paste0("MPoD_Q1_", 1:15),
                   paste0("BFI_Q1_", 1:10))

screener_clean <- screener_clean |>
  dplyr::mutate(dplyr::across(dplyr::any_of(anchored_cols), recode_anchored))

### Helper: rowwise mean with NA handling
row_mean <- function(df, cols, min_valid = NULL) {
  mat <- as.matrix(dplyr::select(df, dplyr::all_of(cols)))
  mat <- apply(mat, 2, as.numeric)
  if (!is.null(min_valid)) {
    n_valid <- rowSums(!is.na(mat))
    result  <- rowMeans(mat, na.rm = TRUE)
    result[n_valid < min_valid] <- NA
    result
  } else {
    rowMeans(mat, na.rm = TRUE)
  }
}

row_sum <- function(df, cols, min_valid = NULL) {
  mat <- as.matrix(dplyr::select(df, dplyr::all_of(cols)))
  mat <- apply(mat, 2, as.numeric)
  if (!is.null(min_valid)) {
    n_valid <- rowSums(!is.na(mat))
    result  <- rowSums(mat, na.rm = TRUE)
    result[n_valid < min_valid] <- NA
    result
  } else {
    rowSums(mat, na.rm = TRUE)
  }
}

### PSS-10 (reverse items: 4,5,7,8 → mapped to PSS_Q1_4, Q1_5, Q1_7, Q1_9)
# Scale: 1-7 (survey used 1=Never to 7=Always)
# Reverse-scored items: 4 (confident), 5 (going your way), 7 (control irritations), 8 (on top of things)
# Note: Q1_8 is the attention check; actual PSS item 8 is Q1_9 due to inserted attention check item
pss_all      <- paste0("PSS_Q1_", 1:11)
pss_forward  <- paste0("PSS_Q1_", c(1, 2, 3, 6, 10, 11))
pss_reverse  <- paste0("PSS_Q1_", c(4, 5, 7, 9))   # Q1_9 = PSS item 8

# (PSS columns already numeric via recode_anchored in 2f preamble)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(pss_reverse),
                  ~ 8 - .x,
                  .names = "{.col}_r"),
    PSS_score = row_sum(dplyr::cur_data(),
                        c(pss_forward, paste0(pss_reverse, "_r")),
                        min_valid = 9)
  )

### PHQ-4 (0-3 each; GAD-2 = items 1-2, PHQ-2 = items 3-4)
phq_cols <- paste0("PHQ4_Q1_", 1:4)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(phq_cols), as.numeric),
    PHQ4_GAD2  = row_sum(dplyr::cur_data(), phq_cols[1:2], min_valid = 2),
    PHQ4_PHQ2  = row_sum(dplyr::cur_data(), phq_cols[3:4], min_valid = 2),
    PHQ4_total = row_sum(dplyr::cur_data(), phq_cols,      min_valid = 3)
  )

### MAIA-2 subscales (items 1-37, scale 0-5 recoded from 1-6 in survey)
# Subscales (1-indexed):
#   Noticing:           1-4
#   Not-Distracting:    5-10  (reverse: 5,6,7,8,9,10 — all reverse)
#   Not-Worrying:       11-15 (reverse: 11,12,15)
#   Attention Reg:      16-22
#   Emotional Awareness:23-27
#   Self-Regulation:    28-31
#   Body Listening:     32-34
#   Trusting:           35-37

maia_cols <- paste0("MAIA_Q1_", 1:37)

# MAIA survey used 1-6; published scale is 0-5; recode
screener_clean <- screener_clean |>
  dplyr::mutate(dplyr::across(dplyr::all_of(maia_cols),
                              ~ as.numeric(.x) - 1,
                              .names = "{.col}_r0"))

maia_r0 <- paste0(maia_cols, "_r0")

# Reverse items within 0-5 range: reversed = 5 - x
maia_rev_idx <- c(5:10, 11, 12, 15)   # item numbers to reverse
maia_rev_cols <- paste0("MAIA_Q1_", maia_rev_idx, "_r0")

screener_clean <- screener_clean |>
  dplyr::mutate(dplyr::across(dplyr::all_of(maia_rev_cols),
                              ~ 5 - .x,
                              .names = "{.col}_rev"))

# Build final item columns for scoring (reversed where needed)
maia_scored <- character(37)
for (i in 1:37) {
  base <- paste0("MAIA_Q1_", i, "_r0")
  rev  <- paste0(base, "_rev")
  maia_scored[i] <- if (i %in% maia_rev_idx && rev %in% names(screener_clean)) rev else base
}

subscale_mean <- function(df, item_indices, scored_cols, min_valid = NULL) {
  cols <- scored_cols[item_indices]
  row_mean(df, cols, min_valid = ceiling(length(cols) * 0.8))
}

screener_clean <- screener_clean |>
  dplyr::mutate(
    MAIA_noticing     = subscale_mean(dplyr::cur_data(), 1:4,   maia_scored),
    MAIA_notdist      = subscale_mean(dplyr::cur_data(), 5:10,  maia_scored),
    MAIA_notworry     = subscale_mean(dplyr::cur_data(), 11:15, maia_scored),
    MAIA_attnreg      = subscale_mean(dplyr::cur_data(), 16:22, maia_scored),
    MAIA_emotaware    = subscale_mean(dplyr::cur_data(), 23:27, maia_scored),
    MAIA_selfreg      = subscale_mean(dplyr::cur_data(), 28:31, maia_scored),
    MAIA_bodylisten   = subscale_mean(dplyr::cur_data(), 32:34, maia_scored),
    MAIA_trusting     = subscale_mean(dplyr::cur_data(), 35:37, maia_scored)
  )

### MPoD-t Decentering (3 subscales, items 1-15, scale 1-7)
#   Observer:     1-5
#   Transcendent: 6-10
#   Non-reactive: 11-15
mpod_cols <- paste0("MPoD_Q1_", 1:15)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(mpod_cols), as.numeric),
    MPoD_observer     = row_mean(dplyr::cur_data(), mpod_cols[1:5],   min_valid = 4),
    MPoD_transcendent = row_mean(dplyr::cur_data(), mpod_cols[6:10],  min_valid = 4),
    MPoD_nonreactive  = row_mean(dplyr::cur_data(), mpod_cols[11:15], min_valid = 4),
    MPoD_total        = row_mean(dplyr::cur_data(), mpod_cols,        min_valid = 12)
  )

### Big Five-10 (scale 1-5; reverse items 1,3,4,5,7)
bfi_cols <- paste0("BFI_Q1_", 1:10)
bfi_rev  <- paste0("BFI_Q1_", c(1, 3, 4, 5, 7))

screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(bfi_cols), as.numeric),
    dplyr::across(dplyr::all_of(bfi_rev), ~ 6 - .x, .names = "{.col}_r"),
    BFI_extraversion    = row_mean(dplyr::cur_data(),
                                   c("BFI_Q1_1_r", "BFI_Q1_6"),    min_valid = 2),
    BFI_agreeableness   = row_mean(dplyr::cur_data(),
                                   c("BFI_Q1_2",   "BFI_Q1_7_r"),  min_valid = 2),
    BFI_conscientiousness = row_mean(dplyr::cur_data(),
                                     c("BFI_Q1_3_r","BFI_Q1_8"),   min_valid = 2),
    BFI_neuroticism     = row_mean(dplyr::cur_data(),
                                   c("BFI_Q1_4_r", "BFI_Q1_9"),    min_valid = 2),
    BFI_openness        = row_mean(dplyr::cur_data(),
                                   c("BFI_Q1_5_r", "BFI_Q1_10"),   min_valid = 2)
  )

### SWLS (3-item version used here; sum, scale 1-7)
swls_cols <- c("SWLS_1", "SWLS_2", "SWLS_3")
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(swls_cols), as.numeric),
    SWLS_total = row_sum(dplyr::cur_data(), swls_cols, min_valid = 3)
  )

### SDT need satisfaction (3 subscales x 3 items, scale 1-5)
#   Autonomy:    SDT_1, SDT_2, SDT_3
#   Competence:  SDT_4, SDT_5, SDT_6
#   Relatedness: SDT_7, SDT_8, SDT_9
sdt_cols <- paste0("SDT_", 1:9)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(sdt_cols), as.numeric),
    SDT_autonomy    = row_mean(dplyr::cur_data(), sdt_cols[1:3], min_valid = 2),
    SDT_competence  = row_mean(dplyr::cur_data(), sdt_cols[4:6], min_valid = 2),
    SDT_relatedness = row_mean(dplyr::cur_data(), sdt_cols[7:9], min_valid = 2)
  )

### Wellbeing multidomain (domain means; 0-10 scales)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::starts_with("WBmulti_"), as.numeric),
    WBmulti_d1_happiness   = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Q5_1", "WBmulti_Q5_2"),        min_valid = 2),
    WBmulti_d2_health      = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Q6_1", "WBmulti_Q6_2"),        min_valid = 2),
    WBmulti_d3_meaning     = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Q7_1", "WBmulti_Q7_2"),        min_valid = 2),
    WBmulti_d4_character   = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Q8_1", "WBmulti_Q8_2"),        min_valid = 2),
    WBmulti_d5_social      = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Q9_1", "WBmulti_Q9_2"),        min_valid = 2),
    WBmulti_d6_financial   = row_mean(dplyr::cur_data(),
                                      c("WBmulti_Domain 6: Financial_1",
                                        "WBmulti_Domain 6: Financial_2"),        min_valid = 2)
  )

### PWB-8 (total mean, scale 1-7)
pwb_cols <- paste0("PWB8_", 1:8)
screener_clean <- screener_clean |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(pwb_cols), as.numeric),
    PWB8_total = row_mean(dplyr::cur_data(), pwb_cols, min_valid = 6)
  )

## 2g. Drop intermediate scoring columns (keep raw items + composites only) --
drop_intermediate <- c(
  paste0(pss_reverse, "_r"),
  maia_r0,
  maia_rev_cols,
  paste0(maia_rev_cols, "_rev"),
  bfi_rev |> (\(x) paste0(x, "_r"))()
)
screener_clean <- screener_clean |>
  dplyr::select(-dplyr::any_of(drop_intermediate))

## 2h. Add logging metadata ---------------------------------------------------
screener_clean <- screener_clean |>
  dplyr::mutate(
    log_script_version = "1.0",
    log_date_cleaned   = as.character(Sys.Date()),
    log_n_total_rows   = nrow(screener_clean)
  )

screener_clean <- strip_emails(screener_clean)
readr::write_csv(screener_clean, file.path(dir_out, "screener_clean.csv"))
cat(sprintf("screener_clean.csv written: %d rows x %d cols\n",
            nrow(screener_clean), ncol(screener_clean)))

# ── 3. Post-Scan Survey (Google Form, single header) ────────────────────────

postscan_raw <- readxl::read_excel(f_postscan, sheet = 1, col_names = TRUE)

cat(sprintf("Post-scan raw: %d rows x %d cols\n", nrow(postscan_raw), ncol(postscan_raw)))

# Post-scan form already contains the P### in "What is the subject ID?" — no
# email join needed. Detect that column robustly (matches "subject ID").
id_col_ps <- names(postscan_raw)[
  grepl("subject\\s*id", names(postscan_raw), ignore.case = TRUE)
][1]

if (is.na(id_col_ps)) {
  stop("Could not detect subject ID column in post-scan file. Check column names.")
} else {
  cat(sprintf("Post-scan ID column detected: '%s'\n", id_col_ps))
}

postscan_clean <- postscan_raw |>
  dplyr::rename(participant_id = dplyr::all_of(id_col_ps)) |>
  dplyr::mutate(participant_id = trimws(toupper(participant_id))) |>
  dplyr::relocate(participant_id, .before = 1) |>
  dplyr::mutate(
    log_script_version = "1.0",
    log_date_cleaned   = as.character(Sys.Date())
  )

cat(sprintf("Post-scan: %d of %d rows have a P### ID\n",
            sum(!is.na(postscan_clean$participant_id) &
                  postscan_clean$participant_id != "NA"),
            nrow(postscan_clean)))

postscan_clean <- strip_emails(postscan_clean)
readr::write_csv(postscan_clean, file.path(dir_out, "postscan_clean.csv"))
cat(sprintf("postscan_clean.csv written: %d rows x %d cols\n",
            nrow(postscan_clean), ncol(postscan_clean)))

# ── 4. Run Notes ─────────────────────────────────────────────────────────────

## 4a. RunCorrection tab -------------------------------------------------------
run_corr_raw <- readxl::read_excel(f_tracker, sheet = "RunCorrection",
                                   col_names = FALSE)

# Structure: col1 = P###, col2 = flag/count, col3 = run note (names the run),
#            col4 = "Text" general notes, col5 = "Scanning Notes" general notes
run_corr <- run_corr_raw[-1, ] |>                       # drop header row
  dplyr::rename(participant_id = 1, run_flag = 2,
                run_note = 3, text_note = 4, scanning_note = 5) |>
  dplyr::filter(!is.na(participant_id)) |>
  dplyr::mutate(
    participant_id = as.character(participant_id),
    run_note       = as.character(run_note),
    text_note      = as.character(text_note),
    scanning_note  = as.character(scanning_note)
  )

# Derive run column from the note text. Matches "loc 1/2", "breath 1/2", "auth".
# Returns the target column name, or NA if no clean single match.
derive_run_col <- function(note) {
  if (is.na(note)) return(NA_character_)
  n <- tolower(note)
  hits <- c(
    run1_loc1_note    = grepl("loc(alizer)?\\s*1", n),
    run2_loc2_note    = grepl("loc(alizer)?\\s*2", n),
    run3_breath1_note = grepl("breath\\s*1", n),
    run4_breath2_note = grepl("breath\\s*2", n),
    run5_auth1_note   = grepl("auth", n)
  )
  matched <- names(hits)[hits]
  if (length(matched) == 1) matched else NA_character_
}

run_cols <- c("run1_loc1_note", "run2_loc2_note", "run3_breath1_note",
              "run4_breath2_note", "run5_auth1_note")

# Build one row per participant with run columns + general notes + review flag
run_corr_wide <- run_corr |>
  dplyr::rowwise() |>
  dplyr::mutate(target_col = derive_run_col(run_note)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # A note exists but couldn't be mapped to a single run = needs review
    flag_review = !is.na(run_note) & is.na(target_col),
    # Combine col4 + col5 into a single general notes field
    general_notes_tracker = dplyr::case_when(
      !is.na(text_note) & !is.na(scanning_note) ~ paste(text_note, scanning_note, sep = " | "),
      !is.na(text_note)                          ~ text_note,
      !is.na(scanning_note)                      ~ scanning_note,
      TRUE                                       ~ NA_character_
    ),
    # Unmapped run notes get appended to general notes so nothing is lost
    general_notes_tracker = dplyr::case_when(
      flag_review & !is.na(general_notes_tracker) ~ paste0(general_notes_tracker, " | [unmapped run note: ", run_note, "]"),
      flag_review                                  ~ paste0("[unmapped run note: ", run_note, "]"),
      TRUE                                         ~ general_notes_tracker
    )
  )

# Place each mapped run_note into its derived run column
for (rc in run_cols) {
  run_corr_wide[[rc]] <- ifelse(run_corr_wide$target_col == rc & !is.na(run_corr_wide$target_col),
                                run_corr_wide$run_note, NA_character_)
}

run_corr_wide <- run_corr_wide |>
  dplyr::select(participant_id, dplyr::all_of(run_cols),
                general_notes_tracker, flag_review)

## 4b. Session notes docx ------------------------------------------------------
parse_session_notes <- function(path) {
  doc <- officer::read_docx(path)
  paras <- officer::docx_summary(doc) |>
    dplyr::filter(content_type == "paragraph") |>
    dplyr::pull(text)
  paras <- paras[nzchar(trimws(paras))]
  
  entries <- list()
  current_id   <- NA_character_
  current_date <- NA_character_
  current_notes <- character(0)
  in_notes <- FALSE
  
  flush_entry <- function() {
    if (!is.na(current_id)) {
      entries[[length(entries) + 1]] <<- tibble::tibble(
        participant_id      = current_id,
        session_date        = current_date,
        general_notes_docx  = paste(current_notes, collapse = " | ")
      )
    }
  }
  
  for (p in paras) {
    p <- trimws(p)
    # Detect participant ID line
    if (grepl("^(Participant\\s+)?ID\\s*:\\s*P\\d+", p, ignore.case = TRUE)) {
      flush_entry()
      current_id    <- regmatches(p, regexpr("P\\d+", p))
      current_date  <- NA_character_
      current_notes <- character(0)
      in_notes      <- FALSE
    } else if (grepl("^(Session\\s+)?Date\\s*:", p, ignore.case = TRUE)) {
      current_date <- trimws(sub("^(Session\\s+)?Date\\s*:\\s*", "", p,
                                 ignore.case = TRUE))
    } else if (grepl("^Notes\\s*$", p, ignore.case = TRUE)) {
      in_notes <- TRUE
    } else if (in_notes && !is.na(current_id)) {
      current_notes <- c(current_notes, p)
    }
  }
  flush_entry()
  dplyr::bind_rows(entries)
}

# officer needed for docx parsing
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
library(officer, quietly = TRUE)

notes_docx <- parse_session_notes(f_notes)
cat(sprintf("Session notes parsed: %d participants\n", nrow(notes_docx)))

## 4c. Cancellation status (from Sheet2 of tracker) ---------------------------
# P076 = cancelled at attendance; P079, P081 = cancelled before attendance
cancellations <- tibble::tibble(
  participant_id     = c("P076", "P079", "P081"),
  cancellation_status = c("cancelled_at_attendance",
                          "cancelled_before_attendance",
                          "cancelled_before_attendance")
)

## 4d. Merge all note sources --------------------------------------------------
# Base: all P### from RunCorrection (most complete list of attempted participants)
run_notes <- run_corr_wide |>
  dplyr::full_join(notes_docx, by = "participant_id") |>
  dplyr::left_join(cancellations, by = "participant_id") |>
  dplyr::mutate(
    cancellation_status = dplyr::if_else(
      is.na(cancellation_status), "attended", cancellation_status
    ),
    # Participants only in docx (not RunCorrection) have NA flag → set FALSE
    flag_review = dplyr::if_else(is.na(flag_review), FALSE, flag_review),
    log_script_version = "1.0",
    log_date_cleaned   = as.character(Sys.Date())
  ) |>
  dplyr::arrange(participant_id) |>
  dplyr::select(
    participant_id, session_date, cancellation_status,
    run1_loc1_note, run2_loc2_note, run3_breath1_note,
    run4_breath2_note, run5_auth1_note,
    general_notes_tracker, general_notes_docx, flag_review,
    log_script_version, log_date_cleaned
  )

run_notes <- strip_emails(run_notes)
readr::write_csv(run_notes, file.path(dir_out, "run_notes.csv"))
cat(sprintf("run_notes.csv written: %d rows\n", nrow(run_notes)))

# ── 5. Summary log ───────────────────────────────────────────────────────────
cat("\n── Summary ──────────────────────────────────────────────────────\n")
cat(sprintf("screener_clean.csv : %d rows, %d cols\n",
            nrow(screener_clean), ncol(screener_clean)))
cat(sprintf("postscan_clean.csv : %d rows, %d cols\n",
            nrow(postscan_clean), ncol(postscan_clean)))
cat(sprintf("run_notes.csv      : %d rows, %d cols\n",
            nrow(run_notes), ncol(run_notes)))
cat(sprintf("Screener P### match rate : %.1f%%\n",
            100 * mean(!is.na(screener_clean$participant_id))))
cat(sprintf("Post-scan P### match rate: %.1f%%\n",
            100 * mean(!is.na(postscan_clean$participant_id))))
cat(sprintf("Cancellations: %d (%d at attendance, %d before)\n",
            nrow(cancellations),
            sum(cancellations$cancellation_status == "cancelled_at_attendance"),
            sum(cancellations$cancellation_status == "cancelled_before_attendance")))
cat(sprintf("Run ended: %s\n", Sys.time()))