# intero2023_behaviour.R
# Behavioural analysis of the change-detection task.
#
# THE DESIGN CONSTRAINS WHAT CAN BE ASKED, so that is settled first.
#
#   14 paced trials per participant (7 per run x 2 runs), of which a MEDIAN OF
#   2 are No Change. A per-participant false-alarm rate estimated from 2 trials
#   can only take the values 0, 0.5 or 1, so per-participant d' and criterion
#   are not estimable and are not reported as such. Everything that depends on
#   the No Change trials is therefore pooled, or fitted hierarchically with
#   participant as a random effect, which is what borrows strength across the
#   sparse cells instead of pretending each participant supports its own SDT.
#
# THE TASK IS A THREE-CATEGORY IDENTIFICATION, AND IS ANALYSED AS ONE.
# Each trial ends in one of Faster / No Change / Slower, so chance is 1/3 for an
# unbiased responder, not 1/2.
#
# An earlier version of this script collapsed the responses into "changed" vs
# "not" and reported a yes/no d' of +0.07 as though it showed no sensitivity.
# That was wrong. Collapsing imposes a detection stage the task does not have,
# and discards the direction information that carries nearly all the signal:
# with zero sensitivity, direction on true-change trials would sit at 0.5, and
# it sits at 0.78. Three-category accuracy is 0.525 against a 1/3 baseline
# (kappa 0.271, X2(4) = 141.7, p = 1e-29) and the ordinal sensitivity term is
# 0.698 probit units per stimulus step, z = 13.2.
#
# The collapsed statistic is retained below as a BIAS diagnostic, which is all
# it ever measured: P(say changed) is ~0.70 whether or not anything changed.
# The reason is visible in the confidence data -- "No Change" is the LOWEST
# confidence response (5.03) rather than a positive judgement that nothing
# happened -- so it behaves as an "I don't know", and a yes/no SDT built on it
# is measuring willingness to commit, not sensitivity.
#
# THE FREE-BREATHING TEST IS THE NOVEL ONE. Free epochs have no instructed
# change, so until the belt was processed there was no ground truth at all.
# actual_period_ratio now supplies one, which makes "did the reported direction
# track what the breathing actually did" answerable for the first time.
#
# Produces: Repo/Data/behaviour_{summary,subject}.csv, Repo/Figures/behaviour.png
#
# Script version: 1.0
# Author: RADlab

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(purrr); library(stringr)
  library(lme4); library(ordinal)
})

dir_root <- "I:/Shared drives/MRI_Intero2023"
dir_out  <- file.path(dir_root, "Repo", "Data")
dir_fig  <- file.path(dir_root, "Repo", "Figures")
dir_R    <- file.path(dir_root, "Repo", "R")
for (f in c("intero_ids.R")) source(file.path(dir_R, f))
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Run started: %s\n", Sys.time()))

b <- read_csv(file.path(dir_out, "breathing_trials.csv"), show_col_types = FALSE,
              progress = FALSE) %>%
  mutate(across(c(confidence, arousal, valence_run, actual_period_ratio,
                  instructed_period_ratio), ~ suppressWarnings(as.numeric(.x))),
         salience = ifelse(is.na(salience) | salience == "n/a", NA_character_, salience))

paced <- b %>% dplyr::filter(trial_type == "paced_breathing", !is.na(response))
free  <- b %>% dplyr::filter(trial_type == "free_breathing",  !is.na(response))

cat("\n=====================================================================\n")
cat("1. WHAT THE DESIGN SUPPORTS\n")
cat("=====================================================================\n")
cat(sprintf("paced trials: %d (%d participants, %d runs)\n", nrow(paced),
            n_distinct(paced$subject), n_distinct(paste(paced$subject, paced$bids_entity))))
nc <- paced %>% dplyr::filter(change_direction == "No Change") %>% count(subject)
cat(sprintf("No Change trials per participant: median %d, max %d  <-- the SDT noise trials\n",
            median(nc$n), max(nc$n)))
cat("  => per-participant d'/criterion are NOT estimable; pooled and hierarchical only.\n")

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("2a. THREE-CATEGORY IDENTIFICATION (the primary analysis, chance = 1/3)\n")
cat("=====================================================================\n")
p3 <- paced %>%
  mutate(stim = factor(change_direction, levels = c("Faster", "No Change", "Slower")),
         resp = factor(response, levels = c("Faster", "No Change", "Slower")),
         correct3 = as.integer(as.character(change_direction) == as.character(response)))
tb3 <- table(p3$stim, p3$resp)
cat("confusion matrix (rows = stimulus, cols = response):\n"); print(addmargins(tb3))
cat("\nrow proportions:\n"); print(round(prop.table(tb3, 1), 3))
cat(sprintf("\noverall accuracy %.3f (%d/%d), binomial vs 1/3: p = %s\n",
            mean(p3$correct3), sum(p3$correct3), nrow(p3),
            signif(binom.test(sum(p3$correct3), nrow(p3), 1/3, alt = "greater")$p.value, 3)))
print(as.data.frame(p3 %>% group_by(stim) %>% summarise(n = n(),
  acc = round(mean(correct3), 3),
  p_vs_chance = signif(binom.test(sum(correct3), n(), 1/3, alt = "greater")$p.value, 3),
  .groups = "drop")))
exp_acc <- sum(rowSums(tb3) / sum(tb3) * colSums(tb3) / sum(tb3))
cat(sprintf("\nguessing with THIS response distribution would give %.3f; observed %.3f\n",
            exp_acc, mean(p3$correct3)))
cat(sprintf("Cohen's kappa = %.3f ; stimulus-response independence X2(%d) = %.1f, p = %s\n",
            (mean(p3$correct3) - exp_acc) / (1 - exp_acc),
            chisq.test(tb3)$parameter, chisq.test(tb3)$statistic,
            signif(chisq.test(tb3)$p.value, 3)))

cat("\n-- ordinal SDT: sensitivity separated from criterion placement --\n")
cat("   responses are ordered on one axis, so a cumulative-link model gives\n")
cat("   sensitivity as the stimulus slope and the No-Change band as the gap\n")
cat("   between the two thresholds.\n")
p3$stim_num <- c(Faster = -1, `No Change` = 0, Slower = 1)[as.character(p3$stim)]
p3$resp_ord <- factor(p3$response, levels = c("Faster", "No Change", "Slower"), ordered = TRUE)
mo <- ordinal::clmm(resp_ord ~ stim_num + (1 | subject), data = p3, link = "probit")
print(round(coef(summary(mo)), 4))
th <- mo$Theta; if (is.null(dim(th))) th <- matrix(th, nrow = 1)
cat(sprintf("   sensitivity = %.3f probit units per step; No-Change band = %.3f SD\n",
            coef(mo)["stim_num"], th[1, 2] - th[1, 1]))

cat("\n-- can they identify NO CHANGE specifically? --\n")
ncr <- p3 %>% dplyr::filter(stim == "No Change")
cat(sprintf("   accuracy on No Change trials %.3f vs chance 1/3, p = %s  <-- at chance\n",
            mean(ncr$correct3), signif(binom.test(sum(ncr$correct3), nrow(ncr), 1/3)$p.value, 3)))
cat(sprintf("   P(respond 'No Change') = %.3f while No Change trials are %.3f of trials\n",
            mean(p3$resp == "No Change"), mean(p3$stim == "No Change")))

cat("\n=====================================================================\n")
cat("2b. THE COLLAPSED YES/NO STATISTIC -- a BIAS diagnostic, not sensitivity\n")
cat("=====================================================================\n")
det <- paced %>%
  mutate(signal = change_direction != "No Change",
         said_change = response != "No Change")
tab <- table(signal = det$signal, said = det$said_change)
H <- sum(det$said_change & det$signal) / sum(det$signal)
FA <- sum(det$said_change & !det$signal) / sum(!det$signal)
zc <- function(p, n) qnorm(pmin(pmax((p * n + 0.5) / (n + 1), 1e-6), 1 - 1e-6))
dprime <- zc(H, sum(det$signal)) - zc(FA, sum(!det$signal))
crit   <- -0.5 * (zc(H, sum(det$signal)) + zc(FA, sum(!det$signal)))
cat(sprintf("  hit rate  P(say change | change)    = %.3f  (%d trials)\n", H, sum(det$signal)))
cat(sprintf("  FA  rate  P(say change | no change) = %.3f  (%d trials)\n", FA, sum(!det$signal)))
cat(sprintf("  pooled d' = %+.3f   criterion = %+.3f\n", dprime, crit))

# bootstrap over participants, because trials within a person are not independent
set.seed(1)
subs <- unique(det$subject)
bs <- replicate(2000, {
  s <- sample(subs, length(subs), replace = TRUE)
  d <- bind_rows(lapply(s, function(x) det[det$subject == x, ]))
  h <- sum(d$said_change & d$signal) / max(1, sum(d$signal))
  f <- sum(d$said_change & !d$signal) / max(1, sum(!d$signal))
  zc(h, sum(d$signal)) - zc(f, sum(!d$signal))
})
cat(sprintf("  d' 95%% CI (participant bootstrap): [%+.3f, %+.3f]\n",
            quantile(bs, .025), quantile(bs, .975)))
cat("  => P(say changed) is ~0.70 regardless of stimulus. Read this as a strong\n")
cat("     liberal bias in WHETHER TO COMMIT to a direction, NOT as absent\n")
cat("     sensitivity -- section 2a shows sensitivity is large. The two coexist\n")
cat("     because 'No Change' is the lowest-confidence response (see section 4),\n")
cat("     so it acts as 'I don't know' rather than 'nothing happened'.\n")

cat("\n-- confidence by outcome, which is what shows that --\n")
print(as.data.frame(paced %>% mutate(kind = case_when(
    change_direction == "No Change" & response == "No Change" ~ "NC correct",
    change_direction == "No Change" ~ "NC error (said a direction)",
    change_direction == response ~ "direction correct",
    response == "No Change" ~ "said No Change on a change trial",
    TRUE ~ "wrong direction")) %>%
  group_by(kind) %>% summarise(n = n(), confidence = round(mean(confidence, na.rm = TRUE), 2),
                               .groups = "drop")))

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("3. DIRECTION: given a directional response, was it right?\n")
cat("=====================================================================\n")
dir_d <- paced %>%
  dplyr::filter(change_direction %in% c("Faster", "Slower"),
                response %in% c("Faster", "Slower")) %>%
  mutate(correct = as.integer(response == change_direction),
         salience = factor(salience, levels = c("low", "high")),
         change_direction = factor(change_direction, levels = c("Faster", "Slower")))
cat(sprintf("  trials: %d (%d participants)\n", nrow(dir_d), n_distinct(dir_d$subject)))
cat(sprintf("  overall correct direction: %.3f\n", mean(dir_d$correct)))
print(as.data.frame(dir_d %>% group_by(change_direction, salience) %>%
  summarise(n = n(), correct = round(mean(correct), 3), .groups = "drop")))

m_dir <- glmer(correct ~ salience * change_direction + (1 | subject),
               data = dir_d, family = binomial,
               control = glmerControl(optimizer = "bobyqa"))
cat("\n  mixed logistic, correct ~ salience * direction + (1|subject):\n")
print(round(summary(m_dir)$coefficients, 4))
m0 <- glmer(correct ~ change_direction + (1 | subject), data = dir_d, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
m1 <- glmer(correct ~ salience + change_direction + (1 | subject), data = dir_d,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("\n  salience main effect (LRT):\n"); print(anova(m0, m1)[2, c("Chisq", "Df", "Pr(>Chisq)")])

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("4. CONFIDENCE: does it track being right? (metacognitive sensitivity)\n")
cat("=====================================================================\n")
cf <- dir_d %>% dplyr::filter(is.finite(confidence))
cat(sprintf("  confidence 1-7, median %.1f\n", median(cf$confidence)))
print(as.data.frame(cf %>% group_by(correct) %>%
  summarise(n = n(), conf = round(mean(confidence), 2), sd = round(sd(confidence), 2),
            .groups = "drop")))
m_conf <- glmer(correct ~ scale(confidence) + (1 | subject), data = cf, family = binomial,
                control = glmerControl(optimizer = "bobyqa"))
cat("\n  correct ~ z(confidence) + (1|subject):\n")
print(round(summary(m_conf)$coefficients, 4))
# type-2 AUC per participant: can confidence separate correct from incorrect?
auc1 <- function(x, y) { a <- x[y == 1]; c0 <- x[y == 0]
  if (!length(a) || !length(c0)) return(NA_real_)
  mean(outer(a, c0, ">") + 0.5 * outer(a, c0, "==")) }
m2 <- cf %>% group_by(subject) %>%
  summarise(n = n(), auc = auc1(confidence, correct), .groups = "drop") %>%
  dplyr::filter(is.finite(auc), n >= 6)
cat(sprintf("\n  type-2 AUC per participant: median %.3f, n=%d, %.0f%% above 0.5\n",
            median(m2$auc), nrow(m2), 100 * mean(m2$auc > 0.5)))
cat(sprintf("  one-sample test vs 0.5: t=%.2f, p=%.4f\n",
            t.test(m2$auc, mu = 0.5)$statistic, t.test(m2$auc, mu = 0.5)$p.value))

cat("\n-- arousal --\n")
ar <- paced %>% dplyr::filter(is.finite(arousal))
print(as.data.frame(ar %>% group_by(change_direction) %>%
  summarise(n = n(), arousal = round(mean(arousal), 2), .groups = "drop")))
cat(sprintf("  arousal-confidence correlation (within participant, mean r): %.3f\n",
            mean((paced %>% group_by(subject) %>%
              summarise(r = suppressWarnings(cor(confidence, arousal,
                        use = "complete.obs")), .groups = "drop"))$r, na.rm = TRUE)))

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("5. FREE BREATHING: does the report track what breathing ACTUALLY did?\n")
cat("=====================================================================\n")
fm <- free %>% dplyr::filter(resp_measured %in% TRUE, is.finite(actual_period_ratio))
cat(sprintf("  free epochs with measured respiration: %d (%d participants)\n",
            nrow(fm), n_distinct(fm$subject)))
cat("  actual_period_ratio > 1 means the period LENGTHENED, i.e. breathing slowed.\n")
print(as.data.frame(fm %>% group_by(response) %>%
  summarise(n = n(), ratio_med = round(median(actual_period_ratio), 3),
            ratio_mean = round(mean(actual_period_ratio), 3), .groups = "drop")))
fmd <- fm %>% dplyr::filter(response %in% c("Faster", "Slower")) %>%
  mutate(said_slower = as.integer(response == "Slower"),
         lr = log(actual_period_ratio))
if (nrow(fmd) > 30) {
  m_free <- glmer(said_slower ~ scale(lr) + (1 | subject), data = fmd, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
  cat("\n  said_slower ~ z(log measured period ratio) + (1|subject):\n")
  print(round(summary(m_free)$coefficients, 4))
  cat("  a POSITIVE slope means reports track the real direction of change.\n")
}

cat("\n-- same test on PACED trials, measured vs reported (not instructed) --\n")
pm <- paced %>% dplyr::filter(resp_measured %in% TRUE, is.finite(actual_period_ratio),
                              response %in% c("Faster", "Slower")) %>%
  mutate(said_slower = as.integer(response == "Slower"), lr = log(actual_period_ratio))
if (nrow(pm) > 30) {
  m_p <- glmer(said_slower ~ scale(lr) + (1 | subject), data = pm, family = binomial,
               control = glmerControl(optimizer = "bobyqa"))
  print(round(summary(m_p)$coefficients, 4))
}

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("5b. WHERE DOES SALIENCE ACT: DETECTION OR DIRECTION?\n")
cat("=====================================================================\n")
cat("Salience is not significant on direction accuracy above, yet it looks\n")
cat("strong on the combined measure. Both are true and they are not in\n")
cat("conflict -- the combined measure counts a 'no change' response as a\n")
cat("direction error, so it carries the detection step inside it.\n\n")
cmb <- paced %>% dplyr::filter(change_direction %in% c("Faster", "Slower")) %>%
  mutate(correct = as.integer(response == change_direction),
         salience = factor(salience, levels = c("low", "high")),
         change_direction = factor(change_direction, levels = c("Faster", "Slower")))
print(as.data.frame(cmb %>% group_by(change_direction, salience) %>%
  summarise(n = n(), correct = round(mean(correct), 3), .groups = "drop")))
c0 <- glmer(correct ~ change_direction + (1 | subject), cmb, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
c1 <- glmer(correct ~ salience + change_direction + (1 | subject), cmb, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
cat("\n  salience LRT on COMBINED accuracy:\n"); print(anova(c0, c1)[2, c("Chisq", "Df", "Pr(>Chisq)")])
sd_ <- paced %>% dplyr::filter(change_direction %in% c("Faster", "Slower")) %>%
  mutate(said_change = as.integer(response != "No Change"),
         salience = factor(salience, levels = c("low", "high")))
print(as.data.frame(sd_ %>% group_by(salience) %>%
  summarise(n = n(), p_say_change = round(mean(said_change), 3), .groups = "drop")))
s0 <- glmer(said_change ~ 1 + (1 | subject), sd_, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
s1 <- glmer(said_change ~ salience + (1 | subject), sd_, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
cat("\n  salience LRT on P(say changed):\n"); print(anova(s0, s1)[2, c("Chisq", "Df", "Pr(>Chisq)")])
cat("  => salience acts on NOTICING, not on identifying which way.\n")

cat("\n=====================================================================\n")
cat("5c. IS THE PACED 'TRACKING' JUST THE PACER THEY WERE WATCHING?\n")
cat("=====================================================================\n")
pmi <- paced %>% dplyr::filter(resp_measured %in% TRUE, is.finite(actual_period_ratio),
        is.finite(instructed_period_ratio), response %in% c("Faster", "Slower")) %>%
  mutate(said_slower = as.integer(response == "Slower"),
         lr = scale(log(actual_period_ratio))[, 1],
         li = scale(log(instructed_period_ratio))[, 1])
cat(sprintf("  trials=%d participants=%d  cor(measured, instructed)=%.3f\n",
            nrow(pmi), n_distinct(pmi$subject), cor(pmi$lr, pmi$li)))
cat("\n  both predictors together:\n")
print(round(summary(glmer(said_slower ~ lr + li + (1 | subject), pmi, family = binomial,
      control = glmerControl(optimizer = "bobyqa")))$coefficients, 4))

cat("\n-- how big is the per-epoch measurement noise vs the imposed signal? --\n")
sp <- b %>% dplyr::filter(resp_measured %in% TRUE, is.finite(actual_period_ratio))
print(as.data.frame(sp %>% group_by(trial_type) %>% summarise(n = n(),
  med = round(median(actual_period_ratio), 3),
  sd_log = round(sd(log(actual_period_ratio)), 3), .groups = "drop")))
print(as.data.frame(sp %>% dplyr::filter(trial_type == "paced_breathing") %>%
  group_by(change_direction) %>% summarise(n = n(),
  med = round(median(actual_period_ratio), 3),
  sd_log = round(sd(log(actual_period_ratio)), 3), .groups = "drop")))
cat("\n  The No Change sd_log IS the noise floor. Read the measured-change\n")
cat("  coefficient against it before concluding anything about interoception:\n")
cat("  a noisy predictor competing with a clean one is attenuated toward zero\n")
cat("  (regression dilution), so a null there is NOT evidence of no tracking.\n")

# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("6. DOES CONDITIONING ON VERIFIED ADHERENCE CHANGE ANYTHING?\n")
cat("=====================================================================\n")
adh <- dir_d %>% dplyr::filter(adherence_verified %in% TRUE)
cat(sprintf("  direction trials with verified adherence: %d (%d participants)\n",
            nrow(adh), n_distinct(adh$subject)))
cat(sprintf("  correct direction: %.3f  (vs %.3f in the full set)\n",
            mean(adh$correct), mean(dir_d$correct)))
deta <- det %>% dplyr::filter(adherence_verified %in% TRUE)
if (sum(!deta$signal) > 5) {
  Ha <- sum(deta$said_change & deta$signal) / sum(deta$signal)
  FAa <- sum(deta$said_change & !deta$signal) / sum(!deta$signal)
  cat(sprintf("  detection hit %.3f  FA %.3f  d' %+.3f  (n_noise=%d)\n", Ha, FAa,
              zc(Ha, sum(deta$signal)) - zc(FAa, sum(!deta$signal)), sum(!deta$signal)))
}

# ---------------------------------------------------------------------------
subj <- paced %>% group_by(subject) %>%
  summarise(n_paced = n(),
            p_say_change = mean(response != "No Change"),
            conf = mean(confidence, na.rm = TRUE),
            arousal = mean(arousal, na.rm = TRUE), .groups = "drop") %>%
  left_join(dir_d %>% group_by(subject) %>%
              summarise(n_dir = n(), p_correct_dir = mean(correct), .groups = "drop"),
            by = "subject") %>%
  left_join(m2 %>% select(subject, type2_auc = auc), by = "subject")
write_csv(subj, file.path(dir_out, "behaviour_subject.csv"))
write_csv(tibble(measure = c("hit", "fa", "dprime", "criterion", "p_correct_direction"),
                 value = c(H, FA, dprime, crit, mean(dir_d$correct))),
          file.path(dir_out, "behaviour_summary.csv"))

png(file.path(dir_fig, "behaviour.png"), width = 1500, height = 460, res = 110)
op <- par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
barplot(c(Hit = H, FA = FA), col = c("steelblue", "firebrick"), ylim = c(0, 1),
        main = sprintf("detection (d' = %+.2f)", dprime), ylab = "P(say changed)")
with(dir_d %>% group_by(change_direction, salience) %>% summarise(p = mean(correct), .groups = "drop"),
     barplot(p, names.arg = paste(change_direction, salience), col = "grey70", ylim = c(0, 1),
             main = "direction correct", las = 2, ylab = "P(correct)"))
abline(h = 0.5, col = "red", lty = 2)
hist(m2$auc, breaks = 15, col = "grey80", border = "white", main = "type-2 AUC (confidence)",
     xlab = "AUC"); abline(v = 0.5, col = "red", lwd = 2)
boxplot(actual_period_ratio ~ response, data = fm, col = "grey85", log = "y",
        main = "free: measured vs reported", ylab = "measured period ratio", xlab = "")
abline(h = 1, col = "red")
par(op); dev.off()

cat(sprintf("\nRun ended: %s\n", Sys.time()))
