# Localizer: response and adherence quality by condition (descriptive)

Logged in `Preregistration/localizer_deviations.md` (2026-11-28) before the
per-participant values were computed.
- Script: `Scripts/intero2023_loc_response_summary.py`.
- Tables: `Analysis/results/EXPLORATORY_LOC_response_summary.csv` (group) and
  `EXPLORATORY_LOC_response_summary_participants.csv` (per participant).

**Scope.** Registered run set: 102 runs, 54 participants, attempt of record.

**How values are reported.** Each measure is computed per participant. The
tables give the group median [IQR] and the number of participants. No tests are
run.

## How responses were logged

- **Intero.** Presses were logged cleanly.
- **Extero, Match and Self.** Each log row (a pacer half-cycle, or a trait word)
  keeps only the **last** key pressed in its response window of about 2 s.
  - The scanner trigger arrives as key `5` every 2 s, so it usually
    overwrites the participant's press.
  - A real button is stored only when it was pressed *after* that window's
    trigger.
  - The stored buttons are therefore a selected subset. They show *whether*
    presses were correct when we can see them, not how well people tracked.
- **Physiology files.** The Biopac files have no button channel.

## Summary by condition

| Condition | Measure | Median [IQR] | N |
|---|---|---|---|
| **Intero** (press at each inhale and exhale onset) | Presses per block | 11 [8.6, 12.9] | 54 |
| | Compliant blocks (≥ 4 presses) | 100% [100, 100] | 54 |
| | Key alternation (1 then 2) | 1.00 [1.00, 1.00] | 53 |
| | Press–breath phase-locking (PLV) | **.79** [.71, .87]; chance .26 | 38 |
| | Press-pair score (cross-validated) | **.79** [.60, .88]; chance .39 | 36 |
| | Press lead before the belt extremum | 0.84 s [0.67, 1.03] | 36 |
| **Extero** (press at each circle phase change) | Trials with a stored button | 24% [17, 31] | 54 |
| | Trials where the trigger overwrote the response | 76% [69, 83] | 54 |
| | Stored buttons that were the correct key | **100%** [100, 100]; pooled 97.5% of 606 | 54 |
| | Latency of stored buttons after phase onset | 0.48 s [0.37, 0.69] (selected; not a tracking measure) | 53 |
| | Breathing phase-locked to the circle (not instructed) | .56 [.36, .68]; virtual-pacer null .32 | 40 |
| **Match** (breathe with the circle, press at phase changes) | Trials with a stored button | 21% [15, 29] | 54 |
| | Trials where the trigger overwrote the response | 77% [71, 83] | 54 |
| | Stored buttons that were the correct key | **100%** [100, 100]; pooled 97.6% of 550 | 53 |
| | Breathing phase-locked to the circle (instructed) | **.87** [.77, .91]; null .32 | 40 |
| **Self** (trait-word judgements) | Trials with a stored answer | 0% [0, 0]; 30 answers in total | 54 |
| **Breathing runs** (reference) | Paced-epoch breathing PLV to the pacer | .82 [.78, .89]; free-epoch virtual .31 | 33 |
| | Run-level adherence (`adhere_r`) | .43 [.20, .57] | 35 |

## Reading

- **Intero.**
  - Participants pressed reliably: a median of 11 presses per 30 s block,
    with nearly every block compliant and strict 1/2 alternation.
  - The presses were tightly locked to the breath (PLV .79, three times
    chance).
  - They led the belt extremum by about 0.85 s, which fits pressing at the
    felt start of each inhale and exhale.
- **Extero and Match.**
  - Where a press is visible, it is almost always the correct key for the
    circle's phase, so participants were doing the task.
  - The visible presses cannot be turned into a tracking-accuracy score:
    presses made before the trigger, including on-time presses, are
    overwritten.
  - Breathing adherence is clear: in Match, breathing locks to the circle
    (PLV .87).
  - In Extero, breathing drifts partly into phase with the circle without
    instruction (.56 against a null of .32). Extero therefore does not fully
    separate breathing from the visual pacer.
- **Self.**
  - Answers were almost entirely overwritten, so the trait-word judgements
    cannot be scored.
- **Breathing runs.**
  - Paced-epoch adherence is as high as in localizer Match (PLV .82).
