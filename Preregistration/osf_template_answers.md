# OSF Preregistration — field-by-field answers (copy-paste ready)

Template: OSF Preregistration (osf.io registration form). Each heading below is
one template field, in template order. Text under it is paste-ready. Items
marked [TO CONFIRM] need the author's input before freezing.

Synced to `intero2023_imaging_prereg.md` **v0.4** (2026-11-25).

---

## Title

Neural correlates of awareness of respiratory change during paced breathing (MRI_Intero2023)

## Description

A single-session fMRI study in which participants breathed with a visual pacer
whose rate increased, decreased, or held constant on each trial (with abrupt vs
gradual change salience), interleaved with free-breathing periods, and reported
the perceived direction of change with confidence and arousal ratings.
Respiration was recorded with a chest belt and expired CO2 with a capnograph.
This registration covers the functional imaging analysis only: all behavioural
and physiological processing is complete and disclosed; no functional image
statistic has been computed. Hypotheses include a direct replication of Farb,
Zuo & Price (2023, eNeuro) — insula/opercular deactivation under exogenous
respiratory pacing — and a conceptual extension of that paper's MAIA–anterior
cingulate result to this task's pacing contrast.

## Contributors

Norman A. S. Farb (sole author).

## License

CC-By Attribution 4.0 International (confirmed 2026-08-20; matches the IEAT
paper's open-access licence).

## Subject

Social and Behavioral Sciences → Psychology → Cognition and Perception;
Life Sciences → Neuroscience and Neurobiology → Cognitive Neuroscience.

## Tags

interoception; respiration; breath awareness; fMRI; anterior cingulate;
insula; PetCO2; RVT; metacognition; replication

---

## Research questions or hypotheses

Each hypothesis is tested in one of three first-level models (A awareness,
B confidence, C salience); awareness and the confidence modulator are never
estimated in the same design.

H1 (replication; Model A; N = 47). Exogenously paced breathing, relative to
free breathing, produces relative DEACTIVATION of the insula/operculum pathway
and somatomotor cortex, replicating the Active Matching result of Farb, Zuo &
Price (2023, eNeuro 10(6); "IEAT"). "Deactivation" means lower activity during
paced than free epochs; the design has no rest baseline. The volitional-control
account predicts the opposite (activation); a reversed result will be reported
as a failed replication, not reframed.

H3 (authoritative awareness test; Model A; N = 23). Change trials answered
with the correct direction ("aware") elicit greater activity than change
trials answered "No Change" ("unaware") in left anterior insula, right
anterior insula and dACC (ROI-restricted, Bonferroni ×3). H3 is authoritative
because confidence and awareness covary behaviourally (within-participant
r = .316 median), so a parametric confidence effect cannot arbitrate between
them. Registered limitation: the 23 participants with ≥ 3 trials in both cells
over-represent less frequent detectors (median unaware share .36 vs .19); the
estimate is not generalised beyond this subsample.

H2 (supportive awareness test; Model B; N = 44). Trial-by-trial confidence on
paced trials correlates positively with activity in the same three ROIs.
Model B's confidence effect absorbs the awareness effect rather than being
adjusted for it. Registered in advance: a null H3 alongside a positive H2 is
NOT evidence that confidence rather than awareness drives the effect, because
that comparison (N = 23 vs 44) is confounded with power.

H4 (Model C; N = 43). High-salience (step) changes elicit greater activity
than low-salience (ramp) changes in the same three ROIs, consistent with the
behavioural finding that salience affects noticing (p = .043) but not
direction discrimination (p = .26). Supported if the effect survives in at
least one ROI; falsified if in none. Salience is not crossed with awareness
(cells of ~1 trial).

H5. Adding physiological regressors (RVT; PetCO2) attenuates but does not
abolish the awareness-related effects. The confirmatory test is an
equivalence criterion on the H2 confidence effect (see Inference criteria).
The same comparison for H3 is descriptive only: 13 of its 23 participants have
PetCO2 and 15 have Biopac-belt RVT, too few for an equivalence claim.

H6 (conceptual extension; Model A; N = 46; two-sided). MAIA interoceptive
sensibility moderates paced-vs-free activity in the anterior cingulate, tested
in a 10 mm sphere at the IEAT's coordinate (MNI −4, +12, +36). In the IEAT,
higher MAIA predicted reduced ACC deactivation on Active Interoception vs
Active Exteroception; this task has no exteroceptive condition, and the IEAT
did not relate MAIA to its pacing contrast, so the sign is not fixed by prior
work and the test is two-sided.

## Foreknowledge of data or evidence

Selection (confirmed 2026-08-20): the option stating that the data exist and
the author has accessed and analyzed some of the data relevant to this
registration. See explanation below for exactly what has and has not been
touched.

## Explanation of foreknowledge

Data collection is complete (September 2023 – June 2024). Everything EXCEPT
the functional imaging has been analysed before this registration: task logs
and event timing; all trial-level behaviour (three-category identification
accuracy 0.525 vs chance 1/3; ordinal sensitivity 0.698 probit units/step;
confidence type-2 AUC 0.661); respiration-belt processing including pacer
adherence; capnography/PetCO2 extraction; RVT; and all quality control and
exclusion decisions. The functional images have not been examined beyond
automated preprocessing (fMRIPrep, and MRIQC image-quality metrics used for
run QC): no first-level model has been fitted, no contrast image exists, no
statistical map has been computed or viewed. Before freezing, the design
matrices themselves were built and checked (rank, collinearity, contrast
efficiency, event coverage, trial counts) — properties of the timing and
confound sheets; no BOLD voxel was read. Two narrow, documented exceptions:
fMRIPrep's global_signal column was used to verify the physiological
regressors carry the physiologically correct sign, and its
framewise_displacement column was used to apply the preregistered motion
exclusion. No task regressor has been related to any imaging measure.
Consequence, stated plainly: H2–H4 are informed by known behavioural results;
they are untested predictions about brain data but are not blind. H1 takes its
direction, and H6 its region, from an independent published sample (IEAT).
The full analysis plan with dated decision history is attached
(intero2023_imaging_prereg.md) and git-versioned.

---

# Research Design

## Study type

Randomized Experiment (conditions randomized within participant at the trial
level; no between-participant manipulation).

## Intention for causal interpretation

Direct inference for the within-participant task manipulations (pacing,
change direction, salience). No causal interpretation for the
individual-difference moderator (H6, MAIA) or for awareness contrasts (H2/H3),
which condition on the participant's response.

## Blinding of experimental treatments

Subject blinding: participants could not anticipate the direction or salience
of an upcoming trial. No researcher or analyst blinding.

## Additional blinding

Analyst blinding to behavioural outcomes is impossible (see Foreknowledge);
the analyst has seen no functional imaging statistic. Post-registration
deviations will be logged with date and whether the affected result had been
seen (deviations.md, git-versioned).

## Study design

Within-participant task fMRI, single session, five BOLD runs (2 localizer,
2 breathing, 1 authenticity); this registration concerns the two breathing
runs. Each breathing run: 7 cycles of [paced-breathing epoch (~20 s, visual
pacer; instructed change Faster / Slower / No Change; for changes, salience
high = step at midpoint vs low = linear ramp; endpoint ±25–35% of period) →
rating window (direction of change, confidence 1–7, arousal 1–7) →
free-breathing epoch (~20 s) → rating window]. These events cover ~94% of
each run, so there is no rest baseline; only between-condition contrasts are
interpreted. Change direction and salience orders randomized within run. Two
protocol versions exist (v1 n≈9 with per-epoch valence ratings and longer
runs; v2 the remainder); protocol version is a registered group-level
covariate (sensitivity S5). Acquisition: 3T Siemens, EPI, TR = 2 s, 2 mm
isotropic; PEPOLAR fieldmaps; fMRIPrep 24.1.1, MNI152NLin2009cAsym.
Concurrent physiological recording: Biopac respiration belt and capnograph
(2 kHz, downsampled to 100 Hz), scanner-triggered; Siemens PMU respiration
(400 Hz) as a validation stream.

## Randomization

Trial-level randomization of change direction (Faster / Slower / No Change)
and salience (high / low) within participant, implemented by the stimulus
software at session time. No between-participant randomization.

---

# Sampling

## Data collection procedures

Data are already collected (September 2023 – June 2024, single site,
University of Toronto). Community adult sample; standard MRI safety screening
via Qualtrics before the session; questionnaire battery including the 37-item
MAIA-2 (coded 1–6, converted to the canonical 0–5 metric).
[TO CONFIRM: recruitment channels, compensation, REB protocol number, and any
inclusion criteria beyond MRI safety.]

## Sample size

Collected: 54 non-pilot participants with usable task logs (105 breathing
runs). After preregistered exclusions (below): 90 breathing runs from 47
participants. Analysable samples: H1 N = 47; H2 N = 44; H3 N = 23; H4 N = 43;
H6 N = 46; sensitivity S1 (PetCO2) N = 25; S2 (Biopac RVT) N = 30; S2b (RVT
with PMU fallback) N = 41; S4 (verified pacer adherence) N = 26; S6 (two-run participants only) N = 43. Trials per
participant (median): 7 aware, 3 unaware, 1 wrong-direction, 2 no-change,
14 free.

## Sample size rationale

The sample was collected before this registration; no further collection will
occur, so size is characterized rather than justified. Smallest detectable
effect at 80% power (one-sample t): N = 47 → d = 0.42 at α = .05 (ROI),
d = 0.64 at cluster-forming p < .001; N = 44 → 0.43 / 0.67; N = 23 →
0.61 / 0.98. H6 (two-sided, N = 46) has 80% power for |r| ≥ 0.40 and ≈ 44–59%
at the IEAT's ROI-level |r| ≈ 0.27–0.32; this is stated so a null is not read
as evidence against the IEAT result. Equivalence bounds for H5 are set by what
these Ns support (see Inference criteria).

## Starting and stopping rules

Collection is complete; no starting or stopping rule remains operative.
[TO CONFIRM: the original stopping rule — target N and/or funding window —
for the record.]

---

# Variables

## Manipulated variables

(1) Respiratory control: paced (exogenous, visual pacer) vs free (endogenous).
(2) Instructed change direction within paced epochs: Faster / Slower / No
Change (pacer period endpoint 0.65–0.80× or 1.20–1.35× baseline).
(3) Change salience: high (flat, then step at epoch midpoint) vs low (linear
ramp across the epoch); both reach the same endpoint, so salience is the
trajectory shape, not the change magnitude.

## Measured variables

Behavioural, per trial: reported direction (Faster / No Change / Slower),
confidence (1–7), arousal (1–7), valence (run-level; per-epoch in v1).
Physiological, per volume and per epoch: belt-derived breath period, rate,
RVT; pacer adherence (run-level correlation of measured vs instructed breath
period); end-tidal CO2 (PetCO2) per breath and per volume. Questionnaire:
MAIA-2 (37 items). Imaging: preprocessed BOLD series (fMRIPrep 24.1.1) and
its confound time series. Demographics and screener variables are not
analysed here beyond the registered MAIA covariate.

## Indices

(1) Awareness conditions: paced_aware (response == instructed direction),
paced_unaware (response == "No Change" on a change trial), paced_wrongdir
(directional, wrong direction; regressor of no interest), paced_nochange.
"No Change" responses are treated as low-awareness, never as correct
rejections (behavioural accuracy on No Change trials is at chance, 0.301).
(2) RVT regressor: per-breath amplitude/period at breath midpoints,
interpolated to TR, z-scored within run, convolved with the Birn (2008)
respiration response function (kernel normalised by ABSOLUTE area; signed-area
normalisation inverts biphasic kernels).
(3) PetCO2 regressor: per-breath end-tidal peaks interpolated to a fine grid,
convolved with the SPM canonical HRF at fixed lags {0, 10, 20, 30 s} entering
as a 4-column basis (the transit lag is bounded below at ~9 s but not
identifiable above in these data; a single fitted lag would be false
precision). The Birn RRF is never applied to PetCO2 (wrong sign: RVT and
end-tidal CO2 move oppositely).
(4) MAIA total: mean of the six subscales excluding Not-Worrying and
Not-Distracting (matching the IEAT), on the 0–5 metric.
(5) Pacer adherence: run-level Pearson r between measured breath period and
the instructed period sampled at each breath onset (per-epoch classification
is underpowered at ~3 breaths/epoch and is not used inferentially).

---

# Analysis Plan

## Statistical models

First level (per run, GLM, canonical HRF + temporal derivative, 6 mm FWHM
smoothing, MNI152NLin2009cAsym, no drift model beyond the cosine confounds).
Three models share one nuisance set and differ only in how paced epochs are
split:
Model A (awareness; H1, H3, H6): paced_aware, paced_unaware, paced_wrongdir,
paced_nochange, free, rating, subrating; no parametric modulator.
Model B (confidence; H2 only): all paced epochs as one regressor carrying the
mean-centred confidence modulator; free, rating, subrating as in A.
Model C (salience; H4): paced change epochs split into high vs low salience;
paced_nochange, free, rating, subrating as in A.
Arousal is not a registered test; if reported, it is fitted in a separate
model of Model B's form and labelled exploratory.
Nuisance (all models): 6 motion + 6 derivatives, framewise displacement, CSF,
white matter, cosine drift (128 s), one spike regressor per volume with
FD > 0.5 mm. Sensitivity models add (S2) RVT ⊗ RRF from the Biopac belt, (S2b)
the same with the Siemens PMU cushion as fallback source, or (S1) the 4-lag
PetCO2 basis, or both. Global signal is never a regressor. Design QC before freezing:
all 90 designs are full rank in every model.
Contrasts: C1 (H1) = free vs the unweighted mean of the paced conditions
(aware, unaware, nochange) present in that run; C2 (H3) = aware − unaware;
C3 (H2) = confidence modulator; C5 (H4) = high − low salience. Only
between-condition contrasts are interpreted; no condition-vs-rest contrast
is reported.
Runs are combined within participant by fixed effects (run as a random factor
is degenerate at 2 levels; runs are behaviourally exchangeable: accuracy
0.800 vs 0.754, p = .22). For C3, only runs in which the confidence modulator
varies enter (excludes sub-106 run-1 and sub-112 run-2).
Group level: one-sample t on participant contrast images for H1–H4; MAIA
total (demeaned) as a group covariate on C1 for H6 (two-sided); protocol
version as covariate in sensitivity S5. Registered descriptive diagnostic:
run-1 vs run-2 spatial correlation of contrast maps and cross-participant ICC
of ROI means for C1 and C3 (behavioural benchmark r ≈ 0.7); H3 is excluded
from this diagnostic (estimable in both runs for only 1 participant).

## Transformations

Confidence and arousal mean-centred within run over modelled trials.
RVT z-scored within run (belt gain is uncalibrated and hand-set per session;
raw amplitude would compare amplifier settings). PetCO2 series demeaned by
convolution; physiological confound columns demeaned within run with missing
volumes set to 0 (mean imputation); a physiological column with < 50%
coverage in a run is omitted for that run (whole-basis rule for PetCO2: the
4 lag columns enter together or not at all). Trailing events that extend past
the end of the acquisition are clipped to the scan end and counted (3 runs).
MAIA responses converted from the administered 1–6 coding to the canonical
0–5 metric. No other transformations.

## Inference criteria

Whole brain: cluster-forming p < .001 uncorrected, cluster-level FWE p < .05,
two-sided.
ROI: small-volume FWE p < .05, Bonferroni-corrected over the ROIs tested per
hypothesis (H2/H3/H4: left anterior insula, right anterior insula, dACC; ×3),
one-sided in the predicted direction.
H6: two-sided p < .05, small-volume corrected in a 10 mm sphere at MNI
(−4, +12, +36), with the Harvard–Oxford dACC mask as a secondary test.
No uncorrected result will be reported as a finding.
H5 equivalence: TOST (α = .05 per side) on the paired within-participant
difference between corrected and uncorrected ROI estimates of the H2 effect;
bounds dz = ±0.61 (PetCO2, H2 N = 25) and ±0.59 (Biopac RVT, H2 N = 27), the smallest
bounds supportable at 80% power. Fixed interpretation grid: TOST passes →
statistically equivalent; TOST fails but corrected effect remains significant
→ attenuated but present; both fail → indeterminate (not "abolished"; the
design cannot separate abolition from power loss at these bounds). The
corrected–uncorrected correlation is reported so the bound can be read in raw
units. The corresponding H3 comparison is descriptive.

## Data inclusion and exclusion

Applied before this registration, from non-imaging data plus fMRIPrep's FD
column only: (1) pilot/QA identifiers excluded; (2) sub-051 (no MRI acquired)
and sub-056 (run-1 truncated at 229 volumes; run-2 fails the motion gate)
excluded; (3) motion gate: run excluded if mean FD > 0.5 mm or > 25% of
volumes exceed FD 0.5 mm — fails 12 of 103 runs and removes 5 participants
entirely (median 81 spike volumes vs 3 in retained runs); a registered check
shows the gate does not select against pacer-adherent breathers (FD–adherence
r = −0.003); (4) a participant enters when at least one breathing run passes
these criteria — run count is not an inclusion criterion (four participants
contribute one run; sensitivity S6 removes them together); (5) per-contrast
estimability: ≥ 3 trials per cell per participant, pooled across runs; for
H2, participants whose confidence is constant across all paced trials
(sub-066, sub-070, sub-114) are excluded from H2 only. Wrong-direction trials
(median 1/participant) are modelled but never contrasted. No outlier
exclusion beyond these rules; no exclusions may be added after imaging
analysis begins.
Note on an earlier draft: sub-082 and sub-104 were at one stage excluded for
"shifted run labels"; that applied only to an obsolete Windows copy of the
data. In the dataset of record their runs map correctly, and both are
included (one gate-passing breathing run each).

## Missing data

A missing physiological regressor excludes a run only from the corresponding
sensitivity model, never from the base model. Physio gaps within a run:
mean-imputed (0 after demeaning) when coverage ≥ 50%, else the column is
omitted for that run. MAIA is complete for 46 of the 47 analysable
participants (sub-065 has none and does not enter H6). Participants
contributing one run (rather than two) enter with that run's estimate. No
imputation of behavioural responses; paced trials without a response (0 in
the current data) would be modelled as a separate no-interest regressor.

## Other planned analysis

S1: ±PetCO2 model comparison (N = 25). S2: ±RVT, Biopac belt (N = 30).
S2b (robustness): ±RVT with PMU fallback where the belt failed (N = 41; PMU
RVT tracks the global signal but agrees poorly with the belt, so it is not
confirmatory; one gate was revised after fallback metrics were seen —
disclosed in prereg §6). S3: CO2-available
vs unavailable participants on the uncorrected H2 contrast — registered as a
cohort check only: CO2 availability is nearly collinear with protocol version
(capnograph outage 2023-10-12 to 2024-01-18), so this cannot isolate CO2.
S4: H2 and H3 restricted to participants with verified pacer adherence on a
retained run (N = 26; H2 25, H3 14). S5: protocol version as group covariate.
S6: every hypothesis on two-run participants only (N = 43; H1 43, H3 23,
H2 41, H4 43, H6 42).
Exploratory (labelled as such regardless of outcome): brainstem/thalamic
effects; free-breathing epochs vs measured respiratory change (noise floor
sd_log = 0.399 comparable to the imposed signal); arousal modulation; MAIA
subscale dissociations; direction-specific (Faster vs Slower) effects;
ACC-seeded PPI with the dorsal attention network.

---

# Other

## Context and additional information

This study descends directly from the IEAT (Farb, Zuo & Price, 2023, eNeuro
10(6):ENEURO.0088-23.2023): its paced epochs correspond to the IEAT's Active
Matching condition and its free epochs to Active Interoception. H1 is a
replication with its direction fixed by that paper; H6 carries the IEAT's
MAIA-sensitive ACC region over to this task's pacing contrast as a two-sided
conceptual extension, because the IEAT's MAIA effect was on interoception vs
exteroception, a contrast this design does not contain. Key behavioural facts
constraining the hypotheses (already analysed, see Foreknowledge):
three-category identification is well above chance (0.525 vs 1/3); when a
direction is reported on a true change it is correct 0.778 of the time, while
No Change trials are at chance and "No Change" is the lowest-confidence
response — it functions as "I don't know", so no detection-stage d′ is
interpreted as sensitivity. PetCO2 is a partial MEDIATOR of the task (pacing
manipulates end-tidal CO2 by design), which is why H5 is framed as
attenuation with an equivalence criterion rather than as a purity check. The
attached full analysis plan (intero2023_imaging_prereg.md, v0.4) contains the
dated decision history; the analysis repository (git) will be linked as an
OSF component, and any deviation after freezing is logged in deviations.md
with date, rationale, and whether the affected result had been seen.
