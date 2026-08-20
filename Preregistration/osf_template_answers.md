# OSF Preregistration — field-by-field answers (copy-paste ready)

Template: OSF Preregistration (osf.io registration form). Each heading below is
one template field, in template order. Text under it is paste-ready. Items
marked [TO CONFIRM] need the author's input before freezing.

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
statistic has been computed. Hypotheses include direct replications of Farb,
Zuo & Price (2023, eNeuro): insula/opercular deactivation under exogenous
respiratory pacing, and MAIA-moderated sparing of anterior cingulate activity.

## Contributors

Norman A. S. Farb (sole author).

## License

[TO CONFIRM — recommend CC-By Attribution 4.0 International, matching the IEAT
paper's open-access licence.]

## Subject

Social and Behavioral Sciences → Psychology → Cognition and Perception;
Life Sciences → Neuroscience and Neurobiology → Cognitive Neuroscience.

## Tags

interoception; respiration; breath awareness; fMRI; anterior cingulate;
insula; PetCO2; RVT; metacognition; replication

---

## Research questions or hypotheses

H1 (replication). Exogenously paced breathing, relative to free breathing,
produces relative DEACTIVATION of the insula/operculum pathway and somatomotor
cortex, replicating the Active Matching result of Farb, Zuo & Price (2023,
eNeuro 10(6); "IEAT"). The volitional-control account predicts the opposite
(activation); a reversed result will be reported as a failed replication, not
reframed.

H2 (primary awareness test). Trial-by-trial confidence on change trials
correlates positively with BOLD activity in bilateral anterior insula and
dorsal anterior cingulate cortex (dACC). N = 43.

H3 (secondary awareness test). Change trials answered with the correct
direction ("aware") elicit greater activity than change trials answered "No
Change" ("unaware") in the same three regions of interest. N = 21;
ROI-restricted. Registered in advance: a null here alongside a positive H2
will be attributed to power, since only large effects (dz ≈ 1.0) are
detectable at this N.

H4. High-salience (step) changes elicit greater salience-network activity
(anterior insula, dACC) than low-salience (ramp) changes, consistent with the
behavioural finding that salience affects noticing (p = .043) but not
direction discrimination (p = .26). N = 41.

H5. Adding physiological regressors (RVT; PetCO2) attenuates but does not
abolish the H2 awareness effects. Adjudicated by a preregistered equivalence
criterion (see Inference criteria), not by threshold-crossing alone.

H6 (replication). Higher self-reported interoceptive sensibility (MAIA)
predicts greater (spared) ACC activity during paced breathing, tested
one-sided in a 10 mm sphere at the IEAT peak (MNI −4, +12, +36). N = 43.

## Foreknowledge of data or evidence

[Select the option closest to: "The data exist and I have accessed and
analyzed some of the data relevant to this registration."]

## Explanation of foreknowledge

Data collection is complete (September 2023 – June 2024). Everything EXCEPT
the functional imaging has been analysed before this registration: task logs
and event timing; all trial-level behaviour (three-category identification
accuracy 0.525 vs chance 1/3; ordinal sensitivity 0.698 probit units/step;
confidence type-2 AUC 0.661); respiration-belt processing including pacer
adherence; capnography/PetCO2 extraction; RVT; and all quality control and
exclusion decisions. The functional images have not been examined beyond
fMRIPrep preprocessing: no first-level model has been fitted, no contrast
image exists, no statistical map has been computed or viewed. Two narrow,
documented exceptions: fMRIPrep's global_signal column was used to verify the
physiological regressors carry the physiologically correct sign, and its
framewise_displacement column was used to apply the preregistered motion
exclusion. No task regressor has been related to any imaging measure.
Consequence, stated plainly: H2–H4 are informed by known behavioural results;
they are untested predictions about brain data but are not blind. H1 and H6
take their direction and coordinates from an independent published sample
(IEAT). The full analysis plan with dated decision history is attached
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
free-breathing epoch (~20 s) → rating window]. Change direction and salience
orders randomized within run. Two protocol versions exist (v1 n≈9 with
per-epoch valence ratings and longer runs; v2 the remainder); protocol version
is a registered group-level covariate (sensitivity S5). Acquisition: 3T
Siemens, EPI, TR = 2 s, 2 mm isotropic; PEPOLAR fieldmaps; fMRIPrep 24.1.1,
MNI152NLin2009cAsym. Concurrent physiological recording: Biopac respiration
belt and capnograph (2 kHz, downsampled to 100 Hz), scanner-triggered;
Siemens PMU respiration (400 Hz) as a validation stream.

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
runs). After preregistered exclusions (below), the analysable samples are:
H1/H2/H6: N = 43; H4: N = 41; H3: N = 21; sensitivity S1 (PetCO2): N = 22;
S2 (RVT): N = 26. Trials per participant (median): 7 aware, 2 unaware,
1 wrong-direction, 2 no-change, 14 free.

## Sample size rationale

The sample was collected before this registration; no further collection will
occur, so size is characterized rather than justified. Smallest detectable
effect at 80% power (one-sample t): N = 43 → d = 0.44 at α = .05 (ROI),
d = 0.67 at cluster-forming p < .001; N = 21 → d = 0.63 / 1.01. H6 at N = 43
has ~67% one-sided power at the IEAT's observed effect (r ≈ 0.32); this is
stated so a marginal replication is read against it. Equivalence bounds for
H5 are set by what these Ns support (see Inference criteria).

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
smoothing, MNI152NLin2009cAsym): event regressors paced_aware, paced_unaware,
paced_wrongdir, paced_nochange, free, rating, subrating; parametric modulators
confidence (primary) and arousal (secondary) on paced epochs, mean-centred
over modelled trials; nuisance: 6 motion + 6 derivatives, framewise
displacement, CSF, white matter, cosine drift (128 s), one spike regressor
per volume with FD > 0.5 mm. Sensitivity models add (S2) RVT ⊗ RRF or
(S1) the 4-lag PetCO2 basis, or both. Global signal is never a regressor.
Runs are combined within participant by fixed effects (run as a random factor
is degenerate at 2 levels; runs are behaviourally exchangeable: accuracy
0.800 vs 0.754, p = .22).
Group level: one-sample t on participant contrast images for H1–H4; MAIA
total as a group covariate for H6 (one-sided); protocol version as covariate
in sensitivity S5. Registered descriptive diagnostic: run-1 vs run-2 spatial
correlation of contrast maps and cross-participant ICC of ROI means for the
two well-powered contrasts (behavioural benchmark r ≈ 0.7); H3 is excluded
from this diagnostic (estimable in both runs for only 1 participant).

## Transformations

Confidence and arousal mean-centred within participant over modelled trials.
RVT z-scored within run (belt gain is uncalibrated and hand-set per session;
raw amplitude would compare amplifier settings). PetCO2 series demeaned by
convolution; physiological confound columns demeaned within run with missing
volumes set to 0 (mean imputation); a physiological column with < 50%
coverage in a run is omitted for that run (whole-basis rule for PetCO2: the
4 lag columns enter together or not at all). MAIA responses converted from
the administered 1–6 coding to the canonical 0–5 metric. No other
transformations.

## Inference criteria

Whole brain: cluster-forming p < .001 uncorrected, cluster-level FWE p < .05.
ROI: small-volume FWE p < .05, Bonferroni-corrected over the ROIs tested per
hypothesis (H2/H3: left anterior insula, right anterior insula, dACC; ×3).
H6: one-sided p < .05, small-volume corrected in a 10 mm sphere at MNI
(−4, +12, +36), with the Harvard–Oxford dACC mask as a secondary test.
No uncorrected result will be reported as a finding.
H5 equivalence: TOST (α = .05 per side) on the paired within-participant
difference between corrected and uncorrected ROI estimates of the H2 effect;
bounds dz = ±0.65 (PetCO2, N = 22) and ±0.59 (RVT, N = 26), the smallest
bounds supportable at ≥ 80% power. Fixed interpretation grid: TOST passes →
statistically equivalent; TOST fails but corrected effect remains significant
→ attenuated but present; both fail → indeterminate (not "abolished"; the
design cannot separate abolition from power loss at these bounds). The
corrected–uncorrected correlation is reported so the bound can be read in raw
units.

## Data inclusion and exclusion

Applied before this registration, from non-imaging data plus fMRIPrep's FD
column only: (1) pilot/QA identifiers excluded; (2) sub-082 and sub-104
excluded for shifted run labels in the derivatives (task-breath1 file holds
172 volumes where the slot expects 242; all 242 other runs reconcile exactly);
(3) motion gate: run excluded if mean FD > 0.5 mm or > 25% of volumes exceed
FD 0.5 mm — removes 11/95 runs and 5 participants (median 81 spike volumes vs
3 in retained runs); a registered check shows the gate does not select
against pacer-adherent breathers (FD–adherence r = −0.003); (4) per-contrast
estimability: ≥ 3 trials per cell per participant. Wrong-direction trials
(median 1/participant) are modelled but never contrasted. No outlier
exclusion beyond these rules; no exclusions may be added after imaging
analysis begins.

## Missing data

A missing physiological regressor excludes a run only from the corresponding
sensitivity model, never from the base model. Physio gaps within a run:
mean-imputed (0 after demeaning) when coverage ≥ 50%, else the column is
omitted for that run. MAIA is complete for all 43 analysable participants.
Participants contributing one run (rather than two) enter with that run's
fixed-effects estimate. No imputation of behavioural responses; paced trials
without a response (0 in the current data) would be modelled as a separate
no-interest regressor.

## Other planned analysis

S1: ±PetCO2 model comparison (N = 22). S2: ±RVT (N = 26). S3: CO2-available
vs unavailable participants on the uncorrected H2 contrast — registered as a
cohort check only: CO2 availability is nearly collinear with protocol version
(capnograph outage 2023-10-12 to 2024-01-18), so this cannot isolate CO2.
S4: H2 restricted to participants with verified pacer adherence (N = 24
pre-gate). S5: protocol version as group covariate.
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
Matching condition and its free epochs to Active Interoception; H1 and H6 are
replications with direction and coordinates fixed by that paper. Key
behavioural facts constraining the hypotheses (already analysed, see
Foreknowledge): three-category identification is well above chance (0.525 vs
1/3) with direction accuracy 0.778 on true changes, while No Change trials
are at chance and "No Change" is the lowest-confidence response — it
functions as "I don't know", so no detection-stage d′ is interpreted as
sensitivity. PetCO2 is a partial MEDIATOR of the task (pacing manipulates
end-tidal CO2 by design), which is why H5 is framed as attenuation with an
equivalence criterion rather than as a purity check. The attached full
analysis plan (intero2023_imaging_prereg.md) contains the dated decision
history; the analysis repository (git) will be linked as an OSF component,
and any deviation after freezing is logged in deviations.md with date,
rationale, and whether the affected result had been seen.
