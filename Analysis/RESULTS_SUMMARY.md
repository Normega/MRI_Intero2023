# Results summary — registered imaging analyses (OSF swc7p, plan v0.8)

Run 2026-11-26 after the swc7p update was approved. The localizer (OSF p5ztc) summary is at the end of this file. Pipeline and analyst
choices: `Analysis/glm_implementation_notes.md`, committed before any group
statistic existed. Full tables: `Analysis/results/`. Corrected p-values are
FWE with the registered Bonferroni. Permutation p has a floor of 1e-4
(10,000 permutations), so ×3 gives 3e-4. **No uncorrected result is a finding.**

## Verdicts

| Hypothesis | N | Verdict | Key numbers |
|---|---|---|---|
| **H1** paced < free, somatomotor masks | 47 | **Supported** | Every mask p ≤ .0003 (×3). Peak t −8.9 (posterior insula), −9.4 (central operculum), −8.1 (sensorimotor). Mask-mean dz −0.92 to −1.13. The opposite tail is null everywhere. |
| **H3** aware > unaware (authoritative) | 45 | **Not supported** | Best mask right anterior insula, z 3.07, p .23 (×3). REML z 1.51. Whole brain null in both tails. |
| **H2** confidence, salience masks (supportive) | 44 | **Not supported** | Best mask right anterior insula, t 3.39, p .17 (×3). Whole brain null. |
| **H4** step > ramp, salience masks | 43 | **Not supported** | Best mask left anterior insula, t 2.98, p .40 (×3). Whole brain null. |
| **H5** physiological correction | — | **Not tested** | Neither H3 nor H2 is supported, so under the registered rule both comparisons are descriptive (below). |
| **H6** MAIA × paced − free, IEAT ACC sphere, two-sided | 46 | **Supported** | t +3.28 at (−2.5, 13.5, 45.5), p_FWE .036, 4 voxels. Higher MAIA means less deactivation under pacing. Whole brain: one positive cluster, 687 voxels, peak (3.5, 35.5, 3.5), pregenual ACC (Harvard–Oxford 86%), p .005. |
| **H7a** coupling paced vs free | 47 | **Not supported** | SAL–SMN dz −0.10, SAL–ECN +0.13. Holm p .72 for both. |
| **H7b** coupling aware vs unaware | 23 | **Not supported** | SAL–SMN dz −0.02, SAL–ECN +0.12. Holm p 1.0 for both. |

**Whole brain, H1.**
- Four clusters survive. The largest is negative: 87,588 voxels, peak t −10.4 at (−52.5, −54.5, 25.5). It is widespread deactivation, as in the IEAT.
- One cluster is positive: 3,381 voxels, peak t 8.1 at (45.5, −62.5, 7.5), right lateral occipitotemporal cortex. Harvard–Oxford labels it lateral occipital / middle temporal cortex, near MT, which fits the moving visual pacer. The negative peak is in left angular / supramarginal gyrus.

## Sensitivity analyses (best mask of the hypothesis's registered set; ×3 or the registered k)

| | H1 | H3 | H2 | H4 | H6 (sphere) |
|---|---|---|---|---|---|
| Primary | .0003 | .23 | .17 | .40 | .036 |
| S1 base / +PetCO2 (N 25; H4 22; H6 24) | .0003 / .0003 | .39 / .68 | **.028** / .07 | **.018** / **.007** | .55 / .63 |
| S2 base / +RVT, Biopac (N 30; H3 28; H2 27; H4 26; H6 29) | .0003 / .0003 | .89 / 1.0 | **.048** / **.011** | .76 / 1.0 | **.036** / **.028** |
| S2b base / +RVT incl. PMU (N 41; H3 39; H2 38; H4 37; H6 40) | .0003 / .0003 | .46 / .69 | .17 / .13 | .49 / .84 | **.012** / **.022** |
| +both, corrected model (N 21; H4 18; H6 20) | .0003 | 1.0 | **.019** | **.007** | .29 |
| S4 adherence (H3 26, H2 25) | — | .83 | .19 | — | — |
| S5 protocol covariate | — | .33 | .16 | — | — |
| S6 two-run only (N 43; H3 42; H2 41; H6 42) | .0003 | .20 | .12 | .40 | **.018** |
| S7 3-trial rule, FLAME / OLS (N 23) | — | .17 / .23 | — | — | — |
| H3 per-run form (v0.6) | — | .12 | — | — | — |
| Harvard–Oxford masks | .0003 | .07 | .13 | .66 | — |
| S3 CO2-available vs not, H2 (two-sided) | — | — | .16 | — | — |


**Reading, as registered.**
- H1 survives every sensitivity analysis and the atlas swap.
- H6 survives S2, S2b, S6 and RVT correction. It does not reach significance in the PetCO2 subset (S1, N = 24). The base model in that subset is null as well, so the loss comes from the subset, not the correction.
- H2 and H4 reach p < .05 in some smaller subsets (S1, S2, +both) but not in the full sample. This divergence is reported as dependence on the subset. It is **not** used to select a result.
- H3 is null in every form.

## H5 comparisons (descriptive only — H5 not tested)

The corrected–uncorrected correlation is .85–.98 throughout.

| Target | Correction | Subset | N | dz of change | TOST p | Grid cell (descriptive) |
|---|---|---|---|---|---|---|
| H3 | PetCO2 | S1 | 25 | −0.11 to +0.07 | .004–.010 | equivalent |
| H3 | RVT (Biopac) | S2 | 28 | −0.35 to −0.42 | .12–.20 | indeterminate |
| H3 | RVT incl. PMU | S2b | 39 | −0.42 to −0.50 | .34–.52 | indeterminate |
| H2 | PetCO2 | S1 | 25 | −0.06 to 0.00 | .003–.005 | equivalent |
| H2 | RVT (Biopac) | S2 | 27 | −0.02 to +0.02 | .003 | equivalent |
| H2 | RVT incl. PMU | S2b | 38 | −0.03 to +0.05 | .003–.005 | equivalent |

RVT correction shifts the aware − unaware estimate downward (dz about −0.4). It does not shift the confidence effect.

## Reliability diagnostic (descriptive; behavioural benchmark r ≈ 0.7)

| Contrast | Within-participant spatial r, run 1 vs run 2 (median [IQR]) | ROI-mean ICC(3,1) across participants |
|---|---|---|
| C1 paced − free (N 43) | 0.34 [0.25, 0.45] | 0.32–0.45 (H1 masks) |
| C3 confidence (N 39) | 0.00 [−0.05, 0.21] | 0.22–0.35 (H2 masks) |

## H7 checks

- **Global-signal gate:** not reached, because neither part is supported.
- **Five-network-mean partial:** H7a SAL–SMN coupling is *lower* under pacing (dz −0.39, Holm p .02). This is opposite to the prediction, and it is a check, not the test.
- **Amplitude:** pacing lowers temporal SD in all five networks (dz −0.55 to −1.21).
- **Co-activation:** network-mean paced − free activation is negative in SAL, SMN and ECN (dz −0.71 to −1.00), consistent with H1.
- **Specifications:** all 8 H7a and all 4 H7b specifications, and the +RVT and +PetCO2 versions, are null.
- **All 10 pairs (exploratory):** none survives FDR. The smallest is H7a SMN–VIS, p_FDR .055.

## Exploratory (labelled as such; cannot support or overturn any hypothesis)

All are logged in `Preregistration/deviations.md` ("After the update"), each
added after the registered results were seen. Atlas labels are Harvard–Oxford
cortical/subcortical and the FSL cerebellar atlas, and are descriptive.

### Awareness

| | Question | N | Result |
|---|---|---|---|
| E1 | Aware vs unaware, H3 masks, both tails (FLAME, GRF, ×6) | 45 | Null. Best: right anterior insula aware > unaware, z 3.07, p .46. |
| Model A-g | Graded awareness (confidence modulator), H3 masks, ×3 | 45 / 38 | Null. Right anterior insula is best each time: signed p .15, rank p .12, aware-only p .086. |
| E3a | Change (aware, unaware, wrong-direction) vs No Change, whole brain, two-sided | 47 | Null. Largest cluster 23 voxels, p_FWE .78. |
| E3b | Aware vs unaware inside E3a's surviving clusters | — | Not run: no cluster survived E3a, so the mask is empty. |
| E4 | Aware vs unaware, whole-brain TFCE, two-sided (5,000 permutations) | 45 | Null. Smallest p: aware > unaware .56, aware < unaware .96. No voxel p < .05. |
| E5 | JZS Bayes factors, H3 masks (r = 0.707) | 45 | Aware − unaware BF01: left anterior insula 5.1, right 2.2, dACC 5.8 (directional BF0+ 3.4, 1.2, 8.1). Confidence (H2, N 44) BF01: 2.6, 2.0, 5.5. |
| E6 | Change vs No Change, H3 masks, two-sided voxel max-t, ×3 | 47 | Null: p_FWE×3 = 1.0 in all three; no voxel p < .05. Mask-mean dz +0.02, +0.10, −0.06; BF01 (r 0.707) 6.2, 5.0, 5.8. |
| E7 | Transient response to the pacer change: FIR (15 × 2 s bins), joint model, H3 masks, two-sided voxel max-t, ×3 | 46 | Null. Primary (step − ramp, [12, 18) s minus [6, 12) s): best dACC t 3.40, p .54; insulae p 1.0. Mask-mean dz +0.05, +0.01, +0.06; BF01 6.0, 6.2, 5.9. Secondary step − No Change [12, 18) s: best right anterior insula t 3.43, p .31, dz +0.14; ramp − No Change [6, 12) s: p 1.0 everywhere. BF01 4.2–6.2 for every test. Whole brain (primary): largest cluster 26 voxels, p_FWE .67. |
| E8 | Report period: the rating window after paced change epochs, split by the report, joint model, H3 masks, two-sided voxel max-t, ×3 | 36 / 32 / 44 | Primary (reported change − reported No Change, i.e. hit + wrong-direction − miss): null. p_FWE×3 1.0, 1.0, .68 (best dACC t −3.09). Mask-mean dz −0.11, −0.17, −0.20 (change < No Change); BF01 4.6, 3.4, 2.9. Whole brain: largest cluster 13 voxels, p_FWE .93. Robustness (28 two-run participants): run-1/run-2 sign agreement 54%, 57%, 68%; r −0.03, −0.02, +0.51. Secondary confidence modulator (N 44): null, p 1.0 everywhere, dz +0.25, +0.15, +0.04, BF01 1.7, 3.8, 5.9. Secondary Faster − Slower (**motor / key-mapping diagnostic**, N 32): Slower > Faster in left anterior insula (t −5.21, p ×3 .004), right anterior insula (t −4.68, p .017), dACC (t −4.20, p .071). dz −0.48 to −0.53, BF10 4.3–7.9. Whole brain: four negative clusters survive: right frontal pole (934 voxels, p .004), medial superior frontal / pre-SMA (487, p .010), pregenual ACC (198, p .037), left precuneus (194, p .038). There is no left-sensorimotor (finger) cluster. |
| E10 | Is E8's Faster − Slower rating effect about the report or about the stimulus? Same rating split as E8; paced change epochs split by stimulus direction × outcome (joint model, H3 masks, two-sided voxel max-t, ×3) | 32 / 22 / 47 | Primary, reported Faster − Slower with the epochs direction-controlled (E10rep, N 32): **survives, almost unchanged**. Slower > Faster in left anterior insula (t −4.40, p ×3 .034), right anterior insula (t −4.75, p .015), dACC (t −3.87, p .14). dz −0.42, −0.50, −0.38 (E8FS in the same 32: −0.48, −0.53, −0.49; ratio of mask means 1.00, 1.05, 0.89; r across participants .97). BF10 2.2, 5.4, 1.4. Whole brain: three negative clusters survive, all right frontal: frontal pole (699 voxels, p .005), frontal pole / inferior frontal pars triangularis (399, p .011), inferior frontal pars triangularis / opercularis (175, p .041). E8's pre-SMA, pregenual ACC and precuneus clusters do not survive (the nearest are paracingulate / medial frontal 104 voxels, p .08, and SMA 77, p .12). Secondary, epoch stimulus direction × noticed ((F hit − S hit) − (F miss − S miss), N 22, ≥ 1 per cell): null, p 1.0 everywhere, dz −0.04, −0.22, +0.01, BF01 4.4, 2.8, 4.5. Diagnostic, epoch stimulus Slower − Faster (N 47): null, p ×3 .96, .65, 1.0, dz −0.11, −0.12, +0.01, BF01 4.9, 4.6, 6.3; whole brain null (largest 160 voxels, p .057). Wrong-direction dissociation (descriptive; 79 epochs; a voxelwise test was infeasible, only 10 participants have both wrong-direction cells): the pooled WLS (46 participants) splits the effect into report F − S −0.085, −0.088, −0.067 and stimulus F − S −0.052, −0.076, −0.058. Every 95% CI includes 0. Rating SF − FS in the 10 participants: −0.02, −0.00, +0.02 (5–6 of 10 negative). |

**Reading.** The salience masks show no detectable response to a change in
pace, whether or not it was noticed (E3a, E6), and moderate evidence against
one (E6 BF01 5–6). The same holds for a brief change-locked response (E7). Step,
ramp and No Change epochs have the same FIR time course in all three masks: flat
to about 10 s, then a fall of about 0.3–0.4% that lasts to the epoch end. Step
sits slightly above ramp by a constant 0.03–0.08% in every bin, including the
bins before any change. The one sharp difference, at 24–26 s, follows epoch
length (No Change epochs last 20.2 s, change epochs 17.5 or 22.7 s) and is
outside the logged windows. The report period does not carry the decision
either (E8): the rating window after a reported change and after a miss do not
differ in the masks or the whole brain (BF01 2.9–4.6, weaker than E6/E7 because
N is 36), and the estimate is not stable across runs. The Faster − Slower
diagnostic does differ (Slower > Faster in both anterior insulae and medial
frontal cortex), but it is not a finger pattern. About 80% of each of its cells
share the preceding stimulus direction, and Slower epochs are longer (23.0 s vs
17.5 s) and breathed more slowly, which Model A does not split. E8 therefore read
it as a probable carry-over from the preceding paced epoch. **E10 tested that and does
not support it.** When the epochs are split by stimulus direction, the rating-period
Slower > Faster effect is essentially unchanged in both anterior insulae (mask means 0.89–1.05
of E8's), and the right frontal pole / inferior frontal clusters remain. The epoch itself
shows no stimulus-direction effect (BF01 4.6–6.3), so there is no epoch response for the
rating effect to have leaked from. E10 cannot settle whether the effect follows the report or
is an after-effect of the stimulus that falls in the rating window. Report and stimulus agree
in about 80% of epochs, and the 79 wrong-direction epochs split the effect about evenly
between the two, with every interval including zero. The best current reading is a
rating-period difference after Slower vs Faster trials that the epoch model does not explain,
of undetermined source (report or stimulus). There is no epoch-level direction × noticed
interaction (E10, BF01 2.8–4.5). There is moderate evidence against aware > unaware in the
left anterior insula and dACC. The right anterior insula is consistently the
best region but inconclusive (BF01 2.2). The CAREN anterior insula and dACC masks
are broad: they extend back to y = −20.5 and y = −44.5.

### Interoceptive trait and confidence

E9 asks whether the H6 effect (MAIA × paced − free in the registered IEAT ACC
sphere) relates to confidence. Behaviourally, MAIA tracks mean paced confidence
(Spearman .23, Pearson .32, n 46) but not direction accuracy (Spearman −.09).
Confidence is the mean over each participant's paced epochs (7–14). Everything
below is exploratory.

| | Question | N | Result |
|---|---|---|---|
| E9a (primary) | C1 (paced − free) in the IEAT ACC sphere vs mean confidence, Spearman, permutation | 47 | Null. ρ −.02, p .92 (Pearson r +.03). JZS BF10 0.12 (s 1; 0.16 at s 0.707, 0.29 at s 0.354), so BF01 about 9: moderate-to-strong evidence of no linear relation. |
| E9b | As E9a, partial on MAIA (ranks) | 46 | Null. Partial ρ −.13, p .40 (Holm ×4: 1.0). Chain, descriptive: MAIA–C1 ρ +.37, MAIA–confidence ρ +.23, C1–confidence ρ −.03. |
| E9c | MAIA × confidence-modulator maps in the sphere, two-sided voxel max-t FWE | 43 / 43 | Null. C3 (Model B, epoch): peak t +1.50, p_FWE .90. E8conf (rating window): peak t +2.41, p_FWE .28 (×2: .56). Holm ×4: 1.0 for both. No voxel p < .05. |
| E9a-LOO | C1 in a leave-one-out H6 cluster vs confidence | 47 | Feasible: no fold empty, and every fold cluster contains the H6 peak (43–1,182 voxels, median 670; Dice with the full cluster .12–1.00, median .94). Null: ρ +.06, p .71 (Holm 1.0). MAIA–C1(LOO) ρ +.47. |
| Descriptive | C1 in the sphere vs type-2 AUC | 31 | ρ −.36, bootstrap 95% CI [−.67, −.01]. Not a test. The sign is opposite to the MAIA chain, the N is small and the AUCs rest on few trials. |

Reliability (descriptive): the sphere C1 has run-1 vs run-2 ICC(3,1) .49 (N 43).
Mean confidence has split-half reliability .91 (Spearman–Brown). Between-person
SD is 0.80 on the 1–7 scale, and 51% of paced ratings are 6 or 7.

**Reading.** The ACC response that tracks MAIA does not track confidence. This
holds between people (E9a, with BF01 about 9, and E9a-LOO), after removing MAIA
(E9b), and within people: MAIA does not moderate how ACC activity follows
trial-by-trial confidence, either during the epoch or during the rating (E9c).
The MAIA–C1 and MAIA–confidence relations are both present in this sample, but
they do not share variance: C1 and confidence are uncorrelated (ρ −.03). On
these data, H6 and the MAIA–confidence link look like separate correlates of
the trait, not a single pathway. Caveats:
- With N 47, only correlations of about |ρ| ≥ .29 are detectable at 50% power.
- The sphere C1 has reliability of about .5, which attenuates any true relation by about 0.7.
- Confidence has a restricted range.
- The behavioural MAIA–confidence correlation itself is modest (Spearman .23, p .12).

### Physiological regressors (E2: whole brain, two-sided, cluster FWE)

| Regressor | Sample | Result |
|---|---|---|
| RVT ⊗ RRF | S2b (N 41) | One positive cluster, 65,613 voxels (about 30% of the mask), peak t 10.4 in right cerebellum Crus II. |
| RVT ⊗ RRF | S2 Biopac (N 30) | 12 positive clusters. The largest are: right cerebellum X/VIII (11,595 voxels, peak t 8.5); left hippocampus / medial temporal (3,913); right planum polare / posterior insula (3,181); precuneus / paracentral (598); SMA / superior frontal (570); right thalamus (239). |
| Net PetCO2 (sum of 4 lags) | S1 (N 25) | Two small clusters: positive left occipital pole (179 voxels, p .012); negative lingual gyrus / vermis (161 voxels, p .015). |

**Reading.** RVT has the expected widespread positive BOLD relation, strongest
in cerebellum, deep grey matter and peri-vascular territory. Correcting for it
barely changes H1 or H6, so neither is a breathing-volume artifact. The net
CO2 term carries little unique spatial variance.

### Physiological correction of the supported hypotheses

H1 is unchanged by every correction: mask dz shifts by 10% or less, and p stays
.0003. H6 survives RVT correction (S2 .028, S2b .022). In S1 it is null in the
base model as well as the corrected one, so the loss there comes from the
subset (see Sensitivity).

### Nuisance model: the CSF regressor is task-locked (E11, exploratory)

Logged before running (`Preregistration/deviations.md`, E11). Script:
`Scripts/intero2023_explore_csf_breathing.py`. Table:
`Analysis/results/EXPLORATORY_E11_csf_paced_free.csv`.

The fMRIPrep `csf` and `white_matter` columns are in every registered model. Each
was regressed (as percent of its run mean) on the run's Model A design without
them. The paced − free (C1) contrast, averaged over runs (N 47):

| Tissue signal | C1 mean (PSC) | dz | % < 0 | sign-flip p |
|---|---|---|---|---|
| CSF | −0.062 | −0.45 | 70% | .0025 |
| White matter | −0.032 | −0.49 | 68% | .0015 |

**Reading.**
- **Both signals are lower during paced than free epochs.** So part of the
  paced − free variance is shared with the nuisance regressors, and the
  registered model removes it from H1's estimate.
- **The localizer shows the same thing.** In the localizer, the same CSF regressor falls in every task block, and dropping
  it roughly doubles the Match − Intero (L2) effect (whole-brain dz −0.30 →
  −0.61; see `LOCALIZER_RESULTS.md`, "localizer model without the CSF and
  white-matter regressors").
- **H1 is therefore a conservative estimate.** It was supported, with dz about −1,
  with this variance removed.
- **What is not settled.** Whether the removed part is neural (the
  global-grey-matter / CSF anti-coupling) or vascular / respiratory cannot be
  decided from these data.
- **What was not done.** No breathing-run refit without the tissue regressors has been run. The
  null awareness tests (H3, H2, E1–E10) contrast conditions within paced
  epochs, so a paced-vs-free-locked nuisance signal should affect them far less.

# Localizer (OSF p5ztc, plan v1.0 at d29136b): summary

Full results: `Analysis/LOCALIZER_RESULTS.md`. Deviations and every exploratory
analysis: `Preregistration/localizer_deviations.md`. Figures:
`Figures/localizer/`.

**Conditions.** Rest is the implicit baseline. The four task conditions are:
- **Intero:** track one's own breath with presses;
- **Extero:** track a visual circle;
- **Match:** breathe with the circle, the paced-breathing analogue;
- **Self:** trait-word judgements.

## Registered verdicts

| Test | N | Verdict | Key numbers |
|---|---|---|---|
| **L1** Intero − Extero ≠ 0, salience masks (two-sided, ×3) | 48 | **Supported in all three masks: Intero < Extero** | AInsL z −4.95, p .0007, dz −0.65; AInsR z −4.99, p .0005, dz −0.40; dACC z −5.30, p .0001, dz −0.55. |
| **L2** Match − Intero < 0, H1 masks (one-sided, ×3) | 48 | **Supported at peak level in COper and SM, not PIns; mask means null** | COper p .016, SM p .0495, 1 voxel each; mask dz −0.04 to −0.13, BF10 0.16–0.23. |
| **L3** MAIA × (Intero − Extero), IEAT ACC sphere | 46 | **Not supported** | r −.02, p .79; BF01 6.2. |
| **L4** reliability gate for participant-specific fROIs | 44 | **Not met in any parcel** | Cross-run C1 is negative everywhere (dz −0.17 to −0.53). The group fallback becomes primary in AInsR and the sphere. |
| **T1** aware > unaware in the C1 fROIs (×3) | 40 | **Not supported** | dz +0.14 (AInsL), +0.33 (AInsR, fallback), −0.04 (dACC); BF01 4.1 / 0.84 / 5.7. |
| **T3b** MAIA × paced − free in the sphere fROI (one-sided) | 40 | **Supported at the criterion; weak, not fROI-specific** | r .26, permutation p .042, BF10 0.57; the rest of the sphere is stronger (r .38). |
| T3 (secondary) paced − free in the C5 fROIs | 41 | Negative in all five, no more than the parcel | dz −0.58 to −0.89. |

**The direction of L1 is the expected one.** Attending to the breath rather than
an external target produces widespread cortical deactivation (Farb et al.,
2023, *eNeuro*). L1 was registered two-sided. The result replicates that
pattern; it is not a reversal. Its suprathreshold voxels are part of one
whole-brain negative cluster (about 72,500 voxels), so the effect is not specific to
the salience masks.

## Robustness

- **Physiological correction.** RVT⊗RRF and PetCO2 keep 85–107% of each
  contrast's whole-brain effect (Intero 94% / 99%, C1 88% / 91%).
- **The CSF regressor.** The CSF regressor is task-locked in the localizer as well. Dropping the tissue regressors has two effects:
  - **C1:** dz roughly halves. The tissue regressors mainly remove noise here, so the registered model is the more sensitive one for L1.
  - **Match contrasts:** they become much stronger (C5 whole-brain dz −0.30 → −0.61). The weak L2 result is partly due to the nuisance model.

  No verdict changes.
- **The periventricular positive band.** It sits on the borders of the lateral ventricles in the task-vs-Rest maps. It survives physiological correction, and it survives removing the CSF regressor.

## Exploratory (stage C; labelled, no verdict changes)

- **Behaviour (C1).**
  - Participants breathed with the pacer in both tasks: 38/40 in localizer Match, 33/33 in breathing-run change detection.
  - Breath tracking was accurate (PLV .79 vs chance .26) but only moderately reliable (ICC .44).
  - Tracking accuracy and adherence relate to no brain measure after FDR. Adherence does not moderate H1.
- **Connectivity (C2).**
  - Under Intero − Extero, salience–DMN coupling is less anticorrelated (network FC dz +0.69; background FC +0.48). This survives partialling out the global signal.
  - The salience–somatomotor change is a global-signal effect.
  - In the breathing runs, left–right anterior insula coupling is lower under pacing (q .022).
- **Pattern analyses (C3).**
  - **Decoding.** Intero vs Extero is decodable across participants (balanced accuracy .72 whole brain; .70–.71 salience, sensorimotor and DMN; .61 insula). This holds after the map mean is removed.
  - **Match.** Match maps classify between the two training classes. In RSA, Match is the condition closest to Intero in 7 of 9 regions.
  - **Breathing-run epochs (T4).** These express the localizer C1 pattern more when paced than free (+0.07). The expression does not track awareness (−0.017 [−0.047, 0.013], p .28) or confidence (p .86).
  - **No posterior → anterior insula gradient** (dz −0.07).
  - **Respiratory phase.** Inhale vs exhale responses look like breathing-depth artefact.

## Synthesis

Attending to the breath produces a robust, distributed cortical deactivation:
- It is decodable across people.
- It is shared with paced breath-matching.
- It comes with reduced salience–DMN anticorrelation.

The localizer does not yield reliable participant-specific "interoceptive"
regions. Neither its regions nor its multivariate pattern carries the awareness
of breathing changes. The breathing-run awareness null (BF01 ≈ 3–6) holds
when the measures are defined by the localizer.
