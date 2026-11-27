# Results summary — registered imaging analyses (OSF swc7p, plan v0.8)

Run 2026-11-26 after the swc7p update was approved. Pipeline and analyst
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

**Reading.** The salience masks show no detectable response to a change in
pace, whether or not it was noticed (E3a, E6), and moderate evidence against
one (E6 BF01 5–6). There is moderate evidence against aware > unaware in the
left anterior insula and dACC. The right anterior insula is consistently the
best region but inconclusive (BF01 2.2). The CAREN anterior insula and dACC masks
are broad: they extend back to y = −20.5 and y = −44.5.

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
