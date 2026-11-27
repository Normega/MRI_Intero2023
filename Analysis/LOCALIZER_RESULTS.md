# Results summary: registered localizer analyses, stage A (OSF p5ztc)

Run 2026-11-27, after OSF p5ztc was approved. Plan of record:
`Preregistration/intero2023_localizer_prereg.md` v1.0 (frozen at d29136b).
Pipeline and analyst choices: `Analysis/localizer_glm_implementation_notes.md`,
committed (400f1c8, c36230a) before any localizer statistic existed.

- Full tables: `Analysis/results/LOC_*`.
- Corrected p-values are FWE with the registered Bonferroni: FLAME 1 with a GRF
  voxel peak p within the mask (two-sided ×2 where registered). Permutation p has
  a floor of 1e-4.
- Effects are in percent signal change. Rest is the implicit baseline.
- **No uncorrected result is a finding.**

**Scope.** Stage A covers:
- localizer first levels (54 participants, 102 runs);
- L1, L2 and L3;
- the L1 breathing-robustness co-report;
- the registered sensitivity analyses;
- descriptive whole-brain maps.

Participant-specific ROIs, L4, T1, T3b and every exploratory analysis are later
stages, and none has been run.

## Verification (design only, before any statistic)
- **Nuisance sheets.** Built by the breathing-sheet rule. Rebuilding all 103
  committed breathing sheets with the same code reproduces them exactly
  (`Analysis/localizer_sheet_validation.csv`).
- **Inputs, all 102 runs** (`Analysis/localizer_glm_verification.csv`):
  - files exist;
  - header volumes = fMRIPrep rows = sheet rows;
  - sheet motion, FD and spike columns match fMRIPrep;
  - events fit after clipping (6 v1 runs clipped by 4–28 s);
  - every per-run design is full rank in every variant;
  - so are all 54 joint designs.
- **End to end, sub-067 run-1.** `fit()` reproduces nilearn `FirstLevelModel`
  (identical design matrix; effect, variance and t max |Δ| < 4e-7;
  `Analysis/localizer_glm_e2e_check.json`).

## Verdicts

| Test | N | Verdict | Key numbers |
|---|---|---|---|
| **L1** Intero − Extero ≠ 0, salience masks (two-sided, ×3) | 48 | **Supported, in all three masks; the direction is Intero < Extero** | AInsL: peak z −4.95 at (−34.5, 27.5, 3.5), p .0007, mask-mean dz −0.65, BF10 503. AInsR: z −4.99 at (39.5, 13.5, −10.5), p .0005, dz −0.40, BF10 5.0. dACC: z −5.30 at (−6.5, 13.5, 39.5), p .0001, dz −0.55, BF10 68. Permutation check p .0024 / .0009 / .0003. |
| **L2** Match − Intero < 0, H1 masks (one-sided, ×3) | 48 | **Supported at peak level in COper and SM; not in PIns. The mask means are null.** | COper: peak z −3.96 at (53.5, 9.5, −2.5), p .016, 1 voxel. SM: z −4.22 at (55.5, −26.5, 47.5), p .0495, 1 voxel. PIns: z −2.80, p .20. Permutation check: .013 / .013 / .53. Mask-mean dz −0.04 (PIns), −0.11 (COper), −0.13 (SM); BF10 0.16–0.23, BF₋₀ 0.19–0.38. |
| **L3** MAIA × (Intero − Extero), IEAT ACC sphere (two-sided) | 46 | **Not supported** | Peak z −1.80, p .79. Permutation p .70. Mask-mean r(MAIA) = −.02; JZS regression BF10 0.16 (BF01 6.2; 0.29 at prior scale .354). |

**How the verdicts read.**
- **L1** is a large and consistent effect, but it is a **deactivation**: Intero <
  Extero, relative to the visual-tracking condition. Its suprathreshold voxels are part of one whole-brain negative C1 cluster of
  72,521 voxels (below), which covers 54% of AInsL, 16% of AInsR and 51% of
  dACC, so it is not specific to the salience masks.
  The direction matches prior work: Farb et al. (2023, *eNeuro*) found
  widespread cortical deactivation when attention moved to the breath rather than
  an exteroceptive target. L1 was registered two-sided with no predicted
  direction. The result therefore replicates that pattern and is not a reversal.
- **L2** meets its registered criterion, a voxel FWE peak in the mask, in two of
  three masks. The effect is a single suprathreshold voxel in each mask. Each sits in a
  larger right-hemisphere C5 cluster centred outside the mask (inferior frontal
  and supramarginal in the whole-brain map), and the mask
  averages are flat, with the Bayes factors favouring the null. It is a weak
  replication of H1 (dz ≈ −1 in the breathing runs), and it is fragile across
  the sensitivity analyses (below).
- **L3**: the direct IEAT replication fails, and there is moderate evidence
  against it.

## L1 breathing-robustness rule (registered co-report)
- **Sample:** 29 of the 48 have Biopac `rvt_rrf_z` on both runs, the registered
  set.
- **Comparison:** base and +RVT⊗RRF in the same participants
  (`LOC_L1_breathing_robustness.csv`).
- **Rule:** robust if the +RVT effect keeps its sign and is ≥ 50% of the base
  effect.

| Mask | Base mask mean (dz) | +RVT mask mean (dz) | +RVT / base | r(base, +RVT) | FLAME base / +RVT, p ×3 | Outcome |
|---|---|---|---|---|---|---|
| AInsL | −0.142 (−0.79) | −0.113 (−0.67) | 0.79 | .84 | .003 / .008 | **breathing-robust** |
| AInsR | −0.104 (−0.48) | −0.078 (−0.39) | 0.75 | .90 | .009 / .014 | **breathing-robust** |
| dACC | −0.142 (−0.60) | −0.119 (−0.54) | 0.84 | .93 | .009 / .019 | **breathing-robust** |

RVT removes 16–25% of the effect, and L1 stays significant in every mask with
it. RVT⊗RRF does not model CO2. The PetCO2 subset is reported below.

## Sensitivity analyses (peak z / p_FWE ×k for each mask; bold p < .05)

| | N (L1 / L2 / L3) | L1 AInsL | L1 AInsR | L1 dACC | L2 PIns | L2 COper | L2 SM | L3 sphere |
|---|---|---|---|---|---|---|---|---|
| Primary (200 s) | 48 / 48 / 46 | −4.95 / **6.5e-04** | −4.99 / **4.9e-04** | −5.30 / **1.4e-04** | −2.80 / .20 | −3.96 / **.016** | −4.22 / **.0495** | −1.80 / .79 |
| High-pass 128 s | 48 / 48 / 46 | −4.15 / **.019** | −4.02 / **.027** | −4.63 / **.0032** | −2.71 / .26 | −4.26 / **.0058** | −3.26 / 1 | −1.67 / .87 |
| High-pass 320 s | 48 / 48 / 46 | −5.41 / **7.3e-05** | −5.10 / **3.0e-04** | −5.19 / **2.4e-04** | −3.00 / .13 | −4.31 / **.0044** | −4.19 / .055 | −2.08 / .68 |
| Without protocol v1 (P51–59) | 40 / 40 / 38 | −4.74 / **.0017** | −4.35 / **.0073** | −4.75 / **.0018** | −3.37 / .053 | −3.77 / **.030** | −3.81 / .21 | −2.01 / .70 |
| Protocol covariate | 48 / 48 / 46 | −4.88 / **8.8e-04** | −5.00 / **4.8e-04** | −5.30 / **1.4e-04** | −2.85 / .18 | −3.97 / **.016** | −4.20 / .054 | −1.93 / .75 |
| ΔRVT covariate | 35 / 35 / 34 | −4.38 / **.0074** | −4.12 / **.018** | −4.84 / **.0012** | −3.62 / **.026** | −3.44 / .081 | −4.33 / **.033** | +1.73 / .81 |
| Extero-entrainment covariate | 35 / 35 / 34 | −4.25 / **.012** | −4.21 / **.013** | −4.50 / **.005** | −2.94 / .15 | −3.42 / .084 | −3.90 / .16 | +1.56 / .76 |
| PetCO2 subset, base | 24 / 24 / 24 | −3.69 / .079 | −3.84 / **.044** | −3.85 / .054 | −2.64 / .29 | −2.93 / .31 | −3.84 / .21 | −1.72 / .75 |
| PetCO2 subset, +PetCO2 | 24 / 24 / 24 | −3.23 / .31 | −3.28 / .24 | −3.69 / .094 | −2.08 / .64 | −3.08 / .23 | −3.39 / .84 | −2.13 / .65 |
| RVT subset, base | 29 / 29 / 28 | −4.60 / **.0027** | −4.27 / **.0089** | −4.33 / **.009** | −2.41 / .37 | −3.03 / .21 | −3.30 / .88 | −1.95 / .69 |
| RVT subset, +RVT | 29 / 29 / 28 | −4.33 / **.0083** | −4.17 / **.014** | −4.15 / **.019** | −2.47 / .35 | −3.27 / .12 | −3.56 / .43 | −1.86 / .75 |

**Reading, as registered.**
- **L1** survives:
  - both high-pass alternatives;
  - dropping v1;
  - the protocol, ΔRVT and Extero-entrainment covariates;
  - +RVT.

  In the PetCO2 subset (N 24) it is weak in the base model already (p .04–.08).
  With +PetCO2 no mask reaches .05, and the mask-mean dz moves from −0.67 /
  −0.38 / −0.51 to −0.50 / −0.28 / −0.47. That loss is mostly a subset and
  power loss, but a partial CO2 contribution cannot be excluded at N 24.
- **L2** depends on the specification:
  - COper survives 128 s, 320 s, dropping v1 and the protocol covariate.
  - SM survives only the ΔRVT covariate.
  - PIns reaches .05 only with the ΔRVT covariate.
  - No physiological subset supports L2, base or corrected.

  These divergences are reported as dependence on the specification; none is
  used to select a result.
- **L3** is null in every form.

FLAME converged in every mask (no fallback to permutation was needed). The
largest participant weight share in any test is 0.036.

## Whole brain (descriptive; z > 3.09 per tail, GRF cluster FWE ×2; Harvard–Oxford labels at peaks)

| Contrast | N | Clusters p_FWE < .05 (+ / −) | Largest clusters |
|---|---|---|---|
| C1 Intero − Extero | 48 | 4 / 3 | **−72,521 voxels**, peak z −7.1 at (−40.5, −28.5, 67.5), left postcentral gyrus: a widespread deactivation that includes the suprathreshold voxels of all three L1 masks. Positive clusters are mostly in deep white matter and periventricular regions (33.5, −66.5, 9.5; −22.5, −44.5, 17.5; 3.5, −32.5, 9.5) and left occipital fusiform / cerebellum VI. Lateral occipital cortex and the occipital pole are negative. |
| C2 Intero − Self | 48 | 4 / 2 | −100,059 voxels, peak (−8.5, 19.5, 55.5), superior frontal / SMA (the Self reading and decision network). Positive: periventricular white matter and cerebellum. |
| C3 (conjunction-null min z of C1 and C2) | 48 | not applicable | The conjunction is mostly negative: 64,254 voxels, peak (43.5, −26.5, 61.5), postcentral. Positive conjunction clusters are in white matter / periventricular regions (906 and 825 voxels). |
| C4 Match − Extero | 54 | 2 / 1 | −110,411 voxels, peak z −8.0 at (−34.5, 13.5, 3.5), **left insula** (52% insular cortex). |
| C5 Match − Intero | 48 | 1 / 10 | Negative, right-lateralised: inferior frontal pars opercularis (51.5, 11.5, 7.5; 3,530 voxels), supramarginal / parietal operculum (61.5, −36.5, 35.5; 2,923), occipital fusiform / lingual, superior frontal, lateral occipital and middle temporal. The L2 COper and SM peaks belong to the first two clusters. |
| C6 Self vs Rest | 54 | 11 / 9 | Positive: left inferior frontal / orbital (−52.5, 23.5, −8.5), SMA / paracingulate, right cerebellum Crus I, occipital pole, left temporal pole and posterior superior temporal. Negative: precuneus / cuneus (13.5, −66.5, 25.5) and left central operculum. |
| Intero vs Rest | 48 | 4 / 1 | −116,572 voxels, peak (41.5, −24.5, 51.5), postcentral. Positive: periventricular white matter and cerebellum VI. |
| Extero vs Rest | 54 | 10 / 6 | Positive: left precentral (button presses), SMA, cerebellum V/VI and Crus I, right lateral occipital / MT, and right supramarginal. Negative: calcarine / precuneus (61,169 voxels) and bilateral superior / middle frontal. |
| Match vs Rest | 54 | 3 / 2 | −141,051 voxels, peak (9.5, −74.5, 33.5), precuneus / cuneus. Positive: periventricular white matter and cerebellum VIII. |
| L3 MAIA × C1 | 46 | 0 / 0 | No cluster survives. |

**Interpretation note.** Every condition with active breath control or
attention to the breath has a very large negative cluster against Rest: Intero,
Match, and C1 and C4 against Extero. Each comes with positive white-matter and
periventricular clusters. This is the pattern expected from deeper breathing
(Intero RVT dz .88 over Extero, step 1).
- That makes the risk declared in plan §11 concrete: C1 carries a global
  respiratory component.
- The registered RVT rule passes. The PetCO2 subset is too small to settle the
  CO2 contribution.
- The salience-mask L1 effect should therefore be read as part of a global
  Intero < Extero difference, not as a salience-specific increase or decrease.
  That global deactivation is the pattern reported by Farb et al. (2023, *eNeuro*).
  The open question is how much of it is attentional and how much comes from
  the deeper Intero breathing. The RVT rule says the attentional part survives
  breathing correction; the PetCO2 subset cannot rule out a CO2 contribution.

## Deviations
None in stage A. `Preregistration/localizer_deviations.md` is unchanged.
