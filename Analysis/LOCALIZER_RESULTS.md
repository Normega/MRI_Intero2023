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

Participant-specific ROIs, L4, T1, T3b, T3, T2 and the C2 / C3 fROI rules are stage B
(section "Stage B" below). Connectivity, gPPI, T4, T5, decoding, RSA and PLV
relations are stage C and have not been run.

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

## Deviations (stage A)
None in stage A.

---

# Stage B: participant-specific ROIs, L4 and the transfer tests (OSF p5ztc)

Run 2026-11-27. Choices: `Analysis/localizer_stageB_implementation_notes.md`,
committed (03ef477) before any fROI or stage-B statistic existed. Scripts:
`Scripts/intero2023_loc_breath4mm.py` (breathing-run refit) and
`Scripts/intero2023_loc_froi.py`. Tables: `Analysis/results/LOC_froi_*`,
`LOC_L4_gate*`, `LOC_T1*`, `LOC_T3b.csv`, `LOC_T3.csv`, `LOC_T2.csv`,
`LOC_exploratory_rules.csv`, `LOC_stageB_sensitivity.csv`.

**Contrast names.** The breathing-run files use the swc7p numbering:
- aware − unaware is `C2joint` (H3);
- paced − free is `C1` (H1/H6);
- confidence is `C3` (H2);
- change − No Change is `E3changejoint`.

"C1/C2/C3/C5" without a qualifier below means the **localizer** contrast.

**Pipeline checks.**
- **Refit.** The breathing-run refit code, run at 6 mm, reproduces the swc7p maps
  for sub-060: max |Δ| 0 for C1, C3 and E3changejoint; 4e-9 for C2joint; dof 417
  in both (`Analysis/localizer_stageB_refit_check.json`).
- **Batch.** All 47 refits (4 mm and unsmoothed) completed without error. The
  only absent maps are the registered ones: no C3 for sub-066/070/114, and no
  C2joint for sub-070/114.
- **Extraction.** The 6 mm parcel means equal an independent NiftiMasker
  extraction.

**Sets.** fROI-eligible 44.
- T1: 40, the fROI set ∩ H3. Excluded: sub-056, 080 and 085 (no usable breathing
  run) and sub-070 (never "No Change").
- T3: 41.
- T3b: 40 (sub-065 has no MAIA).
- T2: confidence 39; change − No Change 41.
- No-v1: 34 / 34.

## Verdicts

| Test | N | Verdict | Key numbers |
|---|---|---|---|
| **L4** reliability gate (per parcel; all three criteria) | 44 | **Not met in any parcel** | (1) fails everywhere: the cross-run C1 effect is negative in all five parcels (dz −0.17 to −0.53; 32–41% of participants > 0). (2) Dice median .09–.14 against chance .10 fails, except PIns (.22 against chance .18). (3) passes in AInsL, dACC and PIns, and fails in AInsR (p .21) and the IEAT sphere (p .13). **The LOSO group fallback is therefore primary in AInsR and the sphere.** |
| **T1** aware > unaware in the C1 fROIs, and > complement (one-sided, ×3; FLAME as H3) | 40 | **Not supported** | Gate-primary fROIs: AInsL (individual) dz +0.14, FLAME z 0.92, p ×3 .54; AInsR (fallback) dz +0.33, z 0.90, p .55; dACC (individual) dz −0.04, z 0.88, p .57. fROI − complement dz: −0.13 [−0.44, 0.19], −0.05 [−0.44, 0.24], −0.05 [−0.35, 0.28]. BF01 4.1 / 0.84 / 5.7 (directional BF0+ 2.6 / 0.43 / 7.2). |
| **T3b** MAIA × paced − free > 0 in the ACC-sphere fROI (one-sided) | 40 | **Supported at the registered criterion; weak, and not specific to the fROI** | Gate-primary (fallback) fROI: r .26, slope t 1.63, permutation p .042. JZS BF10 0.57, so the Bayes factor is inconclusive. The rest of the sphere carries the effect more strongly (complement r .38, p .006; whole sphere r .37, p .007). fROI − complement slope −0.011 [−0.049, 0.024]. |
| T3 (secondary) paced − free < 0 in the C5 fROIs (one-sided, ×5) | 41 | Negative in all five parcels; the C5 fROIs do not concentrate it | fROI dz −0.58 (AInsL), −0.76 (AInsR), −0.89 (dACC), −0.82 (PIns), −0.78 (sphere). Sign-flip p ×5 .002 (AInsL), .0005 (others, the floor). Parcel dz −0.78 to −1.13. The fROI is no more negative than its complement anywhere: (fROI − complement) dz on the negated effect −0.21 to +0.04, every CI includes 0. |

**How the verdicts read.**
- **L4.** The two-run C1 fROIs are not reliable participant-specific regions in
  any parcel.
  - **Selection.** Within the selecting data the fROIs look strongly
    breath-preferring:
    - median fROI C1 t 1.0–1.5;
    - fROI C1 effect > 0 in 91–98% of participants;
    - parcel means −0.09 to −0.17%.
  - **Cross-run.** Defined in one run and measured in the other, the same voxels
    are Intero < Extero in every parcel.
  - **Reading.** The group C1 effect in these parcels is a deactivation (stage A,
    L1), and the top 10% of each person's C1 map is mostly selection noise
    around it, not a stable breath-preferring subregion.
  - **Spatial structure.** Within-participant spatial r is low (median .005–.10)
    but beats between-participant r in three parcels. So a weak individual
    spatial pattern exists, just not one that yields a reliable positive fROI.
  - **What the registered rule changes.** The fallback replaces the individual
    fROI only in AInsR and the sphere, because the fallback triggers on (3)
    alone. In AInsL, dACC and PIns the individual fROIs remain primary although
    (1) failed, and every transfer result below carries that caveat.
- **T1.**
  - **Primary result.** Using participant-specific or group-constrained ROIs does
    not reveal an awareness effect that the anatomical masks missed.
    - In AInsL and dACC the fROI Bayes factors favour the null (BF01 4.1 and
      5.7), as do the same parcels' 4 mm means in these 40 (3.1, 5.8) and E5
      (5.1, 5.8). As registered, BF01 > 3 in the fROIs counts as a result: in
      those two masks, spatial dilution does not explain the H3 null.
    - In AInsR the evidence is inconclusive in every form (gate-primary fROI
      BF01 0.84, parcel 1.15, E5 2.2).
  - **Non-primary AInsR fROIs.** Several fROIs that are not primary do reach the
    fROI component in AInsR:
    - the individual C1 fROI, rejected by the gate's (3): dz +0.43, p ×3 .025,
      BF10 4.1;
    - the exploratory C2 and C3 rules: p ×3 .045 and .030.

    None beats its complement (every paired CI includes 0). The AInsR parcel
    itself has dz +0.30 in these 40 (p ×3 .20). The pattern fits a weak,
    parcel-wide right anterior insula awareness effect, the region that was
    already "best but inconclusive" in H3/E1/E5. It is not a localizer-specific
    gain. A post-hoc check (deviations log) gives this individual AInsR fROI
    dz +0.36 to +0.44 across specifications. FLAME p ×3 is .007–.039 except
    unsmoothed (.43), and the paired CI always includes 0.
- **T3b** passes its registered one-sided criterion, but only just.
  - It passes in the group-defined (fallback) fROI, which is the sphere's most
    C1-positive voxels at group level.
  - The Bayes factor is inconclusive.
  - The effect is weaker in the fROI than in the rest of the sphere.
  - It is not robust to the unsmoothed extraction (p .17).
  - The individual sphere fROI (not primary) is null (p .15), and its MAIA slope
    is smaller than the complement's (CI [−0.083, −0.000]).

  The result is best read as the swc7p H6 effect (MAIA × paced − free across the
  IEAT sphere) reappearing at 4 mm in this subset. It is not evidence that a
  localizer-defined ACC region carries it.
- **T3.** Paced − free is strongly negative in every C5-defined fROI, as in H1. The
  localizer Match < Intero voxels are no more pacing-sensitive in the breathing
  runs than the rest of their parcel.

## T1 in detail (gate-primary fROIs, N 40, breathing C2joint at 4 mm)

| Parcel | fROI used | fROI dz / FLAME z / p ×3 | Complement dz / p ×3 | Parcel (anatomical control, same 40, 4 mm) dz / z / p ×3 | fROI − complement dz [95% CI] | dz(fROI) − dz(complement) [95% CI] | fROI BF01 / BF0+ | Parcel BF01 / BF0+ | E5 (N 45, 6 mm) BF01 / BF0+, and registered H3 |
|---|---|---|---|---|---|---|---|---|---|
| AInsL | individual | +0.14 / 0.92 / .54 | +0.19 / .43 | +0.19 / 1.04 / .44 | −0.13 [−0.44, +0.19] | −0.05 [−0.26, +0.14] | 4.09 / 2.56 | 3.05 / 1.75 | 5.06 / 3.44 (H3 z 2.67, p .58) |
| AInsR | LOSO fallback | +0.33 / 0.90 / .55 | +0.29 / .17 | +0.30 / 1.51 / .20 | −0.05 [−0.44, +0.24] | +0.04 [−0.24, +0.30] | 0.84 / 0.43 | 1.15 / 0.60 | 2.18 / 1.18 (H3 z 3.07, p .23) |
| dACC | individual | −0.04 / 0.88 / .57 | −0.01 / 1.00 | −0.01 / 0.33 / 1.00 | −0.05 [−0.35, +0.28] | −0.03 [−0.23, +0.20] | 5.66 / 7.17 | 5.84 / 6.25 | 5.81 / 8.05 (H3 z 2.25, p 1.00) |

- The largest FLAME weight share is 0.049 / 0.050 / 0.076.
- **FLAME versus the unweighted estimate.** FLAME and the unweighted estimate
  differ in places. In the AInsR fallback fROI the unweighted t is 2.08
  (dz +0.33) and FLAME z is 0.90. In dACC the unweighted t is −0.28 and FLAME z
  is +0.88. The registered estimator is FLAME, and neither form reaches
  significance after ×3.
- **Fallback fROI in AInsL** (not primary there): dz +0.32, p ×3 .17, fROI −
  complement +0.23 [−0.08, +0.55].

## L4 gate in detail (C1 rule, N 44, per-run 4 mm maps)

| Parcel | Cross-run C1: mean % (dz), p one-sided, share > 0 | (1) | Dice median (chance) | (2) | r within / between (median), p one-sided | (3) | L4 | Primary fROI |
|---|---|---|---|---|---|---|---|---|
| AInsL | −0.041 (−0.22), .93, .41 | fail | .119 (.100) | fail | .090 / .013, .0004 | pass | not met | individual |
| AInsR | −0.036 (−0.17), .86, .34 | fail | .113 (.101) | fail | .005 / .013, .21 | fail | not met | LOSO fallback |
| dACC | −0.054 (−0.23), .93, .41 | fail | .142 (.101) | fail | .068 / .016, .003 | pass | not met | individual |
| PIns | −0.104 (−0.53), 1.00, .34 | fail | .220 (.175) | pass | .103 / .019, .008 | pass | not met | individual |
| IEAT sphere | −0.090 (−0.31), .98, .32 | fail | .094 (.100) | fail | .048 / .025, .13 | fail | not met | LOSO fallback |

**Distribution of the fROI C1 values** (two-run 4 mm map, the selecting data, so circular):

| Parcel | fROI voxels | fROI mean C1 t, median [IQR] | fROI mean C1 effect % | Share with fROI C1 effect > 0 | Complement C1 effect % | Parcel C1 effect % |
|---|---|---|---|---|---|---|
| AInsL | 59 | 1.39 [0.99, 1.79] | +0.314 | .98 | −0.180 | −0.131 |
| AInsR | 53 | 1.45 [0.92, 1.92] | +0.348 | .98 | −0.138 | −0.089 |
| dACC | 67 | 1.50 [0.83, 2.01] | +0.342 | .95 | −0.181 | −0.128 |
| PIns | 50 | 1.01 [0.39, 1.31] | +0.222 | .91 | −0.231 | −0.151 |
| IEAT sphere | 53 | 1.40 [0.68, 1.87] | +0.363 | .93 | −0.233 | −0.173 |

The LOSO fallback fROIs (top 10% of the group map) have a mean C1 effect of
+0.005% (AInsR) and −0.05 to −0.11% elsewhere (`LOC_froi_definition_summary.csv`).

## T3b in detail (sphere, breathing C1 at 4 mm, N 40)

| Region | r(MAIA) | Slope (% per MAIA point) / t | p one-sided, permutation (parametric) | JZS BF10 |
|---|---|---|---|---|
| fROI, LOSO fallback (**gate primary**) | .26 | 0.075 / 1.63 | **.042** (.056) | 0.57 |
| Complement of that fROI | .38 | 0.086 / 2.50 | .006 (.008) | 2.64 |
| Whole sphere (same 40) | .37 | 0.085 / 2.44 | .007 (.010) | 2.32 |
| fROI, individual C1 (not primary) | .16 | 0.042 / 1.02 | .15 (.16) | 0.28 |

fROI − complement slope: −0.011 [−0.049, +0.024] for the fallback fROI and
−0.048 [−0.083, −0.000] for the individual fROI.

## Exploratory (labelled; cannot support or overturn L4, T1 or T3b)

**C2 (Intero > Self) and C3 (conjunction) fROI rules** (`LOC_L4_gate.csv`,
`LOC_exploratory_rules.csv`).
- **Gate.** Neither passes L4 in any parcel.
  - (1) fails: the cross-run C2 effect is negative in every parcel (dz −0.14
    to −0.96). For C3, both the C1 and the C2 cross-run effects are negative.
  - (2) fails except PIns.
  - (3) passes everywhere (C2: r within .12–.22 against between .05–.09, p ≤ .0011),
    so the individual fROIs are primary.
- **T1.** Not supported for either rule.
  - The AInsR fROI component reaches p ×3 .045 (C2, dz +0.35) and .030
    (C3, dz +0.39), but fROI − complement is +0.14 [−0.18, +0.47] and
    +0.24 [−0.07, +0.57].
  - AInsL: p ×3 .32 / .47. dACC: p .79 / .66, BF01 5.8 / 5.9.
- **T3b.** Not supported for either rule: individual fROIs p .18 / .18. The MAIA
  slope is smaller in the fROI than in its complement (C2 −0.050 [−0.085, −0.002]).
  The non-primary LOSO versions give p .037 / .030.

**T2** (C1 fROIs, gate-primary; FLAME, two-sided, ×3; `LOC_T2.csv`). Null.

| Contrast | N | AInsL dz / p ×3 / BF01 | AInsR dz / p ×3 / BF01 | dACC dz / p ×3 / BF01 |
|---|---|---|---|---|
| Confidence (breathing C3) | 39 | +0.23 / 1.0 / 2.2 | +0.29 / 1.0 / 1.3 | +0.18 / 1.0 / 3.2 |
| Change − No Change (E3) | 41 | +0.14 / 1.0 / 4.0 | +0.13 / 1.0 / 4.3 | −0.12 / 1.0 / 4.4 |

The uncorrected FLAME p is .51–.88 throughout, and every fROI − complement CI
includes 0.

**T3, H6 analogue** (sphere C5 fROI, N 40): r .29, permutation p .043, BF10 0.80.
fROI − complement slope −0.023 [−0.052, +0.016]. This has the same character as T3b.

## Sensitivity analyses (stage B; `LOC_stageB_sensitivity.csv`)

**T1.** Each cell shows fROI dz / p ×3 / fROI − complement dz [95% CI] / BF01.
The primary fROI is the gate's choice.
**T3b.** fallback fROI r / permutation p / fROI − complement slope [CI].

| Specification | N | AInsL | AInsR | dACC | T3b |
|---|---|---|---|---|---|
| Primary (4 mm, 10%) | 40 | +0.14 / .54 / −0.13 [−0.44, +0.19] / 4.09 | +0.33 / .55 / −0.05 [−0.44, +0.24] / 0.84 | −0.04 / .57 / −0.05 [−0.35, +0.28] / 5.66 | +.26 / .042 / −0.011 [−0.049, +0.024] |
| Unsmoothed extraction | 40 | +0.17 / 1.0 / −0.05 [−0.38, +0.26] / 3.50 | +0.44 / 1.0 / +0.06 [−0.30, +0.34] / 0.22 | −0.13 / 1.0 / −0.13 [−0.35, +0.24] / 4.32 | +.14 / .17 / −0.006 [−0.078, +0.100] |
| 6 mm extraction (added) | 40 | +0.14 / .32 / −0.09 [−0.41, +0.23] / 4.10 | +0.25 / .36 / −0.14 [−0.53, +0.17] / 1.91 | −0.00 / .32 / +0.00 [−0.31, +0.32] / 5.86 | +.32 / .014 / −0.009 [−0.038, +0.016] |
| fROI 5% (= 50 voxels) | 40 | +0.13 / .55 / −0.13 [−0.45, +0.18] / 4.25 | +0.33 / .54 / −0.05 [−0.45, +0.24] / 0.80 | −0.06 / .67 / −0.07 [−0.38, +0.25] / 5.51 | +.26 / .038 / −0.012 [−0.049, +0.022] |
| fROI 20% | 40 | +0.18 / .49 / −0.09 [−0.39, +0.23] / 3.30 | +0.28 / .49 / −0.16 [−0.54, +0.15] / 1.36 | −0.02 / .56 / −0.01 [−0.30, +0.34] / 5.83 | +.27 / .030 / +0.000 [−0.040, +0.043] |
| Without protocol v1 | 34 | +0.20 / .54 / −0.07 [−0.41, +0.28] / 2.98 | +0.30 / .63 / −0.14 [−0.62, +0.17] / 1.44 | −0.05 / .59 / −0.06 [−0.36, +0.32] / 5.20 | +.30 / .027 / −0.005 [−0.043, +0.033] |
| Other fROI form, individual vs LOSO (not primary) | 40 | LOSO: +0.32 / .17 / +0.23 [−0.08, +0.55] / 0.95 | individual: +0.43 / .025 / +0.26 [−0.04, +0.57] / 0.24 | LOSO: −0.04 / 1.0 / −0.05 [−0.42, +0.25] / 5.71 | individual: +.16 / .15 / −0.048 [−0.083, −0.000] |

**Reading.**
- **T1** is null in every specification for the gate-primary fROIs. No
  specification gives a fROI − complement CI above 0. In AInsL and dACC, BF01
  stays above 3 in every row except the no-v1 AInsL row (2.98).
- **T3b** keeps p < .05 at 6 mm, 5%, 20% and without v1. It is lost unsmoothed
  (p .17), and its fROI − complement slope never departs from 0.
- **Post hoc (deviations log), AInsR individual C1 fROI**
  (`LOC_T1_posthoc_AInsR_individual.csv`):
  - dz: +0.36 (6 mm), +0.39 (20%), +0.43 (5%), +0.44 (unsmoothed), +0.44 (no v1);
  - FLAME p ×3: .007, .039, .023, .43 and .034;
  - fROI − complement CI: includes 0 in all.

## Deviations (stage B)
Each is logged in `Preregistration/localizer_deviations.md`.
- **Readings before running:**
  - T3 is formed in the five constraint parcels, ×5.
  - The scope of the sensitivity analyses: the high-pass, +RVT/+PetCO2 and
    covariate checks are not repeated for transfer.
  - The analyst definitions of L4 and the fallback.
- **One addition before running:** the 6 mm extraction, at the PI's instruction.
- **One post-hoc descriptive check:** the AInsR individual fROI, added after T1
  was seen and labelled as such.
- No registered test was changed.

## Exploratory: breathing-corrected versions of every localizer map
Logged in `Preregistration/localizer_deviations.md` before it was run. Script:
`Scripts/intero2023_loc_physio_maps.py`. Table:
`Analysis/results/EXPLORATORY_LOC_physio_corrected_maps.csv`. Figure:
`Figures/localizer/loc_physio_corrected.png`.

Each descriptive map (Intero, Extero, Match and Self vs Rest; C1, C2, C4, C5) was
refitted with FLAME 1 from the stage-A +RVT⊗RRF and +PetCO2 first levels. It is
compared with the base model **in the same participants**:
- RVT subset: N 32 for the Intero contrasts, 36 for the others. This is every
  participant with RVT on every used run, so it is larger than the registered
  two-run L1 co-report set of 29.
- PetCO2 subset: N 24 and 26.

**The deactivations survive breathing correction.**

| Map | Whole-brain mean dz, base → +RVT | Retained (RVT) | Whole-brain mean dz, base → +PetCO2 | Retained (CO2) |
|---|---|---|---|---|
| Intero vs Rest | −1.33 → −1.25 | 94% | −1.15 → −1.16 | 99% |
| Match vs Rest | −1.74 → −1.73 | 94% | −1.76 → −1.39 | 94% |
| Extero vs Rest | −0.71 → −0.68 | 100% | −0.83 → −0.84 | 106% |
| C1 Intero − Extero | −0.67 → −0.65 | 88% | −0.61 → −0.59 | 91% |
| C2 Intero − Self | −1.12 → −1.01 | 95% | −1.14 → −1.23 | 107% |
| C4 Match − Extero | −0.84 → −0.73 | 90% | −0.98 → −0.74 | 85% |

- **Salience-region effects.** In the salience regions, C1 keeps 76–91% of its
  mean with +RVT and 84–94% with +PetCO2. Correction moves effects towards zero
  by a within-person change of dz ≈ 0.2–0.4. It never reverses a sign.
- **Extent.** The extent of the C1 negative clusters shrinks:
  - 45,767 → 37,931 voxels with +RVT;
  - 28,432 → 20,390 voxels with +PetCO2.
  
  C4's extent halves with +PetCO2 (59,137 → 28,374 voxels), but its mean is 85%
  retained.
- **The periventricular positive clusters** in Intero, Match and C2 are still
  present after either correction.
- **Self.** The Self vs Rest positives (medial prefrontal, left inferior frontal,
  temporal) are unchanged.
- **Reading.**
  - Breathing regressors explain only a small share of the condition
    differences. That share is largest for the breath conditions (Intero, Match)
    and in the insula and dACC.
  - The widespread Intero deactivation is mostly not accounted for by
    RVT⊗RRF or PetCO2 as modelled. This fits an attentional account (Farb et
    al., 2023).
  - The caveat is that both regressors are imperfect models of respiratory BOLD,
    and they are fitted alongside block regressors they partly share variance
    with.
  - The stage-A loss of L1 significance in the PetCO2 subset reflects the
    smaller N (24), not removal of the effect: the mask means keep 84–94%.

# Stage C1 (exploratory): breath tracking and pacer adherence

**Everything in this section is exploratory** (plan §4; p5ztc "Exploratory analysis"). No result here is a test of a registered hypothesis.
- Definitions were logged in `Preregistration/localizer_deviations.md` and committed before any stage-C value was computed (a247228). Post hoc additions made after the first tables were seen are logged separately and labelled below.
- Script: `Scripts/intero2023_loc_behaviour.py`. It reuses the step-1 functions in `intero2023_localizer_qc.py`; the split into callable helpers reproduces every step-1 CSV exactly.
- Tables: `Analysis/results/EXPLORATORY_LOC_C1_*.csv`. Figure: `Figures/localizer/loc_C1_breath_tracking_adherence.png`.

**Run sets.**
- **Localizer:** the 102 usable runs (motion gate applied). Phase measures also need a usable belt, exact or good trigger timing with pulses = volumes, and a resolved polarity.
- **Change-detection:** the 90 swc7p runs. 58 runs (33 participants; 25 with both runs) have a usable belt and timing. Of the 32 lost, 29 have a failed or missing belt and 3 fail on timing.

## Were participants breathing along with the pacer, in both tasks?

**Yes, nearly all of those who could be measured.** "Above chance" means the participant's own pacer blocks beat their own null blocks (one-sided Mann–Whitney p < .05). Values are medians [IQR] of participant means.

| Measure | N | PLV | Null | PLV > null | Above chance (p < .05) |
|---|---|---|---|---|---|
| Localizer Match (instructed) | 40 participants / 69 runs | .87 [.77–.91] | .32 | 40/40 | **38/40** participants; 58/69 runs |
| Change-detection paced epochs, registered null (free epochs vs a virtual pacer) | 33 / 58 | .82 [.78–.89] | .31 | 33/33 | **33/33**; 55/58 runs |
| Change-detection paced epochs, added length-matched null (other epochs' pacer traces) | 33 / 58 | .82 | .58 | 33/33 | **32/33** (one-sided Wilcoxon over epochs); 55/58 runs |
| Localizer Extero (unintended entrainment) | 40 / 69 | .56 [.36–.68] | .32 | 34/40 | 16/40; 23/69 runs |

- **Both tasks:** 33 participants have usable belts in both. **32/33 are above chance in both.** The exception (sub-054) has a single localizer run with two Match blocks (PLV .76 vs null .43, p .14). That test cannot reach p < .05 with so few blocks.
- **Following the pacer's changes.** Against the length-matched null, paced PLV is still well above it (median .82 vs .58). Participants followed the pacer's rate changes, not just its baseline rate.
- **Agreement with the registered measure.** `adhere_r` has a median of .43, and ≥ .30 in 22/35 participants.
  - Raw paced PLV tracks `adhere_r` closely: ρ .74 [.52, .87], n 33, BF10 71.
  - The registered PLV − null index does not (ρ .10 [−.25, .43]).
- **Extero.** Unintended entrainment is clear at group level (index > 0 in 34/40). For most individuals it is not significant on four blocks.

**Caveat on the registered null (post hoc; logged).** The free-epoch and Rest/Self/Intero nulls compare free breathing with a virtual pacer at the pacer's rate. They therefore score regular breathing near that rate, not only chance.
- The null PLV is a stable participant trait. The localizer null and the change-detection null correlate ρ .51 [.14, .79], n 33.
- It dominates the null-subtracted indices:
  - change-detection index × its free null: ρ −.77 (n 33);
  - Match index × its null: ρ −.57 (n 40).
- That is why the change-detection index is so reliable (ICC .85) and yet unrelated to `adhere_r`. The index with the length-matched null relates to `adhere_r` at ρ .65 [.39, .81].
- The registered indices are kept and reported first. For relating adherence across tasks, raw paced PLV and the length-matched index are more interpretable.

## Breath-tracking accuracy (localizer Intero)

| Measure | N | Median [IQR] | Chance | Above chance |
|---|---|---|---|---|
| Tracking PLV (registered aggregate) | 38 participants / 66 runs | .79 [.71–.87] | .26 (shift mean) | **38/38** participants (p < .05 vs 200 shifts); 61/66 runs |
| Pair score, cross-validated (post hoc, disclosed) | 36 / 61 | .79 [.60–.88] | .39 | 36/36 > chance mean; 61/61 runs |

Run 1 vs run 2 reliability, ICC(3,1):

| Measure | ICC | n | Step-1 QC value |
|---|---|---|---|
| Tracking PLV, registered aggregate | .44 | 28 | — |
| Tracking PLV, step-1 mean-of-blocks form | .44 | 28 | — |
| Tracking PLV, with the motion-failed sub-108 run 1 kept | .43 | 29 | .43 (29) |
| Pair score | .55 | 25 | .56 (26) |

- The step-1 values are reproduced exactly when that run is kept.
- For comparison:
  - Match index ICC .37 (29);
  - Extero index .24 (29);
  - change-detection paced PLV .51 (25);
  - `adhere_r` .49 (26).

## Manipulation check (rate, depth, RVT; step-1 values, not recomputed)

From `EXPLORATORY_LOC_C1_manipulation_check.csv` (= `breathing_by_condition.csv`), N 41 participants:

| | Rest | Self | Intero | Extero | Match |
|---|---|---|---|---|---|
| Rate, brpm (median) | 13.7 | 13.4 | 13.4 | 12.2 | 12.6 |
| Relative depth | 0.75 | 0.79 | 1.40 | 1.04 | 1.80 |
| RVT z | −0.42 | −0.35 | 0.53 | −0.20 | 0.74 |

- **Intero vs Extero:**
  - rate does not differ (dz .31, p .13);
  - Intero breathing is deeper (depth dz .59, RVT dz .88, both p < .001).
- **Match is the deepest condition:** RVT vs Extero dz 1.16; vs Intero dz .38, p .018.

## Exploratory relations

- **Statistics:** Spearman ρ (permutation p), bootstrap 95% CI, Pearson r, and JZS BF10 (s = 1).
- **Correction:** BH-FDR across the 18 listed relations. **None survives FDR** (smallest q .14).

| Relation | n | ρ [95% CI] | r | BF10 | q (FDR 18) |
|---|---|---|---|---|---|
| Tracking × MAIA | 33 | −.00 [−.35, .36] | −.04 | 0.14 | .99 |
| Tracking × detection rate | 34 | .21 [−.11, .48] | .18 | 0.22 | .73 |
| Tracking × mean confidence | 34 | .02 [−.34, .39] | −.06 | 0.14 | .99 |
| **Tracking × breathing-run adherence (`adhere_r`)** | 33 | **.46 [.11, .74]**, p .008 | .42 | 2.73 | .14 |
| Tracking × C1 in the fROIs (cross-run; AInsL, AInsR, dACC, PIns, sphere) | 31 | −.23 to .09; every CI spans 0 | −.29 to .00 | 0.14–0.48 | ≥ .73 |
| Tracking × C1 in the masks (joint 6 mm; same five) | 34 | −.25 to −.05; every CI spans 0 | −.18 to −.07 | 0.14–0.22 | ≥ .73 |
| Localizer Match index × `adhere_r` | 34 | .15 [−.23, .49] | .19 | 0.23 | .78 |
| `adhere_r` × H1 C1 (paced − free), PIns / COper / SM mask means | 35 | .27 [−.11, .60] / .17 [−.18, .49] / .02 [−.32, .37] | .24 / .08 / −.03 | 0.35 / 0.14 / 0.13 | .73 / .78 / .99 |

- **Strongest relation:** tracking × `adhere_r`. Better breath trackers in the localizer followed the change-detection pacer more closely. It does not survive FDR, and the BF10 of 2.7 is anecdotal.
- **Tracking is unrelated to MAIA, detection, confidence and C1.**
  - Most BF10 values are 0.14–0.22, i.e. BF01 of about 4.5–7: moderate evidence for no linear relation.
  - The fROI and mask rows are less conclusive: AInsR fROI BF10 0.48, PIns fROI 0.28.
  - The secondary pair score gives the same picture: |ρ| ≤ .21 with MAIA, detection, confidence and `adhere_r` (it was not related to C1).
  - The AInsR gate-primary LOSO fROI: ρ −.09.
- **Adherence across tasks.**
  - The registered Match index is only weakly related to `adhere_r`.
  - Same-scale rows (sensitivity, outside the family):
    - Match index × change-detection index: ρ .50 [.18, .72], BF10 13. Part of this is the shared null-PLV trait.
    - Raw Match PLV × raw change-detection paced PLV (post hoc): ρ .64 [.37, .82], BF10 20.
    - Raw Match PLV × `adhere_r` (post hoc): ρ .36 [−.01, .66].
  - People who followed the pacer closely in one task tended to do so in the other.
- **Other descriptive and post hoc rows:**
  - tracking × Extero entrainment: ρ .32 [.02, .56], n 38, BF10 1.1;
  - tracking × change-detection paced PLV (post hoc): ρ .49 [.14, .75], BF10 4.2.

## Breathing-run adherence as a moderator of H1 (paced − free)

This is the swc7p H6 group model (nilearn permutation, 10,000, seed 20261126) with demeaned `adhere_r`. It is run two-sided in each H1 mask, Bonferroni ×3, on N 35 (the H1 set ∩ adherence available).

| Moderator | N | PIns peak t, p_FWE ×3 | COper | SM | Mask-mean ρ (PIns / COper / SM) |
|---|---|---|---|---|---|
| `adhere_r` (listed) | 35 | +2.98, .99 | +2.70, 1.0 | −3.60, 1.0 | .27 / .17 / .02 |
| Change-detection index (sensitivity) | 33 | −2.24, 1.0 | −1.95, 1.0 | −2.98, 1.0 | −.05 / −.07 / −.13 |
| Change-detection paced PLV (post hoc) | 33 | +3.93, .53 | −3.22, 1.0 | −4.46, 1.0 | .26 / .18 / .05 |

- **No moderation.** No voxel survives, and the mask-mean relations are small with CIs spanning 0 (BF10 0.13–0.35 for `adhere_r`).
- The PIns mask mean leans positive: closer pacer following goes with less paced-breathing deactivation (ρ .27 [−.11, .60]).
- H1's deactivation (mask C1 dz −0.82 to −1.06 in this subset) therefore does not depend detectably on how closely people followed the pacer.

## Reading

- Participants were breathing along with the pacer in both tasks:
  - 38/40 in the localizer Match blocks;
  - 33/33 in the change-detection paced epochs;
  - 32/33 in both.
- They also tracked their own breath far above chance (38/38), with modest run-to-run reliability (ICC about .44).
- Individual differences in tracking do not relate to MAIA, change detection, confidence or the localizer C1 response in this sample.
- The one moderate signal links tracking with pacer adherence. It is exploratory, does not survive FDR, and may reflect general task engagement or breathing regularity.
- **Power.** The relation Ns are 31–38 and the measures' reliabilities are about .4–.5. Only correlations of about |ρ| ≥ .45 were detectable with reasonable power, so small true relations cannot be excluded.

## Deviations (stage C1)

- **Analyst definitions and readings, logged before running.** The main ones:
  - the motion gate is applied, so tracking N is 38 participants / 66 runs (the registered "67 runs" includes the motion-failed run);
  - the change-detection belt is located from the raw `.acq` with the step-1 rule;
  - the FDR family;
  - the H1 moderation model.
- **One added sensitivity null, logged before running:** the length-matched null.
- **Post hoc additions after the first tables were seen** (raw paced PLV, null-term diagnostics, paced PLV as moderator). They are logged and labelled and lie outside the FDR family.
- No registered result was changed.
