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

## Exploratory: localizer model without the CSF and white-matter regressors
Logged before running in `Preregistration/localizer_deviations.md`.
- Script: `Scripts/intero2023_loc_notissue.py`
- Table: `Analysis/results/EXPLORATORY_LOC_notissue_comparison.csv`
- Figure: `Figures/localizer/loc_notissue_comparison.png`

The 6 mm joint base model was refitted for all 54 participants with `csf` and `white_matter` removed. Nothing else changed. FLAME 1 was run in the stage-A sets (48 / 54).

**Why test this.** The fMRIPrep CSF regressor is itself task-locked. Regressed on the task design, it falls in every task block relative to Rest: dz −0.30 (Intero), −0.60 (Extero), −0.53 (Match), −0.36 (Self); N 48. The white-matter signal is mostly flat: Match dz −0.41, Self +0.46, Intero and Extero ≈ 0.

**Result 1: the CSF covariate does not create the periventricular positives.** They remain without the tissue regressors, at the edges of the lateral ventricles. Mean effect in the registered model's positive clusters, with tissue regressors → without:
- Intero vs Rest: +0.156 → +0.150 PSC;
- C1: +0.148 → +0.163;
- Match vs Rest: +0.137 → +0.076.

Their extent shrinks (Intero 8,151 → 4,831 voxels; Match 9,374 → 2,524), mainly because between-subject variance rises.
- The Harvard–Oxford lateral-ventricle mask as a whole falls in both models. The positive signal is confined to the ventricle borders.
- One possibility is the known anti-coupling between global grey-matter BOLD and CSF signal. This was not tested here.

**Result 2: for Intero and C1, the tissue regressors mainly reduce noise.**

| Map | Cortex mean PSC (with → without) | Cortex dz (with → without) |
|---|---|---|
| Intero vs Rest | −0.233 → −0.223 | −1.20 → −0.71 |
| C1 Intero − Extero | −0.110 → −0.088 | −0.71 → −0.28 |

- In the L1 masks, C1 means change by ≤ 0.03 PSC, but dz roughly halves: AInsL −0.65 → −0.28, AInsR −0.41 → −0.14, dACC −0.55 → −0.22.
- The negative C1 extent falls from 73,074 to 33,622 voxels.
- The registered model is therefore the more sensitive one for L1. Its estimate is not inflated.

**Result 3: for the Match contrasts, the CSF regressor absorbed part of the effect.** Without the tissue regressors, the Match deactivation deepens:
- cortex −0.277 → −0.373 PSC;
- Match − Extero (C4) cortex −0.171 → −0.250.

Match − Intero (C5, the L2 / H1-replication contrast) is much stronger:

| Measure | With tissue regressors | Without |
|---|---|---|
| Whole-brain dz | −0.30 | −0.61 |
| PIns dz | −0.04 | −0.39 |
| SM dz | −0.13 | −0.50 |
| Negative extent | 12,755 voxels | 68,082 voxels |

- So the weak, peak-level L2 result is partly due to the registered nuisance model. Paced breathing produces a task-locked global and CSF signal change, and the CSF regressor removes part of it.
- Whether the removed part is neural or vascular cannot be decided from these data.
- The same CSF regressor is in the registered breathing-run models (swc7p), where H1 was supported regardless.

**Result 4: Self vs Rest becomes more positive overall without the tissue regressors.** Whole-brain dz moves from −0.09 to +0.25.

The registered model stays primary. None of the registered verdicts changes.

# Stage C2 (exploratory): connectivity

**Everything in this section is exploratory** (plan §6, §8 T5, §9). No result here tests a registered hypothesis, and none changes an earlier result.
- Choices and readings were committed before any connectivity value was computed (9228853): `Analysis/localizer_stageC2_implementation_notes.md` and `Preregistration/localizer_deviations.md`.
- Scripts: `Scripts/intero2023_loc_connectivity.py` (network FC, background FC, T5) and `Scripts/intero2023_loc_gppi.py` (gPPI).
- Tables: `Analysis/results/EXPLORATORY_LOC_C2_*.csv`.
- Figures: `Figures/localizer/loc_C2_network_fc.png` and `loc_C2_gppi_maps.png`.

**Method in one paragraph.** Network FC uses the registered H7 pipeline (swc7p), imported unchanged:
- unsmoothed BOLD, band-passed at 0.01–0.1 Hz;
- the median signal of each CAREN network;
- block windows shifted by 4 s, with FD > 0.5 mm volumes dropped;
- each window demeaned, then windows concatenated within condition (median 59–60 volumes per condition);
- Fisher z per pair, then a paired t on each condition difference.

There are two versions:
- **network FC:** nuisance regression only, so evoked responses stay in the signal;
- **background FC:** the localizer task regressors are also removed. This is the form identical to H7.

The plan lists both. How they were read is logged in the deviations file.

**Checks.**
- **H7 reproduction.** The T5 code recomputes the H7 network values for all 41 T5 participants exactly (max |Δz| 2e-11; `Analysis/localizer_stageC2_verification.json`).
- **gPPI design.** Each condition boxcar, convolved on the microtime grid, matches nilearn's regressor (r ≥ .9997).
- **Deconvolution.** On a simulated neural signal it recovers r .87–.997, depending on noise (`Analysis/localizer_gppi_design_check.json`). On the real seeds, the reconvolved fit is median r .77–.88 by seed (`EXPLORATORY_LOC_C2_gppi_firstlevel_qc.csv`).
- **Volume rule.** The pre-set rule (≥ 2 windows and ≥ 20 volumes per condition) removed two one-run participants from the Intero contrasts, sub-058 and sub-108, so N is 46. It removed sub-058 from C4, so N is 53.

## Network FC: primary pairs (plan §6)

Values are dz of the paired Fisher-z difference, with two-sided p. Holm ×2 applies within C1.

| Contrast, pair | Network FC (nuisance only) | + global signal partialled | Background FC (task removed) | + global signal partialled |
|---|---|---|---|---|
| C1 Intero − Extero, SAL–SMN (N 46) | +0.45, p .0035 (Holm .0035) | +0.09, p .54 | +0.29, p .058 (Holm .058) | +0.12, p .44 |
| C1 Intero − Extero, SAL–DMN (N 46) | +0.69, p 3e-5 (Holm 5e-5) | +0.58, p .0003 | +0.48, p .0020 (Holm .0040) | +0.51, p .0013 |
| C2 Intero − Self, SAL–DMN (N 46) | +0.75, p 7e-6 | +0.30, p .046 | +0.37, p .016 | +0.25, p .10 |

- **SAL–DMN.** Coupling is higher during Intero than Extero. This holds:
  - in both analyses;
  - with the global signal partialled out;
  - with the five-network mean partialled out (network FC +0.44, p .0045; background +0.41, p .008).

  By the pre-set global-signal rule it can be called a network-coupling difference.
  - After the global-signal partial, SAL and DMN are anticorrelated in both conditions. Intero is less anticorrelated: background −0.20 vs −0.38 in mean z.
  - The difference therefore means "less anticorrelated during Intero", not "more positively coupled".
- **SAL–SMN.** The C1 difference is present in network FC. It is weaker in background FC (p .058). It disappears in both once the global signal is partialled out (p .54 / .44), so it is **not separable from the global signal**.
- **C2 SAL–DMN.** Strong in network FC. Smaller in background FC (p .016). After the global-signal partial: p .046 in network FC and .10 in background FC.
- **Specifications.** The same pattern holds for median or mean signal and for a 4 or 6 s shift. C1 SAL–DMN p ≤ .002 in every specification, both analyses, with and without the global-signal partial. C1 SAL–SMN with the global signal partialled out has p ≥ .44 in every specification, both analyses.

## Network FC: all 10 pairs (BH-FDR over 10; `EXPLORATORY_LOC_C2_netfc.csv`)

| Contrast (N) | Network FC: pairs FDR < .05 | Background FC: pairs FDR < .05 | Background FC, global signal partialled |
|---|---|---|---|
| C1 Intero − Extero (46) | 6, all Intero > Extero: SAL–SMN, SAL–ECN, SAL–DMN, SMN–ECN, SMN–DMN, ECN–DMN | 4, all Intero > Extero: SAL–DMN, SMN–ECN, SMN–DMN, ECN–DMN (largest ECN–DMN +0.64, p 8e-5) | 5. Intero > Extero: SAL–DMN, SMN–DMN, ECN–DMN (+0.80, p 2e-6). Intero < Extero: SAL–VIS −0.40, SMN–VIS −0.42 |
| C2 Intero − Self (46) | 8, all Intero > Self | 6, all Intero > Self | 0 |
| C4 Match − Extero (53) | 3, all Match > Extero: SAL–SMN, SAL–ECN, SMN–ECN | 0 | 0 |
| C5 Match − Intero (46) | 2, both Match < Intero: SAL–DMN −0.44, ECN–DMN −0.55 | 1: ECN–DMN −0.48 | 2: SAL–DMN −0.43, ECN–DMN −0.57 |

- The robust pattern involves the DMN. Its coupling with SAL, SMN and ECN is higher during free-breath tracking (Intero) than during either pacer condition: C1 against Extero, and C5 (Match lower than Intero).
- It survives task removal and the global-signal partial.
- The C2 (vs Self) and C4 (Match vs Extero) differences do not survive the global-signal partial in background FC.

## Amplitude, co-activation, and what they mean for the coupling differences (H7 checks)

**Amplitude** is the temporal SD per network, C1 Intero − Extero:
- Network FC: SD is higher during Intero in SAL (dz +0.54, p .001), SMN (+0.40, p .010) and ECN (+0.44, p .005). DMN and VIS do not differ (+0.01, +0.20).
- Background FC: no network differs (dz −0.25 to +0.21, all p ≥ .10).

So once the evoked responses are removed, the coupling differences are not accompanied by amplitude differences.

**Co-activation** is the network-mean stage-A C1 effect. Intero is lower than Extero in SAL (dz −0.57), SMN (−0.77), ECN (−0.68) and DMN (−0.91), all p < .001; VIS −0.08. This is the stage-A deactivation. The coupling increases occur alongside lower mean activity.

## Physiology and the CSF regressor

The physiology regressor is added, and the model is compared with the base model in the same participants. Values are dz (p).

| Pair, analysis | RVT set (N 31): base → +RVT | PetCO2 set (N 23): base → +PetCO2 |
|---|---|---|
| C1 SAL–DMN, network FC | +0.65 (.0011) → +0.54 (.0053) | +0.56 (.013) → +0.41 (.063) |
| C1 SAL–DMN, background FC | +0.48 (.012) → +0.44 (.020) | +0.32 (.14) → +0.36 (.094) |
| C1 SAL–SMN, network FC | +0.34 (.067) → +0.30 (.10) | +0.49 (.027) → +0.40 (.068) |
| C1 SAL–SMN, background FC | +0.24 (.20) → +0.19 (.31) | +0.22 (.29) → +0.26 (.23) |
| C2 SAL–DMN, network FC | +0.83 (.0001) → +0.79 (.0001) | +0.73 (.0019) → +0.64 (.0056) |
| C2 SAL–DMN, background FC | +0.59 (.0027) → +0.66 (.0010) | +0.43 (.0496) → +0.50 (.025) |

- **RVT⊗RRF** changes the primary pairs little.
- **PetCO2.** In the small PetCO2 set, the C1 SAL–DMN network-FC effect falls from dz 0.56 to 0.41. Its background-FC effect was already weak in that subset. A partial CO2 contribution cannot be excluded at N 23.
- **CSF regressor (`nocsf` sensitivity).** The fMRIPrep CSF regressor is task-locked: it falls in task blocks relative to Rest. Dropping it shifts the numbers slightly:
  - the network-FC amplitude differences grow (SAL dz +0.64, SMN +0.44, ECN +0.42);
  - background C1 SAL–SMN becomes p .016 (dz +0.37);
  - with the global signal partialled out, that pair stays null (+0.21, p .17);
  - C1 SAL–DMN is unchanged (network FC +0.69; background FC +0.53, p .0009).

  Removing the regressor therefore takes out a little condition-locked variance. It does not change any reading above.

## gPPI (seeds: gate-primary fROIs and a-priori masks; whole brain; FLAME 1, cluster FWE at z > 3.09, two-sided)

- 25 clusters survive, in 11 of the 20 seed × contrast maps.
- 11 of them, in 6 maps, would also survive a further ×20 across maps (p < .0025, marked †).
- Positive = stronger coupling with the seed in the first condition.
- `EXPLORATORY_LOC_C2_gppi_clusters.csv`; Harvard–Oxford labels at the peak.

| Seed (N) | Contrast | Clusters: sign, voxels, peak MNI, label, p_FWE |
|---|---|---|
| fROI AInsL (44) | C1 | +158, (39.5, 31.5, −6.5), R frontal orbital / IFG pars triangularis, 9e-6†; +61, (33.5, 11.5, 63.5), R MFG, .021 |
| fROI AInsL (44) | C4 | +73, (49.5, 33.5, 25.5), R MFG / frontal pole, .012 |
| fROI AInsR (LOSO; 44) | C1 | +215, (51.5, 29.5, −4.5), R IFG pars triangularis, 8e-6†; +122, (39.5, 35.5, −12.5), R frontal pole / orbital, .0014†; +78, (53.5, 17.5, 15.5), R IFG pars opercularis, .024; +78, (−2.5, −36.5, 55.5), medial pre/postcentral, .024; +77, (41.5, −16.5, −28.5), R temporal fusiform, .026 |
| fROI AInsR (LOSO; 44) | C4 | +73, (−56.5, 5.5, 15.5), L precentral, .046 |
| fROI dACC (44) | C1 | +86, (23.5, 57.5, 23.5), R frontal pole, .008; +63, (45.5, −60.5, 25.5), R lateral occipital / angular, .047 |
| fROI dACC (44) | C4 | −133, (−4.5, −76.5, −6.5), L lingual, 8e-5†; +91, (25.5, 57.5, 25.5), R frontal pole, .0019† |
| fROI PIns (44) | C4 | −196, (1.5, −80.5, −0.5), lingual / intracalcarine, .0002†; −160, (−4.5, 47.5, −12.5), medial frontal, .0011†; +94, (−62.5, −48.5, 33.5), L supramarginal / angular, .031 |
| fROI sphere (LOSO; 44) | C1 | +113, (−56.5, 7.5, 15.5), L precentral / IFG pars opercularis, .0064; +91, (45.5, −44.5, 57.5), R supramarginal / superior parietal, .023 |
| mask AInsR (48) | C1 | +150, (51.5, 45.5, −4.5), R frontal pole, .0014† |
| mask dACC (54) | C4 | −97, (1.5, −86.5, −2.5), lingual / intracalcarine, .021; −91, (37.5, −88.5, 15.5), R lateral occipital / occipital pole, .030 |
| mask PIns (54) | C4 | −331, (−4.5, 47.5, −14.5), medial frontal, 1e-7†; −234, (−18.5, −96.5, −10.5), L occipital pole, 5e-6†; −171, (−6.5, −76.5, −6.5), L lingual, .0001†; −114, (23.5, −72.5, −4.5), R occipital fusiform, .0029 |

These maps have no surviving cluster:
- fROI PIns C1 (smallest p .17);
- fROI sphere C4 (.077);
- masks AInsL (C1 .47, C4 .61), dACC C1 (.76), PIns C1 (.29) and sphere (C1 .069, C4 .050);
- mask AInsR C4 (.17).

- **C1 (Intero − Extero).** Both anterior-insula seeds couple more with right ventrolateral prefrontal cortex during Intero: right IFG, orbital and frontal pole. So does the a-priori AInsR mask.
  - For the AInsR seeds this cluster lies next to the seed, so local spread at 6 mm may contribute.
  - For the AInsL fROI it is contralateral.
  - The dACC and sphere fROIs give smaller frontal-pole, precentral and parietal clusters.
  - Every C1 cluster is positive.
- **C4 (Match − Extero).** The posterior-insula and dACC seeds (fROI and mask) couple **less** with visual cortex (lingual, intracalcarine, occipital pole) and with ventromedial prefrontal cortex during Match. Both conditions show the same pacer circle and take presses.
- **fROI vs mask.** The fROI seeds give surviving C1 clusters more often than the a-priori masks: 4 of 5 fROI seeds against 1 of 5 masks. The fROI seeds are defined on C1, so this is expected and is not independent evidence.
- **Within-network summaries** (mean gPPI effect per CAREN network, seed voxels excluded; `EXPLORATORY_LOC_C2_gppi_networks.csv`). None of the 100 survives FDR (over 50 per contrast). The smallest is sphere fROI C1 → ECN, dz +0.46, p .0038, q .19. 8 of 100 have uncorrected p < .05.
- **+RVT check (rule fixed before running; `EXPLORATORY_LOC_C2_gppi_rvt.csv`).** Every seed with a surviving cluster was refitted with RVT⊗RRF in the participants with RVT (29 fROI, 32 or 36 mask). In each of the 25 clusters, the mean effect keeps its sign. +RVT retains 63–164% of the base mean in the same participants. The clusters were selected on the base maps, so the base means are inflated. The check concerns the change with RVT, not the size.
- **Deconvolution caveat.** GCV regularises heavily (median λ 0.20–0.57 × s_max²). In 25 of the 950 seed-runs, λ sits at the top of the grid. The deconvolved signal is then shrunk towards zero, and those runs carry little PPI information.

## T5: fROI-node coupling in the breathing runs, compared with H7 (`EXPLORATORY_LOC_C2_T5*.csv`)

**Nodes.** The five gate-primary fROIs, coupled with the CAREN networks (with the fROI voxels removed) and with each other. That gives 35 pairs, with BH-FDR over 35. The pipeline is H7 Model A.

| | T5a paced − free (N 41) | T5b aware − unaware (N 21) |
|---|---|---|
| fROI pairs FDR < .05 | 1: **AInsL–AInsR lower under pacing**, dz −0.58, p .0006, q .022 (mean z 0.10 vs 0.20) | 1: **AInsL fROI–DMN higher when aware**, dz +0.92, p .0004, q .015 (0.28 vs 0.00) |
| Global signal partialled | −0.51, p .0022 (q .077) | +0.67, p .0063 (q .22) |
| Five-network mean partialled | −0.48, p .0040 (q .14) | +0.36, p .12 |
| Other specifications (mean signal, 6 s shift) | q ≤ .032 in all four | q ≤ .015 in all four. AInsL–VIS also reaches q .031 at 6 s |
| RVT set: base → +RVT | N 24: −0.67 (.0035) → −0.65 (.0040) | N 13: +0.92 (.0061) → +0.75 (.019) |
| PetCO2 set: base → +PetCO2 | N 21: −0.51 (.031) → −0.37 (.11) | N 12: +0.85 (.013) → +0.60 (.063) |
| fROI amplitude (SD) | lower under pacing in all five fROIs, dz −0.33 to −0.50, p .003–.041 | dACC −0.46 (p .047); others p ≥ .24 |
| **H7 network pairs, same participants** (from the registered H7 files) | SMN–VIS +0.48, p .0041, q .041 (in the full 47, q .055); SAL–VIS +0.38, p .021; the registered pairs SAL–SMN −0.20 (p .21) and SAL–ECN +0.10 (p .52) are null | none FDR < .05; SAL–DMN +0.53, p .024 (q .21); SAL–SMN −0.02, SAL–ECN +0.12 |

**How T5 compares with H7.**
- **Paced − free.** The fROI level shows one difference that H7's network level did not have: the two anterior-insula fROIs decouple under paced breathing.
  - It comes with lower fROI amplitude under pacing, as H7 found for whole networks.
  - It holds with the global signal partialled out (p .0022), so by the pre-set rule it is a coupling difference. It does not pass FDR after either partial.
  - It weakens in the PetCO2 subset.
  - Pacing lowers both coupling and signal amplitude here. The H7 caution about amplitude-driven correlation changes applies in reverse: lower variance with stable noise would also lower measured correlation.
- **Aware − unaware.** The AInsL fROI–DMN result has the same direction as H7b's SAL–DMN (p .024, not FDR), but is larger.
  - It rests on 21 participants.
  - It weakens after the global-signal partial (p .0063, q .22) and the five-network partial (p .12).
  - It shrinks with +RVT and +PetCO2 in subsets of 13 and 12.
  - It should be treated as a lead to check, not a finding.
  - Neither T5 result is the registered H7 pair (SAL–SMN, SAL–ECN). At the fROI level, those pairs remain null.

## Reading

- **Network FC.** During free-breath tracking (Intero), the default-mode network is less anticorrelated with the salience network, and more coupled with the somatomotor and executive networks, than during either pacer task (Extero, Match).
  - This survives removal of the evoked responses, the global-signal and five-network partials, and RVT.
  - It occurs alongside the stage-A deactivation, with no amplitude difference in background FC.
  - The plan's other primary pair, C1 SAL–SMN, is a global-signal effect.
  - Intero vs Self (C2) and Match vs Extero (C4) differences do not survive the global-signal partial in background FC.
- **gPPI.**
  - The anterior-insula seeds couple more with right ventrolateral prefrontal cortex during Intero than Extero.
  - Posterior-insula and dACC seeds couple less with visual and ventromedial prefrontal cortex during Match than Extero.
  - The clusters are modest (61–331 voxels), are not corrected across the 20 maps beyond the ×20 count above, and are stable to +RVT.
- **T5.** The localizer fROIs add two exploratory leads to H7's null network results: anterior-insula decoupling under pacing, and AInsL–DMN coupling with awareness (N 21). Neither is a registered pair, and the awareness result is fragile.
- **Caveats** (plan §11).
  - C1 also differs in agency, press predictability and breathing (Intero breathing is deeper, RVT dz .88).
  - Corrections for RVT⊗RRF and PetCO2 are partial models of respiratory BOLD, and the PetCO2 subsets are small.
  - "Network FC" and "background FC" follow the reading logged before running.

## Deviations (stage C2)

All were logged in `Preregistration/localizer_deviations.md` before running (one entry, 9228853). They are:
- the network-FC / background-FC reading;
- added contrasts C4 and C5;
- the volume rule;
- the specifications without FIR;
- the added `nocsf` sensitivity;
- the gPPI seeds, sets, deconvolution, multiplicity and +RVT rule;
- T5 nodes and voxel exclusions;
- the global-signal reading rule.

The gPPI figure shows only maps with a cluster that also survives ×20. That is a display choice; all 25 clusters are in the table. No post hoc analysis was added, and no earlier result was changed.

# Stage C3 (exploratory): pattern expression, decoding, RSA, gradients, respiratory phase

**Everything in this section is exploratory** (plan §8 T4, §9, §10). No result here tests a registered hypothesis, and none changes an earlier result.
- Choices were committed before any stage-C3 value was computed (200a864): `Analysis/localizer_stageC3_implementation_notes.md` and `Preregistration/localizer_deviations.md` (entry "Stage C3").
- Post hoc notes made after the tables were seen are logged separately and labelled below.
- Scripts: `Scripts/intero2023_loc_t4.py` and `intero2023_loc_t4_mixed.R` (T4), `intero2023_loc_physio_relations.py`, `intero2023_loc_decoding.py`, `intero2023_loc_rsa.py`, `intero2023_loc_gradient.py`, `intero2023_loc_respphase.py`, `intero2023_loc_c3_figures.py`.
- Tables: `Analysis/results/EXPLORATORY_LOC_C3_*.csv`. Figures: `Figures/localizer/loc_C3_t4_decoding.png` and `loc_C3_rsa_gradient.png`.

**Checks.**
- **LSS designs.** For sub-060, all 28 LSS designs are full rank, and the target column plus the remaining same-condition column reproduce the Model A column exactly (`Analysis/localizer_stageC3_t4_design_check.json`). All 44 participants' LSS fits completed without error.
- **RSA refit.** The per-run localizer betas refitted in the ROI voxels equal the saved stage-A per-run 4 mm maps (max |Δ| 7.6e-7 over 44 participants).

## T4: does the breathing-run epoch pattern look like the localizer C1 pattern, and does that track awareness or confidence?

**Method in brief.**
- **Betas.** LSS betas for every paced (595) and free (594) epoch of 44 participants: the localizer Intero set ∩ the swc7p set. The registered Model A design is used, at 4 mm, in CAREN SAL ∪ SMN (13,189 voxels).
- **Index.** Each epoch's index is its projection on the participant's own two-run localizer C1 map. That map is Intero − Extero, mostly negative.
  - **Primary, dot:** the dot product with the unit-length map, divided by √V, in PSC.
  - **Secondary, corr:** the spatial correlation.
- **How the sign reads.** A positive dot index means the epoch is lower where Intero is lower than Extero, so it looks Intero-like, *including the overall deactivation*. A positive corr index means the same spatial shape, whatever the epoch's mean level.
- **Baseline caveat.** The breathing runs have no clean baseline, so only comparisons between epochs within participant are interpretable.
- **Model:** the registered `index ~ confidence + outcome + (1 | participant)`, on 518 paced epochs with outcome aware, unaware or No Change. Wrong-direction epochs are excluded.

**Result: no awareness or confidence effect** (`EXPLORATORY_LOC_C3_T4_mixed.csv`). Estimates are in index units, with 95% CI and Satterthwaite p.

| Index, template | Confidence (per point) | Aware − unaware | Aware − No Change | Unaware − No Change | Outcome F |
|---|---|---|---|---|---|
| **dot, own (primary)** | −0.001 [−0.014, 0.012], p .86 | −0.017 [−0.047, 0.013], p .28 | −0.012 [−0.045, 0.021], p .47 | +0.004 [−0.033, 0.042], p .82 | F(2, 477) 0.68, p .51 |
| corr, own | +0.005 [−0.006, 0.016], p .40 | −0.008 [−0.034, 0.017], p .52 | −0.003 [−0.031, 0.025], p .84 | +0.005 [−0.026, 0.037], p .74 | F(2, 478) 0.20, p .82 |
| dot, LOSO group map | +0.000 [−0.017, 0.018], p .96 | −0.021 [−0.061, 0.019], p .31 | +0.001 [−0.044, 0.046], p .96 | +0.022 [−0.029, 0.073], p .39 | F(2, 484) 0.59, p .56 |
| corr, LOSO group map | +0.008 [0.001, 0.015], p .022 | −0.011 [−0.028, 0.005], p .18 | +0.003 [−0.015, 0.022], p .71 | +0.015 [−0.006, 0.035], p .16 | F(2, 480) 1.24, p .29 |

- **Sensitivity (dot, own).** The outcome contrasts stay null:
  - with run intercepts: aware − unaware −0.001 [−0.026, 0.025], p .96;
  - with confidence split into within- and between-participant parts: within −0.002 [−0.015, 0.012], between +0.009 [−0.051, 0.069];
  - with wrong direction as a fourth outcome: every contrast p ≥ .26.
- **Per network.** SAL alone and SMN alone are also null: outcome F p .71 and .46 (dot), aware − unaware p .57 and .22.
- **The one p < .05** is the confidence slope for corr with the LOSO template (+0.008, p .022). It is one of 16 confidence and contrast estimates in that table's primary and secondary rows, and it is not seen with the participant's own map. It is not read as a finding.
- **The index does vary, between epoch types.** In the added paced vs free check, paced epochs express the C1 pattern more than free epochs:
  - dot +0.070 [0.056, 0.083], p 1e-22 (free epochs' level +0.056 [0.014, 0.099]);
  - corr +0.055 [0.044, 0.067], p 2e-20.

  So the C1 pattern is closer to the paced epochs than to the free epochs in both level and spatial shape. That matches the localizer, where breath-focused conditions (Intero, Match) share the same deactivation pattern (RSA below: Intero and Match are the closest pair in most regions). It does not mean free breathing is "less interoceptive".

**Reliability (the noise ceiling;** `EXPLORATORY_LOC_C3_T4_reliability.csv`).

| Measure | dot, own | corr, own |
|---|---|---|
| Template split-half: index from the run-1 vs run-2 localizer C1 map, across paced epochs within participant (41 participants, 567 epochs) | r .22 (Spearman–Brown .36) | r .31 (.47) |
| Participant mean paced index, breathing run 1 vs run 2 (42) | r .65 (.78) | r .49 (.66) |
| Participant aware − unaware, breathing run 1 vs run 2 (27 with both cells in both runs) | r .08 (.15) | r .12 (.22) |
| ICC (participant intercept share), paced epochs | .55 (SD between .149, within .135) | .50 |

- The participant's **average** expression of C1 is a reliable trait (cross-run r .65).
- The **epoch-to-epoch** variation that an awareness effect would have to ride on is mostly noise. The template itself is only modestly reliable (split-half r .22), as the L4 gate already showed for individual C1 maps.
- The participant-level aware − unaware difference does not replicate across runs (r .08).
- The aware − unaware CI excludes differences beyond −0.047 and +0.013 (about a third and a tenth of the within-participant epoch SD of .135). Smaller effects cannot be excluded with this index.

## Physiology-informed between-participant relations (§10.2; brief)

Stage C1 already related breath tracking to MAIA, behaviour, adherence and C1 in the fROIs and masks; those are not repeated. The four new relations (BH-FDR ×4; `EXPLORATORY_LOC_C3_physio_relations.csv`) are all null:

| Relation | n | ρ [95% CI] | BF10 | q |
|---|---|---|---|---|
| Tracking PLV × C1 SAL–DMN background-FC difference (stage C2) | 33 | −.00 [−.40, .38] | 0.14 | .98 |
| Tracking PLV × mean paced-epoch T4 index | 32 | .12 [−.25, .47] | 0.21 | .98 |
| Tracking PLV × aware − unaware T4 index | 30 | .12 [−.30, .49] | 0.76 | .98 |
| Breathing-run adherence × mean paced-epoch T4 index | 32 | .01 [−.35, .35] | 0.14 | .98 |

How well people track their breath is unrelated to the C2 coupling effect and to T4 expression. Beyond these, §10.2 is covered by stage C1.

## LOSO decoding of Intero vs Extero (`EXPLORATORY_LOC_C3_decoding.csv`)

Per-run 4 mm Intero and Extero maps (each vs Rest) of 48 participants (92 runs, 184 maps) were classified with a linear SVM (C = 1), leaving one participant out. Chance is .50. The p values are from 5,000 within-run label-swap permutations; the floor is .0002. CIs are participant-bootstrap.

| Mask | Raw PSC maps: balanced accuracy [95% CI], p | Pattern only (map mean and scale removed) | Null 95th percentile |
|---|---|---|---|
| Whole brain | .72 [.66, .78], .0002 | .71 [.65, .77], .0002 | .55–.56 |
| SAL | .70 [.65, .76], .0002 | .66 [.61, .72], .0002 | .55–.56 |
| SMN | .71 [.65, .77], .0002 | .73 [.68, .79], .0002 | .55 |
| DMN | .71 [.65, .77], .0002 | .67 [.62, .74], .0002 | .56 |
| Insula (AInsL ∪ AInsR ∪ PIns) | .61 [.56, .67], .001 | .60 [.55, .66], .005 | .56 |

- **Every one of the 10 tests survives Bonferroni ×10.** The largest corrected p is .046 (insula, pattern).
- **Decoding does not depend only on the global deactivation.** With each map's mean and scale removed, accuracy changes by at most .04.
- **The insula decodes worst** (.60–.61), although the mask means there differ strongly (L1). In the insula the condition difference is mainly a level shift that is shared across participants and carries little transferable spatial pattern.
- **Applied to other conditions** (post hoc reading, logged). The classifier trained on Intero vs Extero calls 49–58% of held-out Match maps "Intero" in the raw features (55–63% pattern). For comparison, it calls 25–36% of held-out Extero maps and only 9–26% of held-out Self maps "Intero".
  - So Match falls between the two training classes.
  - Self is classified with Extero, so the classifier is not simply separating "Intero" from "anything else".
  - The pre-set cross-decoding balanced accuracies (.58–.66, p ≤ .028) mix the Match share with the Extero accuracy. For that reason they are not read as generalisation to Match.

## Cross-run crossnobis RSA (descriptive; `EXPLORATORY_LOC_C3_rsa.csv`)

- **Data.** 44 participants. Ledoit–Wolf multivariate noise normalisation; distances are per voxel, and their expected value is 0 when two conditions do not differ. Rest enters as the zero pattern.
- **Model RDMs.** The breath, visual and self models are only mildly collinear: r(visual, self) .53, r(breath, visual) −.09, r(breath, self) −.17, VIF ≤ 1.43. All three nonetheless share the large Self and Rest distances.

**Distances.**
- **Intero vs Extero** is reliably above 0 in four of the five a-priori parcels: AInsL dz .61, p .0002; dACC .43, p .006; PIns .35, p .025; IEAT sphere .58, p .0004. It is not above 0 in **AInsR** (dz .07, p .63). In SAL, SMN, DMN and VIS it is dz .71–.77, p < .0001.
- **Intero vs Match is the smallest distance in seven of the nine a-priori parcels and networks** (AInsL, dACC, PIns, IEAT sphere, SAL, SMN, DMN). In AInsR, Intero vs Extero is smaller still (×100: 0.2 vs 0.9), and in VIS, Extero vs Match is (2.6 vs 3.4). In the parcels it is not distinguishable from 0 (dz .08–.20, p ≥ .18), except AInsR (.38, p .016). In the networks it is small but above 0 (dz .31–.57, p .0005–.049).
- **Distances involving Self or Rest are the largest,** especially in VIS (×100: 9.5–24.1, against 2.6–6.2 among Intero, Extero and Match).
- **The individual (circular) fROIs** in AInsL, dACC and PIns have *negative* Intero–Extero distances (dz −0.57 to −0.95). Those fROIs were selected on both runs' C1, so this pair is not interpretable there (logged).

**Model fits.** These are per-participant Pearson r between the 10 distances and each model, as a group mean.
- **Visual** fits every region; it is strongest in VIS (.65) and the networks (.44–.48).
- **Self** fits DMN (.39) and VIS (.53). It is weak or absent in AInsL (.03) and PIns (.01).
- **Breath** is positive in every a-priori parcel and network: .10 (AInsR) to .29 (SMN), all p < .05, VIS .07.
  - It is the best-fitting model in PIns (.23 vs visual .18) and ties with visual in AInsL (.18 vs .18).
- Because every model shares the Self and Rest pairs, these fits describe the geometry; they do not separate the explanations.
- **Reading:** in most of the insula and cingulate, Intero and Match form one cluster, Extero sits close by, and Self and Rest are far away. This is what a breath-focus dimension predicts, and it fits the T4 and decoding results.

## Insula posterior → anterior gradient of C1 (`EXPLORATORY_LOC_C3_gradient*.csv`)

**No linear gradient.** The per-participant slope of the joint 4 mm C1 on MNI y in the CAREN insula (1,400 voxels; N 48) is −0.003 PSC per 10 mm [−0.016, 0.010], dz −0.07, p .62 (sign-flip .63).

| Specification | Slope (PSC per 10 mm) [95% CI] | dz | p (sign-flip) |
|---|---|---|---|
| Primary, both hemispheres | −0.003 [−0.016, 0.010] | −0.07 | .63 |
| Left / right | −0.006 / −0.001 | −0.11 / −0.02 | .46 / .89 |
| + x and z covariates | −0.008 [−0.021, 0.006] | −0.17 | .25 |
| Harvard–Oxford insula | −0.007 [−0.020, 0.006] | −0.16 | .27 |
| 6 mm map | −0.003 [−0.014, 0.009] | −0.07 | .65 |
| Intero vs Rest (descriptive) | +0.047 [0.033, 0.061] | +0.96 | .0001 |
| Extero vs Rest (descriptive) | +0.050 [0.035, 0.065] | +0.95 | .0001 |

- Both conditions are progressively less deactivated (relative to Rest) towards the anterior insula, by the same amount, so their difference has no slope.
- The C1 profile is not monotonic (post hoc description of the pre-set profile table). It is least negative in the mid insula (−0.06 PSC at y −6) and most negative at both ends (−0.27 at y −30, −0.23 at y +30).
- **The projection onto a principal cortical gradient was not done:** no gradient map is available locally (logged).

## Inhale- vs exhale-locked responses by condition (descriptive; `EXPLORATORY_LOC_C3_respphase.csv`)

**Feasibility.** The HRF-convolved inhale- and exhale-onset model was not fitted. It failed the pre-set rule: r(inhale, exhale) .91–.92 and VIF 9.2–29.5, because at about 12 breaths/min the HRF removes the alternation. This is logged.

**Model used.** Per-condition, unconvolved cos/sin breath-phase regressors were added to the registered per-run model, at 6 mm, on mask means:
- 69 runs, 40 participants (usable belt, timing and polarity);
- Intero contrasts use 34 Intero-set participants;
- inhale − exhale = the mean signal over the inhale half-cycle minus the exhale half-cycle;
- BH-FDR within each family.

**Results.**
- **The signal is lower during inhalation than exhalation in the paced and pacer conditions.** This is strongest in Match: SM dz −1.54, COper −0.91, dACC −0.64, PIns −0.58, sphere −0.53, all q ≤ .009. It is also present in Extero (SM −0.69, q .001; dACC −0.43, AInsR −0.45, q ≤ .038) and in Intero in SM only (−0.77, q .001). Rest and Self show no inhale − exhale difference (q ≥ .15).
- **Condition differences.** Match differs from Extero in SM (−0.66, q .004) and COper (−0.55, q .018). Match differs from Intero in SM and COper (−0.55 and −0.56, q .022). **Intero vs Extero does not differ** in any mask (q ≥ .10).
- **Amplitude** of the phase-locked component is larger in Intero than in Rest in AInsL, AInsR, dACC, SM and the sphere (dz .48–.65, q .008–.029). It is larger in Match than in Extero in SM (.58, q .008) and PIns (.40, q .047).
- **Reading.** Unconvolved respiratory-phase BOLD is dominated by physiological artefact (chest-motion B0 changes, head motion, CSF pulsation). It scales with breathing depth, and depth is greatest in Match and Intero (stage C1). The pattern follows depth and pacing; it is strongest in sensorimotor cortex, and attending to the breath without pacing (Intero vs Extero) does not change it. These are **not** neural inhale- or exhale-locked responses, and the data cannot isolate such responses at this TR.

## Deviations (stage C3)

**Logged before running (one entry, 200a864).** Readings and choices:
- the T4 set (44);
- the dot product as primary and correlation as secondary;
- the outcome coding (wrong direction excluded, and added as a sensitivity);
- the added paced vs free check, run intercepts and reliability measures;
- the four §10.2 relations;
- the decoding, RSA and gradient specifications;
- the gradient projection not done;
- the respiratory-phase onset model infeasible, so the phase model was used.

**Post hoc (logged after the tables were seen):**
- the cross-decoding balanced accuracy is not read as generalisation; the shares are reported instead, and the figure shows them;
- the negative fROI Intero–Extero distances are not interpreted;
- the gradient profile description;
- two code fixes that did not change any value by choice: an empty F-test row, and a dtype error in the phase group step.

No registered result was changed.

# Exploratory: cortical gradients and the breath-tracking correlate (Q1, Q4)

**Both analyses are exploratory and post hoc.** Every earlier result had been seen. Each analysis was logged in `Preregistration/localizer_deviations.md` and committed (3e55705) before any alignment or brain–behaviour value was computed. A later row records the post hoc notes.
- Script: `Scripts/intero2023_loc_cortgrad.py`. It is separate from the stage-C3 insula slope in `intero2023_loc_gradient.py`.
- Tables: `Analysis/results/EXPLORATORY_LOC_cortgrad_*.csv` and `EXPLORATORY_LOC_Q4_*.csv`.
- Figure: `Figures/localizer/loc_cortgrad.png`.

**Downloaded maps** (kept outside the repo):
- Schaefer-400 7-network: MNI152NLin2009cAsym res-02 from templateflow, and fsLR 32k from CBIG (commit 634f676).
- Margulies et al. (2016) gradients G1–G3, fsLR 32k, via neuromaps 0.0.7 (OSF 4mw3a).
- Installed: neuromaps 0.0.7, brainspace 0.3.0, brainsmash 0.11.0, statsmodels 0.15.0.
- This reverses the stage-C3 note that the gradient projection could not be done.

## Q1: where the localizer maps sit on the principal cortical gradients (tier A)

**Method in brief.**
- **Maps.** Nine per-participant effect maps (joint 6 mm, PSC):
  - Intero, Extero, Match and Self vs Rest;
  - C1, C2, C4 and C5;
  - the breathing-run paced − free map (swc7p H1).
- **Parcellation.** Each map is averaged in 394 Schaefer-400 cortical parcels. Six temporal-pole parcels are dropped because less than half of each lies in the brain mask.
- **Gradients.** Each parcel takes the mean of Margulies G1 (unimodal → transmodal), G2 (somatomotor → visual) and G3 (default → frontoparietal/attention) over its fsLR vertices. Parcels correspond by label. Volumetric and surface centroids agree to a median of 3.8 mm.
- **Tests.**
  - *Group-map alignment:* Spearman ρ, tested against 10,000 Váša spins. BrainSMASH variogram surrogates are the sensitivity null. Holm ×9 within each gradient; G1 and G2 are primary and G3 descriptive.
  - *Condition differences:* a paired sign-flip test on each participant's Fisher-z alignment, for 4 contrasts × G1/G2, Holm ×8.

**Positive control: passed, narrowly.** Self vs Rest (N 54) aligns with the transmodal end of G1: ρ +.30, spin p .045, which meets the pre-set rule (ρ > 0, p < .05).
- Under the BrainSMASH null, the same alignment gives p .081.
- The control therefore holds with the spin null only. It is weak evidence that the pipeline recovers the expected Self pattern.
- In the profile, Self is the only map that rises above Rest at the transmodal end (top decile +0.045 PSC).

**Group-map alignment, Spearman ρ** (`EXPLORATORY_LOC_cortgrad_alignment.csv`). Each cell gives ρ, then spin p Holm ×9, then [BrainSMASH p Holm ×9]. Bold = p < .05 under both nulls.

| Map | N | G1 (primary) | G2 (primary) | G3 (descriptive) |
|---|---|---|---|---|
| Intero vs Rest | 48 | −.21, .31 [1.0] | +.12, 1.0 [1.0] | **+.42, .0009 [.0009]** |
| Extero vs Rest | 54 | +.17, .76 [1.0] | −.11, 1.0 [1.0] | **+.39, .0009 [.0009]** |
| Match vs Rest | 54 | −.12, .76 [1.0] | −.03, 1.0 [1.0] | **+.21, .036 [.021]** |
| Self vs Rest (C6) | 54 | +.30, .22 [.49] | −.24, 1.0 [1.0] | −.03, .73 [.71] |
| **C1 Intero − Extero** | 48 | **−.48, .002 [.004]** | +.43, .16 [.17] | +.15, .26 [.16] |
| C2 Intero − Self | 48 | −.31, .11 [.31] | +.34, .16 [.51] | **+.37, .0009 [.0009]** |
| **C4 Match − Extero** | 54 | **−.41, .009 [.037]** | +.18, 1.0 [1.0] | **−.32, .002 [.0009]** |
| C5 Match − Intero | 48 | +.21, .15 [1.0] | −.33, **.005** [.60] | **−.53, .0009 [.0009]** |
| Breathing paced − free | 47 | +.06, .76 [1.0] | −.09, 1.0 [1.0] | **+.32, .0009 [.0012]** |

**Condition differences in participant alignment** (`EXPLORATORY_LOC_cortgrad_paired.csv`). Δ = mean difference in Fisher z, with a 95% participant-bootstrap CI.

| Contrast | N | G1: Δ [CI], dz, p Holm ×8 | G2: Δ [CI], dz, p Holm ×8 |
|---|---|---|---|
| C1 Intero − Extero | 48 | −0.18 [−0.24, −0.12], −0.85, **.0008** | +0.16 [0.10, 0.21], +0.79, **.0008** |
| C2 Intero − Self | 48 | −0.24 [−0.33, −0.15], −0.78, **.0008** | +0.24 [0.18, 0.30], +1.10, **.0008** |
| C4 Match − Extero | 54 | −0.14 [−0.21, −0.07], −0.55, **.0012** | +0.05 [−0.01, 0.11], +0.24, .18 |
| C5 Match − Intero | 48 | +0.05 [−0.01, 0.11], +0.22, .18 | −0.09 [−0.15, −0.03], −0.42, **.013** |

The participant alignments themselves are weak. Mean z lies between −0.16 and +0.14 on G1 and G2, so a single participant's map follows a gradient only loosely.

**What the profiles show** (`EXPLORATORY_LOC_cortgrad_profiles.csv`; PSC by gradient decile).
- **G1.**
  - Intero and Match are deactivated in every decile, and most deeply at the transmodal end: Intero from −0.23 (lowest decile) to −0.30 (top); Match from −0.24 to −0.30.
  - Extero is flat (−0.10 at both ends).
  - So C1 is negative everywhere and most negative in the transmodal decile: −0.19, against −0.05 to −0.13 in deciles 1–8.
  - **Reading.** The breath conditions deactivate the default/transmodal end more than the visual task does. Shifting attention to the breath moves the pattern away from the transmodal pole relative to Extero and Self. It does not move it toward that pole.
- **G2.**
  - Intero is flat along G2 (−0.19 to −0.28).
  - Extero is most deactivated at the visual pole (−0.18 in the top decile, matching the calcarine and precuneus negatives of stage A).
  - C1 is therefore near 0 at the visual pole (−0.01) and negative elsewhere.
  - **Reading.** The positive C1 × G2 shift reflects Extero's visual-pole deactivation. It does not reflect an Intero shift toward the somatomotor pole.
  - Match is somewhat more deactivated at the visual pole than Intero (C5 G2). That alignment is not robust to BrainSMASH (Holm .60).
- **G3** (descriptive). Every task condition and the breathing-run paced − free map are least deactivated at the frontoparietal/attention end (Intero −0.32 in the lowest decile vs −0.13 in the top). They are most deactivated at the default and sensorimotor end. Relative to Rest, the attention networks are the least suppressed.
- **Breathing runs.** Paced − free has no G1 or G2 alignment (ρ +.06, −.09; spin p .65, .60), only the shared G3 pattern. The swc7p deactivation is spread evenly along the hierarchy.

**Against the plan's expectations.**
- The plan's informative outcome was "Intero and Match shifted toward the somatomotor pole of G2 relative to Extero, with no G1 difference". It was **not** observed: the C1 difference is on G1, and the G2 effect comes from Extero. This reading is logged as post hoc.
- Breath attention does separate from self-reference along G1 (C2 G1 dz −0.78). But it is not "unrelated to G1": Intero sits lower on G1 than Extero too.
- As the plan expected, the gradients mainly restate stage A: the breath-focused deactivation is widespread, and deepest in sensorimotor and transmodal cortex.

## Q4: does breath tracking have a brain correlate beyond pacer-following and breathing regularity?

**What this can and cannot test.** There is no valid circle-tracking score. The Extero / Match logs keep only the last key in each window, and the scanner trigger usually overwrites the press (deviations log, 2026-11-28). So this analysis **does not isolate a "breath-specific" tracking effect beyond visual tracking.** It asks whether breath-tracking accuracy relates to the brain response once pacer-following (Match PLV), breathing regularity (Intero `period_cv`) and protocol version are held constant.

**Model.**
- One FLAME 1 model per map: C1 and Intero vs Rest, joint 6 mm.
- N 34, the Intero set ∩ all covariates (4 on v1).
- Covariates, all demeaned: tracking PLV (median .80, range .40–.94), Match PLV, `period_cv` and protocol.
- VIF ≤ 2.83, so Match PLV was kept.
- Tracking term tested two-sided:
  - voxel FWE in the five a-priori parcels, ×5;
  - whole-brain cluster FWE at z > 3.09.
- **Check.** The zero-order Spearman ρ between tracking and the C1 mask means (−.05 to −.25) reproduces the stage-C1 rows for the same 34.

**Result: null everywhere** (`EXPLORATORY_LOC_Q4_roi.csv`, `_maskmeans.csv`, `_wholebrain.csv`).

| Parcel | C1: peak z, p ×5 | C1 mask mean: partial r, p (perm.), BF01 | Intero vs Rest: peak z, p ×5 | Intero mask mean: partial r, p (perm.), BF01 |
|---|---|---|---|---|
| AInsL | −1.89, 1.0 | −.11, .55, 2.2 | +1.92, 1.0 | +.06, .75, 2.5 |
| AInsR | −2.14, 1.0 | −.19, .29, 1.6 | +1.75, 1.0 | −.02, .92, 2.6 |
| dACC | −2.53, 1.0 | −.18, .35, 1.7 | +2.32, 1.0 | +.06, .75, 2.5 |
| PIns | −1.97, 1.0 | −.18, .34, 1.7 | +2.81, .59 | +.02, .90, 2.6 |
| IEAT sphere | −2.53, 1.0 | −.21, .25, 1.4 | −1.83, 1.0 | +.01, .97, 2.6 |

- No voxel survives in any parcel; the smallest corrected p is .59. No whole-brain cluster survives either: the largest is 24 voxels, cluster p_FWE 1.0.
- The mask-mean Bayes factors lean toward the null but are weak:
  - BF01 1.4–2.2 for C1;
  - BF01 2.5–2.6 for Intero vs Rest;
  - JZS r = .707. At r = .354, BF10 is 0.66–0.92 for C1 and 0.57–0.59 for Intero.
- Adjusting for the covariates does not reveal a tracking correlate hidden by suppression. The C1 partial correlations (−.11 to −.21) are of the same size and sign as the zero-order ones.
- **Power.** The adjusted tracking score has a reliability of about .5 at best. With N 34, the test could detect only observed |r| ≳ .46. This result therefore cannot exclude small true relations; it adds no positive evidence for one.

## Deviations (Q1, Q4)
- **Logged before running (3e55705):** both full specifications, and the label-free checks (install, downloads, parcel coverage, label and centroid correspondence, gradient orientation, Q4 N and VIFs).
- **Post hoc (logged after the tables were seen):**
  - a figure-code fix (no value changed);
  - the positive control's narrow margin, and its p .081 under BrainSMASH;
  - the reading that the plan's G2 pattern was not observed;
  - the C5 G2 alignment marked as not robust.
- **Not run:** tier B (participant gradients); a surface refit of the first levels (the plan's sensitivity).
- No registered result was changed.

# Exploratory: does the breathing-depth difference explain the Intero deactivation across participants?

**Exploratory and post hoc.** Every earlier localizer result had been seen. The full specification was logged in `Preregistration/localizer_deviations.md` and committed (aafc71b) before any intercept or slope was computed. A later row records the post hoc notes.
- Script: `Scripts/intero2023_loc_depth_c1.py`.
- Tables: `Analysis/results/EXPLORATORY_LOC_depth_{checks,models,subjects}.csv`.
- Figure: `Figures/localizer/loc_depth_c1.png`.

**Question.** Intero breathing is deeper than Extero breathing (RVT dz .88). Within-run RVT regressors remove 16–25% of L1. Here the question is asked between participants. Does a participant's Intero − Extero depth difference predict their C1? And what C1 does the regression predict at equal depth (x = 0, the intercept)?

**Method in brief.**
- **Depth measure, x.** It is scale-free, because belt gain differs between participants. Per run, x is the log of mean block RVT in Intero over mean block RVT in Extero. Block RVT = median breath amplitude / median breath period, from the step-1 block table. The participant value is the mean over usable runs, and zero means equal depth.
  - The stored `rvt_z` is z-scored within run, so it cannot form a ratio.
  - x correlates .75 with the Intero − Extero `rvt_z` difference.
  - Secondary measure: the same ratio using breath amplitude only (depth).
- **Outcome.** The registered joint 6 mm C1 mask mean in AInsL, AInsR and dACC. The extraction reproduces `LOC_L1_subjects.csv` exactly. Secondary outcome: the cortical grey-matter mean.
- **Model.** OLS, C1 ~ 1 + x, with x uncentred and HC3 SEs. The primary test is the intercept, Holm ×3. The slope is tested Holm ×3 separately.
- **Measurement-error correction (regression calibration).** The true-score variance is the run 1 × run 2 covariance of x. Each participant's reliability is the Spearman–Brown value for their number of runs. CIs come from a participant bootstrap that re-estimates the reliability in every draw (10,000 draws, seed 20261128).

**N.** 31 of the 48 L1 participants (22 with both runs).
- 13 have no usable belt run.
- 4 have no Extero block with ≥ 2 detected breaths.
- The 31 have mean C1 −0.115 / −0.084 / −0.131 PSC (dz −0.72 / −0.46 / −0.62). This is close to the full L1 sample.

**The depth measure is noisy.**
- Intero is deeper in 23 of 31 participants: mean x 0.53, i.e. about 1.7 × the Extero RVT, SD 0.77. x = 0 lies 0.7 SD below the mean, inside the data (8 participants have x ≤ 0).
- Run 1 vs run 2: r .21 (95% CI −.23 to .58; n 22). Spearman–Brown reliability is .35 for a two-run mean and .21 for one run. The depth-only ratio gives .24 / .39.
- The registered `rvt_z` difference is not reliable at all (r −.05).
- One-run participants vary twice as much as single runs of two-run participants (var 1.10 vs 0.55). Their single runs are probably noisier still.

**C1 on the depth ratio (N 31; `EXPLORATORY_LOC_depth_models.csv`).** CIs are 95% HC3.

| Mask | Mean C1 | Intercept (C1 at equal depth) [95% CI] | p, Holm ×3 | Intercept / mean | Slope per log unit [95% CI] | p, Holm ×3 |
|---|---|---|---|---|---|---|
| AInsL | −0.115 | −0.076 [−0.148, −0.004] | .039, **.12** | 0.66 | −0.073 [−0.125, −0.020] | .009, **.026** |
| AInsR | −0.084 | −0.042 [−0.128, +0.044] | .32, **.32** | 0.50 | −0.078 [−0.143, −0.014] | .019, **.037** |
| dACC | −0.131 | −0.087 [−0.183, +0.010] | .076, **.15** | 0.66 | −0.083 [−0.152, −0.014] | .020, **.037** |
| Cortical GM (secondary) | −0.105 | −0.069 [−0.131, −0.008] | .028 (unadjusted) | 0.66 | −0.066 [−0.107, −0.026] | .002 |

- **Slope.** A larger Intero − Extero depth difference goes with a larger deactivation, in every mask (Holm p .026–.037; Spearman ρ −.42 to −.46). The same holds in cortical grey matter.
- **Intercept.** At equal depth the predicted C1 is 50–66% of the mean and still negative in every mask. No mask passes Holm. AInsL is the only mask with an uncorrected p below .05 (.039).
- **Influence.** Leave-one-out intercepts range from −0.092 to −0.060 (AInsL), −0.058 to −0.016 (AInsR) and −0.100 to −0.060 (dACC). The maximum leverage is .22.
- **Depth only (secondary).** Intercepts −0.068 / −0.034 / −0.079 (Holm .12 / .39 / .16). Slopes −0.096 / −0.102 / −0.107 (Holm .003 / .006 / .005).
- **Two-run participants only (sensitivity, n 22).** Intercepts −0.051 / −0.020 / −0.045 (p .23 / .72 / .44). Slopes −0.117 / −0.130 / −0.132 (p .029 / .082 / .046).

**Corrected for the unreliability of x.**

| Mask | Corrected intercept [bootstrap 95% CI] | Corrected slope [CI] |
|---|---|---|
| AInsL | +0.034 [−0.103, +0.724] | −0.29 [−1.66, −0.08] |
| AInsR | +0.078 [−0.085, +0.918] | −0.31 [−2.01, −0.07] |
| dACC | +0.041 [−0.146, +0.905] | −0.33 [−2.12, −0.06] |
| Cortical GM | +0.027 [−0.084, +0.619] | −0.25 [−1.51, −0.08] |

- Correcting the attenuated slope moves the point intercepts to about zero or slightly positive.
- The correction is barely identified. In 3,065 of 10,000 bootstrap draws the reliability fell below the 0.10 floor. The upper CI bounds are therefore set largely by the floor.
- The corrected values bound the problem; they are not an estimate. The uncorrected intercept overstates the C1 that remains at equal depth by an amount these data cannot fix. The plausible range runs from about the uncorrected value to zero or above.

**Match − Extero (C4; N 37, 24 with both runs; Match/Extero RVT ratio, mean 0.78).**
- C4 at equal depth stays clearly negative: intercepts −0.157 [−0.237, −0.077], −0.176 [−0.280, −0.071] and −0.178 [−0.274, −0.082] (Holm p .001 / .002 / .001). That is 75–87% of the mean C4.
- The slopes are smaller: −0.067 / −0.035 / −0.072 (Holm .11 / .40 / .11).
- Depth only: intercepts Holm ≤ .001 in all three masks.
- Two-run only (n 24): intercepts −0.099 / −0.140 / −0.109 (p .041 / .068 / .031).
- Corrected intercepts are −0.058 [−0.192, +0.879], −0.122 [−0.351, +0.571] and −0.078 [−0.245, +1.003]; 1,054 draws were floored.
- C4 also contrasts breath control with visual tracking, and its depth difference is larger. Yet more of it remains at equal depth than for C1.

**PetCO2 (descriptive, no correction).**
- **C1 subset (n 20).** End-tidal CO2 is *higher* in Intero than in Extero blocks: +0.28 %CO2, dz 0.63, p .011. It is unrelated to x (ρ .10).
  - C1 ~ ΔPetCO2: slopes −0.16 / −0.23 / −0.19 (p .07 / .04 / .12).
  - Intercepts −0.081 / −0.017 / −0.077 (p .05 / .65 / .15).
- **C4 subset (n 23).** ΔPetCO2 +0.40 %CO2, dz 0.71, p .003, ρ with x .71.
  - C4 intercepts at ΔPetCO2 = 0 stay negative (p .001 / .015 / .004).
- **Reading.** Deeper breaths sample more alveolar gas, so a higher end-tidal reading need not mean higher arterial CO2. In the small C1 subset, larger CO2 differences go with larger deactivations. That fits a CO2 contribution, but it is not a test of one.

**Reading.**
- **The slope.** Across participants, a larger Intero − Extero breathing-depth difference goes with a larger C1 deactivation in all three registered masks. The between-participant link between depth and the deactivation is real.
- **The intercept.** At equal depth the predicted C1 shrinks to 50–66% of its mean. It is still negative in every mask, but not Holm-significant in any. The depth ratio is so unreliable (two-run reliability about .35) that even this intercept is probably an overstatement. The corrected intercept is near zero, with a CI too wide to say anything.
- **What can be concluded.** These data do not show that the Intero deactivation survives at equal breathing depth. They do not show that it disappears either.
- **Two caveats limit any causal reading:**
  - deeper breathing may itself be part of attending to the breath;
  - a between-participant slope can also reflect shared engagement (people who engage more may both breathe more deeply and deactivate more).
- **How this fits the earlier results.** The within-run RVT co-report (L1 breathing-robust, 75–84% retained) remains the registered answer. This between-participant analysis is weaker evidence for a non-respiratory component. The Match − Extero contrast keeps most of its effect at equal depth.

## Deviations (depth × C1)
- **Logged before running (aafc71b):** the full specification, with the label-free checks (coverage, distributions, run 1 vs run 2 reliability).
- **Post hoc (logged after the tables were seen):**
  - the reading of the corrected estimates as uninformative, because 31% of bootstrap draws (C1) hit the reliability floor;
  - a figure legend wording change.
- No registered result was changed.

# Synthesis: what the localizer shows about attention to the breath

**What it shows.**
- **Whole-brain deactivation.** Directing attention to one's own breath, rather than to a visual pacer, produces a large, reproducible and mostly whole-brain *deactivation*. It is strongest in sensorimotor, insular and cingulate cortex, and it largely survives RVT and PetCO2 correction (L1; physio-corrected maps).
- **A distributed pattern.** The difference is a distributed pattern, not only a level shift. It is decodable across participants from whole brain and from each network after the map mean is removed (balanced accuracy .66–.73), and it is reliable within participants at the mask level (crossnobis).
- **A breath-focus dimension.** Paced breath-matching (Match) shares that pattern. Intero and Match are the closest pair in most regions, a classifier puts Match between Intero and Extero, and breathing-run paced epochs express the localizer C1 pattern more than free epochs do.
- **DMN coupling.** At the network level, breath attention comes with the DMN being less anticorrelated with the salience network (C2).

**What it does not show.**
- **No reliable participant-specific region.** It does not identify a reliable participant-specific "interoceptive attention" region: the top C1 voxels do not replicate across runs (L4).
- **No awareness signal.** It gives no evidence that this pattern carries awareness of breathing changes:
  - T1 and T2 were null in the fROIs;
  - T4 was null at the level of single epochs, for outcome and confidence alike, with a reliable participant-level index but a noisy epoch-level one;
  - tracking accuracy is unrelated to any brain measure.
- **No IEAT replication or spatial gradient.** It does not replicate the IEAT MAIA effect (L3), and C1 has no posterior → anterior insula gradient.
- **No neural respiratory-phase responses.** Respiratory-phase-locked signals behave like breathing-depth artefact, not like attention-dependent responses.
- **Attention is not isolated.** The contrasts cannot fully separate attention from breathing depth, agency and press predictability (plan §11). The deactivation is the localizer's most robust neural correlate of attending to the breath, and its role in breath *awareness* remains untested by any positive result here.
