#!/usr/bin/env python
"""
intero2023_h5_reliability.py -- H5 equivalence (TOST) and the registered run-1 vs
run-2 reliability diagnostic (prereg v0.8 §3 H5, §12.4, §8).

Run after intero2023_group.py has produced the base and corrected ROI tests in
each physiological subset (it reads Analysis/results/*_roi.csv for the grid's
"corrected effect remains significant" cell).

H5: per participant, the ROI mean of the contrast from the base and the corrected
model; TOST (alpha .05 per side) on the equal-weighted paired difference with
bounds in dz units (raw bound = dz bound x SD of the differences).
  H3 (confirmatory if H3 is supported): C2joint; PetCO2 S1 N = 25 +/-0.61;
     Biopac RVT S2 N = 28 +/-0.58; S2b N = 39 +/-0.49 (robustness)
  H2 (supportive): C3; S1 25 +/-0.61; S2 27 +/-0.59; S2b 38 +/-0.49
Reliability: C1 (Model A) and C3 (Model B), two-run participants: within-
participant spatial correlation of the run-1 and run-2 effect maps (group mask),
and the cross-participant ICC(3,1) and Pearson r of ROI means. Descriptive.
Outputs: Analysis/results/H5_tost.csv, Analysis/results/reliability.csv
"""
import os, json
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
import intero2023_glm_common as g
import intero2023_group as G

RD = os.path.join(g.REPO, "Analysis", "results")
ROIS3 = ["caren_AInsL", "caren_AInsR", "caren_dACC"]
BOUNDS = {("H3", "co2", "S1"): (0.61, 25), ("H3", "rvt", "S2"): (0.58, 28), ("H3", "rvt", "S2b"): (0.49, 39),
          ("H2", "co2", "S1"): (0.61, 25), ("H2", "rvt", "S2"): (0.59, 27), ("H2", "rvt", "S2b"): (0.49, 38)}


def roi_vals(files, roi):
    mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
    return np.array([mk.transform(f).mean() for f in files])


def tost(d, b):
    n = len(d); m = d.mean(); sd = d.std(ddof=1); se = sd / np.sqrt(n); D = b * sd
    p_lo = stats.t.sf((m + D) / se, n - 1); p_hi = stats.t.cdf((m - D) / se, n - 1)
    return dict(mean_diff=m, sd_diff=sd, dz=m / sd, bound_dz=b, bound_raw=D,
                p_lower=p_lo, p_upper=p_hi, p_tost=max(p_lo, p_hi))


def sig_corrected(target, var, subset, roi):
    f = os.path.join(RD, f"{target}_var-{var}_set-{subset}_roi.csv")
    d = pd.read_csv(f)
    if "test" in d:
        d = d[d.test == "predicted"]
    return float(d[d.roi == roi].p_fwe_bonf.iloc[0])


def supported(target):
    d = pd.read_csv(os.path.join(RD, f"{target}_var-base_set-all_roi.csv"))
    if "test" in d:
        d = d[d.test == "predicted"]
    return sorted(d[d.p_fwe_bonf < 0.05].roi)


def h5():
    sets, _, _ = G.subject_sets()
    rows = []
    sup = {t: supported(t) for t in ("H3", "H2")}
    for (target, var, subset), (b, n_reg) in BOUNDS.items():
        con = G.CON[target]
        subs = G.subjects_for(target, var, subset, sets)
        subs = [s for s in subs if os.path.exists(G.fl(s, "base", con, "effect"))]
        assert len(subs) == n_reg, (target, var, subset, len(subs), n_reg)
        for roi in ROIS3:
            y0 = roi_vals([G.fl(s, "base", con, "effect") for s in subs], roi)
            y1 = roi_vals([G.fl(s, var, con, "effect") for s in subs], roi)
            r = tost(y1 - y0, b)
            p_corr = sig_corrected(target, var, subset, roi)
            cell = ("statistically equivalent" if r["p_tost"] < 0.05 else
                    "attenuated but present" if p_corr < 0.05 else "indeterminate")
            rows.append(dict(target=target, correction=var, subset=subset, roi=roi, N=len(subs),
                             robustness_only=subset == "S2b", target_supported_in_roi=roi in sup[target],
                             r_corrected_uncorrected=float(np.corrcoef(y0, y1)[0, 1]),
                             mean_uncorrected=float(y0.mean()), mean_corrected=float(y1.mean()),
                             p_corrected_effect_bonf=p_corr, grid_cell=cell, **r))
    d = pd.DataFrame(rows)
    # registered decision rule (prereg §3 H5)
    if sup["H3"]:
        rule = "H3 supported: H5 verdict is the H3 TOST cell in the ROI(s) where H3 was supported; H2 supportive"
    elif sup["H2"]:
        rule = "H3 not supported, H2 supported: H5 verdict is the H2 TOST cell, about the confidence measure only"
    else:
        rule = "Neither H3 nor H2 supported: H5 not tested; both comparisons descriptive"
    d["decision_rule"] = rule
    d.to_csv(os.path.join(RD, "H5_tost.csv"), index=False)
    print(rule)


def icc31(a, b):
    X = np.column_stack([a, b]); n, k = X.shape
    ms_r = k * ((X.mean(1) - X.mean()) ** 2).sum() / (n - 1)
    ms_e = (((X - X.mean(1, keepdims=True) - X.mean(0, keepdims=True) + X.mean()) ** 2).sum()) / ((n - 1) * (k - 1))
    return (ms_r - ms_e) / (ms_r + (k - 1) * ms_e)


def reliability():
    man = g.run_list()
    two = man.groupby("subject").size().pipe(lambda x: x[x == 2]).index
    gm = NiftiMasker(mask_img=g.GROUP_MASK).fit()
    rows = []
    for con, rois in (("C1", G.SETS["H1"]), ("C3", G.SETS["H2"])):
        subs = [s for s in two if all(os.path.exists(os.path.join(g.OUT, "firstlevel", s,
                f"{s}_run-{k}_var-base_con-{con}_effect.nii.gz")) for k in (1, 2))]
        f = lambda s, k: os.path.join(g.OUT, "firstlevel", s, f"{s}_run-{k}_var-base_con-{con}_effect.nii.gz")
        sp = [np.corrcoef(gm.transform(f(s, 1)).ravel(), gm.transform(f(s, 2)).ravel())[0, 1] for s in subs]
        rows.append(dict(contrast=con, measure="spatial_r_within_participant", roi="group_mask", N=len(subs),
                         median=float(np.median(sp)), q25=float(np.percentile(sp, 25)), q75=float(np.percentile(sp, 75))))
        for roi in rois:
            a = roi_vals([f(s, 1) for s in subs], roi); b = roi_vals([f(s, 2) for s in subs], roi)
            rows.append(dict(contrast=con, measure="roi_mean_across_participants", roi=roi, N=len(subs),
                             icc31=float(icc31(a, b)), pearson_r=float(np.corrcoef(a, b)[0, 1])))
    pd.DataFrame(rows).to_csv(os.path.join(RD, "reliability.csv"), index=False)
    print("wrote reliability.csv")


if __name__ == "__main__":
    import sys
    {"h5": h5, "reliability": reliability}[sys.argv[1]]()
