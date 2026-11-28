#!/usr/bin/env python
"""
intero2023_loc_gradient.py -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §10.5):
posterior -> anterior gradient of C1 (Intero - Extero) along the insula, as a per-participant
slope on MNI y. Everything this script writes is exploratory.

Choices: Analysis/localizer_stageC3_implementation_notes.md §5 (committed before any value).

  python Scripts/intero2023_loc_gradient.py

Set: localizer Intero set (48). Map: joint 4 mm C1 effect (PSC). Region: CAREN AInsL u AInsR u
PIns (PIns split by the sign of x for per-hemisphere slopes). Slope: OLS of voxel C1 on y, PSC
per 10 mm (positive = C1 larger / less negative anteriorly). Group: one-sample t (two-sided), dz,
95% CI, sign-flip permutation p (10,000, seed 20261126). Sensitivity: + x and z covariates;
Harvard-Oxford insula; 6 mm map. Descriptive: Intero vs Rest and Extero vs Rest slopes; y-bin profile.
The projection onto a principal cortical gradient is not done (no gradient map available locally).
Tables: Analysis/results/EXPLORATORY_LOC_C3_gradient*.csv
"""
import os, sys, warnings
import numpy as np, pandas as pd, nibabel as nib
import nibabel.affines  # noqa: F401
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
N_PERM, SEED = 10000, 20261126


def region(names):
    m = None
    for n in names:
        v = nib.load(os.path.join(g.REPO, "ROIs", f"{n}.nii")).get_fdata() > 0
        m = v if m is None else (m | v)
    return m


def coords(m, aff):
    ijk = np.argwhere(m)          # C order, same as d[m]
    return nib.affines.apply_affine(aff, ijk)


def slope(y, xyz, cov=False):
    Y = xyz[:, 1] / 10.0
    X = np.c_[np.ones(len(Y)), Y] if not cov else np.c_[np.ones(len(Y)), Y, xyz[:, 0] / 10, xyz[:, 2] / 10]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b[1]


def signflip_p(v, rng):
    v = np.asarray(v); obs = abs(v.mean())
    null = np.abs((rng.choice([-1, 1], size=(N_PERM, len(v))) * v).mean(1))
    return (1 + (null >= obs - 1e-15).sum()) / (N_PERM + 1)


def main():
    s, _, _ = L.participant_sets()
    subs = s["intero"]
    aff = nib.load(g.GROUP_MASK).affine
    REG = {"caren": ["caren_AInsL", "caren_AInsR", "caren_PIns"], "ho": ["ho_AInsL", "ho_AInsR", "ho_PIns"]}
    specs = [  # name, region, smoothing, contrast, covariates, hemisphere
        ("primary", "caren", 4, "C1", False, "both"),
        ("left", "caren", 4, "C1", False, "left"),
        ("right", "caren", 4, "C1", False, "right"),
        ("covariates_xz", "caren", 4, "C1", True, "both"),
        ("harvard_oxford", "ho", 4, "C1", False, "both"),
        ("smooth_6mm", "caren", 6, "C1", False, "both"),
        ("intero_vs_rest", "caren", 4, "intero", False, "both"),
        ("extero_vs_rest", "caren", 4, "extero", False, "both"),
    ]
    masks = {k: region(v) for k, v in REG.items()}
    XYZ = {k: coords(m, aff) for k, m in masks.items()}
    cache, rows, prof = {}, [], []
    for x in subs:
        for name, reg, sm, con, cov, hemi in specs:
            key = (x, reg, sm, con)
            if key not in cache:
                cache[key] = nib.load(L.fl_name(x, "joint", sm, "base", con, "effect")).get_fdata(dtype=np.float32)[masks[reg]].astype(float)
            y, xyz = cache[key], XYZ[reg]
            sel = np.ones(len(y), bool) if hemi == "both" else (xyz[:, 0] < 0 if hemi == "left" else xyz[:, 0] > 0)
            rows.append(dict(subject=x, spec=name, slope_per10mm=slope(y[sel], xyz[sel], cov), mean=y[sel].mean(),
                             n_vox=int(sel.sum())))
        y = cache[(x, "caren", 4, "C1")]; yy = XYZ["caren"][:, 1]
        bins = np.arange(np.floor(yy.min() / 4) * 4, yy.max() + 4, 4)
        bi = np.digitize(yy, bins)
        for b in np.unique(bi):
            prof.append(dict(subject=x, y_bin_center=bins[b - 1] + 2, C1=y[bi == b].mean(), n_vox=int((bi == b).sum())))
    S = pd.DataFrame(rows)
    rng = np.random.default_rng(SEED)
    out = []
    for name, d in S.groupby("spec", sort=False):
        v = d.slope_per10mm.values; t = stats.ttest_1samp(v, 0)
        se = v.std(ddof=1) / np.sqrt(len(v)); tc = stats.t.ppf(.975, len(v) - 1)
        out.append(dict(spec=name, N=len(v), n_vox=int(d.n_vox.iloc[0]), slope_mean=v.mean(), slope_sd=v.std(ddof=1),
                        ci_lo=v.mean() - tc * se, ci_hi=v.mean() + tc * se, dz=v.mean() / v.std(ddof=1),
                        t=t.statistic, p=t.pvalue, p_signflip=signflip_p(v, rng), share_positive=float((v > 0).mean()),
                        region_mean=d["mean"].mean()))
    T = pd.DataFrame(out)
    T.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient.csv"), index=False)
    S.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient_subjects.csv"), index=False)
    P = pd.DataFrame(prof).groupby("y_bin_center").agg(C1_mean=("C1", "mean"), C1_sd=("C1", "std"),
                                                       N=("C1", "size"), n_vox=("n_vox", "first")).reset_index()
    P["C1_se"] = P.C1_sd / np.sqrt(P.N)
    P.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient_profile.csv"), index=False)
    print(T.round(4).to_string()); print(P.round(3).to_string())


if __name__ == "__main__":
    main()
