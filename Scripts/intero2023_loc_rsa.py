#!/usr/bin/env python
"""
intero2023_loc_rsa.py -- localizer stage C3 (EXPLORATORY, descriptive; OSF p5ztc plan §10.4):
cross-run crossnobis RSA of the five localizer conditions (Intero, Extero, Match, Self, and the
Rest baseline as the zero pattern). Everything this script writes is exploratory.

Choices: Analysis/localizer_stageC3_implementation_notes.md §4 (committed before any value).

  python Scripts/intero2023_loc_rsa.py subject <sub>    # per-run refit in ROI voxels -> RDMs
  python Scripts/intero2023_loc_rsa.py group            # tables

Per run: the registered per-run 4 mm localizer model (intero2023_loc_common.run_designs, base,
200 s DCT; AR(1) via nilearn run_glm) refitted in the union of the ROI voxels; betas of intero,
extero, match, self; whitened residuals. Noise normalisation per run: Ledoit-Wolf shrinkage of the
residual covariance, applied as Sigma^(-1/2) through the residual SVD (exact). Crossnobis
d(i, j) = (b_i^1 - b_j^1)'(b_i^2 - b_j^2) / P.
Outputs: $OUT/loc_c3/rsa/<sub>.json; Analysis/results/EXPLORATORY_LOC_C3_rsa*.csv
"""
import os, sys, json, itertools, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
OD = os.path.join(g.OUT, "loc_c3", "rsa")
FD = os.path.join(g.OUT, "froi")
COND = ["intero", "extero", "match", "self", "rest"]
LAB = ["I", "E", "M", "S", "R"]
PAIRS = list(itertools.combinations(range(5), 2))
PARCELS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere"]
NETS = ["net_SAL", "net_SMN", "net_DMN", "net_VIS"]
PRIMARY_FROI = {"caren_AInsL": "C1", "caren_AInsR": "C1loso", "caren_dACC": "C1",
                "caren_PIns": "C1", "ieat_ACCsphere": "C1loso"}      # stage-B gate result


def rmask(n):
    return nib.load(os.path.join(g.REPO, "ROIs", f"{n}.nii")).get_fdata() > 0


def roi_sets(sub):
    """ROI name -> boolean 3-D mask."""
    R = {p: rmask(p) for p in PARCELS + NETS}
    z = np.load(os.path.join(FD, "froi_masks.npz"))
    for p, rule in PRIMARY_FROI.items():
        sel = z[f"{sub}|{rule}|0.10|{p}"]
        m = np.zeros(R[p].shape, bool); m[R[p]] = sel
        R[f"froi_{p}"] = m
    return R


def model_rdms():
    breath = np.array([1, 0, 1, 0, 0]); vis = np.array([0, 0, 0, 1, 2]); slf = np.array([0, 0, 0, 1, 0])
    v = lambda f: np.array([float(f(i, j)) for i, j in PAIRS])
    return {"breath": v(lambda i, j: breath[i] != breath[j]),
            "visual": v(lambda i, j: vis[i] != vis[j]),
            "self": v(lambda i, j: slf[i] != slf[j])}


def whiten(B, R):
    """B: k x P betas, R: n x P residuals. Returns B Sigma^(-1/2), Sigma = LW-shrunk residual covariance."""
    from sklearn.covariance import ledoit_wolf_shrinkage
    Rc = R - R.mean(0)
    n, P = Rc.shape
    a = ledoit_wolf_shrinkage(Rc, assume_centered=True, block_size=2000)
    U, s, Vt = np.linalg.svd(Rc / np.sqrt(n), full_matrices=False)
    mu = (s ** 2).sum() / P
    lam = (1 - a) * s ** 2 + a * mu
    BV = B @ Vt.T
    return (B - BV @ Vt) / np.sqrt(a * mu) + (BV / np.sqrt(lam)) @ Vt, float(a)


def run_fit(sub, run_row, vox):
    """Betas (4 x V) and whitened residuals (n x V) of the registered per-run model in `vox` voxels."""
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level import run_glm
    from nilearn.glm.first_level.first_level import mean_scaling
    ref = nib.load(g.GROUP_MASK)
    masker = NiftiMasker(mask_img=nib.Nifti1Image(vox.astype(np.uint8), ref.affine),
                         smoothing_fwhm=4, standardize=False).fit()
    n = int(run_row.n_vol)
    Y = masker.transform(g.bold_path(sub, L.entity(run_row.run))); assert Y.shape[0] == n
    Y, _ = mean_scaling(Y, axis=0)
    (_, dm), = L.run_designs(sub, pd.DataFrame([run_row]))
    X = dm.values.astype(float)
    labels, results = run_glm(Y, X, noise_model="ar1")
    cols = list(dm.columns)
    B = np.zeros((4, Y.shape[1])); Rw = np.zeros_like(Y, dtype=np.float64)
    for lab, res in results.items():
        sel = labels == lab
        B[:, sel] = res.theta[[cols.index(c) for c in ("intero", "extero", "match", "self")]]
        Rw[:, sel] = res.whitened_residuals
    return B, Rw


def subject(sub):
    runs = L.run_list(); runs = runs[runs.subject == sub]
    assert len(runs) == 2, sub
    R = roi_sets(sub)
    union = np.zeros(next(iter(R.values())).shape, bool)
    for m in R.values():
        union |= m
    fits = [run_fit(sub, r, union) for r in runs.itertuples(index=False)]
    # check against the saved per-run 4 mm effect maps
    chk = {}
    for k, (B, _) in enumerate(fits):
        for j, c in enumerate(("intero", "extero", "match", "C6")):
            saved = nib.load(L.fl_name(sub, f"run{k + 1}", 4, "base", c, "effect")).get_fdata(dtype=np.float32)[union]
            chk[f"run{k + 1}_{c}"] = float(np.abs(saved - B[j]).max())
    out = dict(subject=sub, check_max_abs_diff_vs_saved=chk, rois={})
    for name, m in R.items():
        sel = m[union]
        W, alphas = [], []
        for B, Rw in fits:
            Bw, a = whiten(B[:, sel], Rw[:, sel])
            W.append(np.vstack([Bw, np.zeros((1, sel.sum()))]))         # Rest = zero pattern
            alphas.append(a)
        P = sel.sum()
        d = [float((W[0][i] - W[0][j]) @ (W[1][i] - W[1][j]) / P) for i, j in PAIRS]
        out["rois"][name] = dict(n_vox=int(P), lw_shrinkage=alphas, crossnobis=d)
    os.makedirs(OD, exist_ok=True)
    json.dump(out, open(os.path.join(OD, f"{sub}.json"), "w"), indent=1)
    print(f"{sub}: max |beta - saved| {max(chk.values()):.2e}", flush=True)


def group():
    s, _, _ = L.participant_sets()
    subs = [x for x in s["froi"] if os.path.exists(os.path.join(OD, f"{x}.json"))]
    J = {x: json.load(open(os.path.join(OD, f"{x}.json"))) for x in subs}
    MR = model_rdms()
    mm = pd.DataFrame(MR)
    collin = dict(corr=mm.corr().round(4).to_dict())
    Xm = np.c_[np.ones(len(mm)), mm.values]
    vif = {}
    for j, k in enumerate(mm.columns, start=1):
        y = Xm[:, j]; Z = np.delete(Xm, j, axis=1)
        b, *_ = np.linalg.lstsq(Z, y, rcond=None); r2 = 1 - ((y - Z @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        vif[k] = 1 / max(1 - r2, 1e-12)
    collin["vif"] = vif
    dist, fitr, chk = [], [], []
    for x in subs:
        chk.append(max(J[x]["check_max_abs_diff_vs_saved"].values()))
        for roi, v in J[x]["rois"].items():
            d = np.array(v["crossnobis"])
            for (i, j), dv in zip(PAIRS, d):
                dist.append(dict(subject=x, roi=roi, pair=f"{LAB[i]}-{LAB[j]}", d=dv, n_vox=v["n_vox"]))
            for k, mv in MR.items():
                fitr.append(dict(subject=x, roi=roi, model=k, pearson=np.corrcoef(d, mv)[0, 1],
                                 spearman=stats.spearmanr(d, mv).statistic))
    D = pd.DataFrame(dist); F = pd.DataFrame(fitr)
    rows = []
    for (roi, pair), d in D.groupby(["roi", "pair"], sort=False):
        t = stats.ttest_1samp(d.d, 0)
        rows.append(dict(roi=roi, kind="distance", item=pair, N=len(d), mean=d.d.mean(), sd=d.d.std(ddof=1),
                         dz=d.d.mean() / d.d.std(ddof=1), t=t.statistic, p=t.pvalue, n_vox_median=d.n_vox.median()))
    for (roi, model), d in F.groupby(["roi", "model"], sort=False):
        for meas in ("pearson", "spearman"):
            y = d[meas].dropna(); t = stats.ttest_1samp(y, 0)
            rows.append(dict(roi=roi, kind=f"model_{meas}", item=model, N=len(y), mean=y.mean(), sd=y.std(ddof=1),
                             dz=y.mean() / y.std(ddof=1), t=t.statistic, p=t.pvalue,
                             p_wilcoxon=stats.wilcoxon(y).pvalue))
    T = pd.DataFrame(rows)
    T.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_rsa.csv"), index=False)
    D.round(6).to_csv(os.path.join(OD, "distances_long.csv"), index=False)
    json.dump(dict(model_rdm_collinearity=collin, n_participants=len(subs), max_abs_beta_diff_vs_saved=max(chk),
                   pairs=[f"{LAB[i]}-{LAB[j]}" for i, j in PAIRS], models={k: v.tolist() for k, v in MR.items()}),
              open(os.path.join(RD, "EXPLORATORY_LOC_C3_rsa_models.json"), "w"), indent=1)
    print("N", len(subs), "max beta diff", max(chk)); print(collin)
    print(T[T.kind == "model_pearson"].round(3).to_string())


if __name__ == "__main__":
    if sys.argv[1] == "subject":
        subject(sys.argv[2])
    else:
        group()
