#!/usr/bin/env python
"""
intero2023_loc_depth_c1.py -- EXPLORATORY (logged in Preregistration/localizer_deviations.md,
row "Added exploratory analysis: does the Intero - Extero breathing-depth difference explain C1
across participants?"). Logged and committed before any intercept or slope was computed.

Question: across participants, does the size of the Intero - Extero breathing-depth difference
predict the size of the C1 (Intero - Extero) deactivation, and is C1 still present at zero
depth difference (the regression intercept)?

  python Scripts/intero2023_loc_depth_c1.py check    # label-free: coverage, distributions, reliability
  python Scripts/intero2023_loc_depth_c1.py run      # the logged models
  python Scripts/intero2023_loc_depth_c1.py figure

Depth inputs: derivatives/localizer_qc_cache/physio_stage.pkl (step-1 per-block breathing
measures, BR) and Analysis/localizer_qc/belt_usability.csv; Intero block compliance from
derivatives/localizer_qc_cache/blocks.csv; PetCO2 from Data/regressors/*_petco2.tsv.
Outcomes: the registered L1 per-participant mask means (Analysis/results/LOC_L1_subjects.csv)
and the stage-A FLAME input stacks $OUT/group/LOC_L1 and LOC_C4 (cope.nii.gz, joint 6 mm).
Tables: Analysis/results/EXPLORATORY_LOC_depth_*.csv. Figure: Figures/localizer/loc_depth_c1.png.
No wall-clock time is written to any table.
"""
import os, sys, json, pickle, warnings
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")

REPO = g.REPO
RD = os.path.join(REPO, "Analysis", "results")
FIG = os.path.join(REPO, "Figures", "localizer")
CACHE = os.path.join(g.ROOT, "derivatives", "localizer_qc_cache")
GD = os.path.join(g.OUT, "group")
PFX = "EXPLORATORY_LOC_depth_"
SKIP = 3.0                      # step-1 block window: onset + 3 s to block end (ENTRAIN_SKIP_S)
MIN_PRESSES = L.MIN_PRESSES     # compliant Intero blocks = the blocks modelled as `intero`
MASKS = ["caren_AInsL", "caren_AInsR", "caren_dACC"]
GM = "cortex_gm"                # secondary: Harvard-Oxford cerebral cortex (sub maxprob thr25) in the group mask
HO = os.path.join(os.environ.get("FSLDIR", "/home/normega/fsl"), "data", "atlases", "HarvardOxford",
                  "HarvardOxford-sub-maxprob-thr25-2mm.nii.gz")
HO_CORTEX = [2, 13]
N_BOOT, SEED = 10000, 20261128
LAMBDA_MIN = 0.20               # corrected estimates not identifiable if two-run reliability is below this
LAMBDA_FLOOR = 0.10             # floor for the two-run reliability inside bootstrap draws (count reported)
CONS = {"C1": ("intero", "extero", "intero", "LOC_L1"),     # numerator, denominator, set, group dir
        "C4": ("match", "extero", "all", "LOC_C4")}


# --------------------------------------------------------------------------- depth measures
def block_table():
    """Step-1 per-block breathing measures on usable runs (registered run list x usable belt)."""
    _, _, BR, _, _ = pickle.load(open(os.path.join(CACHE, "physio_stage.pkl"), "rb"))
    U = pd.read_csv(os.path.join(L.QCDIR, "belt_usability.csv"))
    ok = U[U.usable_rate_depth][["subject", "run"]]
    reg = L.run_list()[["subject", "run"]]
    ok = ok.merge(reg, on=["subject", "run"])                  # motion gate etc. (102 registered runs)
    blk = pd.read_csv(os.path.join(CACHE, "blocks.csv"))[["subject", "run", "block_index", "trial_type",
                                                          "onset", "duration", "n_presses"]]
    B = BR.merge(ok, on=["subject", "run"]).merge(blk, on=["subject", "run", "block_index"], how="left")
    assert (B.condition == B.trial_type).all()
    B = B[B.condition.isin(["intero", "extero", "match"])].copy()
    B["compliant"] = (B.condition != "intero") | (B.n_presses.fillna(0) >= MIN_PRESSES)
    B = B[B.compliant]
    # block RVT = median breath amplitude / median breath period, in run units (amp_run cancels
    # in any within-run ratio): depth_rel * rate / 60
    B["rvt_rel"] = B.depth_rel * B.rate_brpm / 60.0
    B = B[np.isfinite(B.rvt_rel) & np.isfinite(B.depth_rel) & (B.depth_rel > 0)]
    return B


def run_ratios(B, num, den):
    """Per run: log(mean block value in `num` / mean in `den`), for RVT and depth."""
    rows = []
    for (s, r), d in B.groupby(["subject", "run"]):
        a, b = d[d.condition == num], d[d.condition == den]
        if len(a) and len(b):
            rows.append(dict(subject=s, run=r, n_num=len(a), n_den=len(b),
                             lr_rvt=np.log(a.rvt_rel.mean() / b.rvt_rel.mean()),
                             lr_depth=np.log(a.depth_rel.mean() / b.depth_rel.mean()),
                             d_rvtz=a.rvt_z_mean.mean() - b.rvt_z_mean.mean()))
    return pd.DataFrame(rows)


def participant_x(R):
    P = R.groupby("subject").agg(k=("run", "size"), lr_rvt=("lr_rvt", "mean"), lr_depth=("lr_depth", "mean"),
                                 d_rvtz=("d_rvtz", "mean"), n_num=("n_num", "sum"), n_den=("n_den", "sum"))
    return P.reset_index()


def reliability(R, col):
    W = R.pivot(index="subject", columns="run", values=col).dropna()
    r = np.corrcoef(W[1], W[2])[0, 1] if len(W) > 2 else np.nan
    return dict(n_two_run=len(W), r_run1_run2=r, spearman_brown=2 * r / (1 + r),
                err_var_run=np.var(W[1] - W[2], ddof=1) / 2, icc31=icc31(W[1].values, W[2].values))


def icc31(a, b):
    X = np.c_[a, b]; n, k = X.shape
    gm = X.mean(); rm = X.mean(1); cm = X.mean(0)
    msr = k * ((rm - gm) ** 2).sum() / (n - 1)
    mse = ((X - rm[:, None] - cm[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse)


def petco2_runs(B, num, den):
    """Per run: mean petco2_raw (lag 0, demeaned within run, %CO2) over each condition's block
    windows (onset + 3 s to block end), num - den. Runs with has_petco2 in the manifest."""
    man = L.run_list().set_index(["subject", "run"])
    rows = []
    for (s, r), d in B.groupby(["subject", "run"]):
        if not bool(man.loc[(s, r), "has_petco2"]):
            continue
        f = os.path.join(L.REGDIR, f"{s}_{L.entity(r)}_petco2.tsv")
        if not os.path.exists(f):
            continue
        p = pd.read_csv(f, sep="\t", na_values="n/a").petco2_raw.values
        vt = g.TR * np.arange(len(p))
        val = {}
        for c in (num, den):
            m = np.zeros(len(p), bool)
            for b in d[d.condition == c].itertuples():
                m |= (vt >= b.onset + SKIP) & (vt < b.onset + b.duration)
            v = p[m]
            val[c] = np.nanmean(v) if np.isfinite(v).sum() >= 3 else np.nan
        if np.isfinite(val[num]) and np.isfinite(val[den]):
            rows.append(dict(subject=s, run=r, d_petco2=val[num] - val[den]))
    return pd.DataFrame(rows)


def check():
    """Label-free: coverage, distributions and run-1 vs run-2 reliability of the depth measures."""
    s, _, _ = L.participant_sets()
    B = block_table()
    out = []
    for con, (num, den, st, _) in CONS.items():
        R = run_ratios(B, num, den)
        R = R[R.subject.isin(s[st])]
        P = participant_x(R)
        rec = dict(contrast=con, set=st, n_set=len(s[st]), n_with_physio=len(P),
                   n_two_run=int((P.k == 2).sum()), n_runs=len(R))
        for c in ("lr_rvt", "lr_depth", "d_rvtz"):
            x = P[c]
            rec.update({f"{c}_mean": x.mean(), f"{c}_sd": x.std(ddof=1), f"{c}_median": x.median(),
                        f"{c}_min": x.min(), f"{c}_max": x.max(), f"{c}_n_pos": int((x > 0).sum()),
                        f"{c}_skew": stats.skew(x)})
            rr = reliability(R, c)
            rec.update({f"{c}_{k}": v for k, v in rr.items()})
            W = R.pivot(index="subject", columns="run", values=c).dropna()
            ev, vt = calib_params(W[1].values, W[2].values)
            rec.update({f"{c}_true_var": vt, f"{c}_rel_two_run_sb": vt / (vt + ev / 2),
                        f"{c}_rel_one_run": vt / (vt + ev),
                        f"{c}_var_one_run_participants": P[P.k == 1][c].var(ddof=1),
                        f"{c}_var_two_run_runlevel": np.var(np.r_[W[1], W[2]], ddof=1)})
        rec["r_lr_rvt_lr_depth"] = np.corrcoef(P.lr_rvt, P.lr_depth)[0, 1]
        rec["r_lr_rvt_d_rvtz"] = np.corrcoef(P.lr_rvt, P.d_rvtz)[0, 1]
        Q = petco2_runs(B, num, den)
        Q = Q[Q.subject.isin(P.subject)]
        rec["n_petco2"] = Q.subject.nunique()
        out.append(rec)
    C = pd.DataFrame(out)
    C.to_csv(os.path.join(RD, PFX + "checks.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 200)
    print(C.T.to_string())


# --------------------------------------------------------------------------- outcomes
def outcomes(con):
    """Per-participant mask means of the joint 6 mm effect (the FLAME input stack)."""
    import nibabel as nib
    from nilearn import image
    num, den, st, gd = CONS[con]
    d = os.path.join(GD, gd)
    subs = json.load(open(os.path.join(d, "loc_meta.json")))["subjects"]
    ref = nib.load(g.GROUP_MASK); bm = ref.get_fdata() > 0
    cope = nib.load(os.path.join(d, "cope.nii.gz")).get_fdata(dtype=np.float32)
    assert cope.shape[-1] == len(subs)
    masks = {}
    for r in MASKS:
        masks[r] = bm & (image.resample_to_img(os.path.join(REPO, "ROIs", f"{r}.nii"), ref, interpolation="nearest",
                                                force_resample=True, copy_header=True).get_fdata() > 0)
    ho = image.resample_to_img(HO, ref, interpolation="nearest", force_resample=True, copy_header=True).get_fdata()
    masks[GM] = bm & np.isin(ho, HO_CORTEX)
    Y = pd.DataFrame({r: cope[m].mean(0) for r, m in masks.items()}, index=subs)
    Y.index.name = "subject"
    if con == "C1":   # must reproduce the registered extraction
        ref_ = pd.read_csv(os.path.join(RD, "LOC_L1_subjects.csv")).pivot(index="subject", columns="roi",
                                                                           values="mask_mean")
        for r in MASKS:
            assert np.allclose(Y.loc[ref_.index, r].values, ref_[r].values, atol=1e-5), r
    return Y.reset_index(), {r: int(m.sum()) for r, m in masks.items()}


# --------------------------------------------------------------------------- models
def ols_hc3(x, y):
    X = np.c_[np.ones(len(x)), x]
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    V = XtXi @ (X.T * (e / (1 - h)) ** 2) @ X @ XtXi
    se = np.sqrt(np.diag(V))
    return b, se, h


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    adj = np.empty(m); run = 0
    for i, j in enumerate(o):
        run = max(run, (m - i) * p[j]); adj[j] = min(1, run)
    return adj


def calib_params(x1, x2):
    """From two-run participants: run-level error variance and true-score variance (cov of runs)."""
    ev = np.var(x1 - x2, ddof=1) / 2
    vt = np.cov(x1, x2)[0, 1]
    return ev, vt


def calibrate(x, k, ev, vt):
    """Regression calibration: E[T | x] = mu + lambda_i (x - mu), lambda_i = vt / (vt + ev / k_i)."""
    li = vt / (vt + ev / k)
    return x.mean() + li * (x - x.mean()), li


def fit_one(x, y, k, R2, rng):
    n = len(x); df = n - 2
    b, se, h = ols_hc3(x, y)
    t = b / se; p = 2 * stats.t.sf(np.abs(t), df); tc = stats.t.ppf(.975, df)
    W = R2.pivot(index="subject", columns="run", values="x")
    subs = R2.subject.drop_duplicates().values
    assert len(subs) == n
    two = np.array([s_ in W.dropna().index for s_ in subs])
    x1 = np.array([W.loc[s_, 1] if t_ else np.nan for s_, t_ in zip(subs, two)])
    x2 = np.array([W.loc[s_, 2] if t_ else np.nan for s_, t_ in zip(subs, two)])
    ev, vt = calib_params(x1[two], x2[two])
    sb = vt / (vt + ev / 2)                       # Spearman-Brown reliability of a two-run mean
    ident = bool(sb >= LAMBDA_MIN)
    xh, li = calibrate(x, k, ev, vt)
    bc, _, _ = ols_hc3(xh, y)
    # participant bootstrap of the corrected estimates (ev, vt re-estimated in every draw)
    boot = np.full((N_BOOT, 2), np.nan); nfloor = 0
    for i in range(N_BOOT):
        j = rng.integers(0, n, n)
        tj = two[j]
        if tj.sum() < 3:
            continue
        evb, vtb = calib_params(x1[j][tj], x2[j][tj])
        vmin = LAMBDA_FLOOR / (1 - LAMBDA_FLOOR) * evb / 2      # two-run reliability floored at LAMBDA_FLOOR
        if vtb < vmin:
            vtb = vmin; nfloor += 1
        xhb, _ = calibrate(x[j], k[j], evb, vtb)
        boot[i] = ols_hc3(xhb, y[j])[0]
    # leave-one-out range of the intercept (influence)
    loo = np.array([ols_hc3(np.delete(x, i), np.delete(y, i))[0][0] for i in range(n)])
    return dict(N=n, n_two_run=int(two.sum()), mean_x=x.mean(), sd_x=x.std(ddof=1), mean_y=y.mean(),
                dz_y=y.mean() / y.std(ddof=1),
                intercept=b[0], intercept_se_hc3=se[0], intercept_lo=b[0] - tc * se[0],
                intercept_hi=b[0] + tc * se[0], intercept_t=t[0], intercept_p=p[0],
                intercept_over_mean=b[0] / y.mean(),
                slope=b[1], slope_se_hc3=se[1], slope_lo=b[1] - tc * se[1], slope_hi=b[1] + tc * se[1],
                slope_t=t[1], slope_p=p[1], r_xy=np.corrcoef(x, y)[0, 1],
                rho_xy=stats.spearmanr(x, y)[0], max_leverage=h.max(),
                intercept_loo_min=loo.min(), intercept_loo_max=loo.max(),
                err_var_run=ev, true_var=vt, r_run1_run2=np.corrcoef(x1[two], x2[two])[0, 1],
                rel_two_run_sb=sb, rel_one_run=vt / (vt + ev), lambda_mean=li.mean(),
                corr_identifiable=ident,
                intercept_corr=bc[0], intercept_corr_lo=np.nanpercentile(boot[:, 0], 2.5),
                intercept_corr_hi=np.nanpercentile(boot[:, 0], 97.5),
                slope_corr=bc[1], slope_corr_lo=np.nanpercentile(boot[:, 1], 2.5),
                slope_corr_hi=np.nanpercentile(boot[:, 1], 97.5),
                boot_n_floored=nfloor, boot_n_valid=int(np.isfinite(boot[:, 0]).sum()))


def run():
    s, _, _ = L.participant_sets()
    B = block_table()
    rows, subj = [], []
    for con, (num, den, st, _) in CONS.items():
        R = run_ratios(B, num, den); R = R[R.subject.isin(s[st])]
        P = participant_x(R)
        Y, nvox = outcomes(con)
        D = P.merge(Y, on="subject")
        Q = petco2_runs(B, num, den); Q = Q[Q.subject.isin(D.subject)]
        PQ = Q.groupby("subject").d_petco2.mean().rename("d_petco2").reset_index()
        D = D.merge(PQ, on="subject", how="left")
        D.insert(0, "contrast", con)
        subj.append(D)
        for xm, fam in (("lr_rvt", "primary" if con == "C1" else "secondary_C4"),
                        ("lr_depth", "secondary_depth" if con == "C1" else "secondary_C4_depth")):
            R2 = R[R.subject.isin(D.subject)][["subject", "run", xm]].rename(columns={xm: "x"})
            for r in MASKS + [GM]:
                o = fit_one(D[xm].values, D[r].values, D.k.values.astype(float), R2, np.random.default_rng(SEED))
                rows.append(dict(contrast=con, x=xm, family=fam, region=r, n_vox=nvox[r], **o))
        # sensitivity: two-run participants only (every x is a two-run mean), primary measure
        D2 = D[D.k == 2]
        R2 = R[R.subject.isin(D2.subject)][["subject", "run", "lr_rvt"]].rename(columns={"lr_rvt": "x"})
        for r in MASKS + [GM]:
            o = fit_one(D2.lr_rvt.values, D2[r].values, D2.k.values.astype(float), R2, np.random.default_rng(SEED))
            rows.append(dict(contrast=con, x="lr_rvt", family="sens_two_run_only", region=r, n_vox=nvox[r], **o))
        # PetCO2 (descriptive): C ~ dPetCO2, and dPetCO2 vs log RVT ratio
        E = D[D.d_petco2.notna()]
        for r in MASKS + [GM]:
            x, y = E.d_petco2.values, E[r].values
            b, se, _ = ols_hc3(x, y); df = len(x) - 2; tc = stats.t.ppf(.975, df)
            rows.append(dict(contrast=con, x="d_petco2", family="descriptive_petco2", region=r, n_vox=nvox[r],
                             N=len(x), mean_x=x.mean(), sd_x=x.std(ddof=1), mean_y=y.mean(),
                             dz_y=y.mean() / y.std(ddof=1),
                             intercept=b[0], intercept_se_hc3=se[0], intercept_lo=b[0] - tc * se[0],
                             intercept_hi=b[0] + tc * se[0], intercept_t=b[0] / se[0],
                             intercept_p=2 * stats.t.sf(abs(b[0] / se[0]), df),
                             slope=b[1], slope_se_hc3=se[1], slope_lo=b[1] - tc * se[1], slope_hi=b[1] + tc * se[1],
                             slope_t=b[1] / se[1], slope_p=2 * stats.t.sf(abs(b[1] / se[1]), df),
                             r_xy=np.corrcoef(x, y)[0, 1], rho_xy=stats.spearmanr(x, y)[0],
                             petco2_dz=x.mean() / x.std(ddof=1),
                             petco2_p_onesample=stats.ttest_1samp(x, 0).pvalue,
                             rho_petco2_lr_rvt=stats.spearmanr(x, E.lr_rvt.values)[0]))
    T = pd.DataFrame(rows)
    fam = T.family.isin(["primary", "secondary_depth", "secondary_C4", "secondary_C4_depth"])
    for (con, xm), d in T[T.region.isin(MASKS) & fam].groupby(["contrast", "x"]):
        T.loc[d.index, "intercept_p_holm3"] = holm(d.intercept_p.values)
        T.loc[d.index, "slope_p_holm3"] = holm(d.slope_p.values)
    T.to_csv(os.path.join(RD, PFX + "models.csv"), index=False)
    S = pd.concat(subj, ignore_index=True)
    S.to_csv(os.path.join(RD, PFX + "subjects.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
    cols = ["contrast", "x", "region", "N", "mean_y", "intercept", "intercept_lo", "intercept_hi", "intercept_p",
            "intercept_p_holm3", "slope", "slope_lo", "slope_hi", "slope_p", "slope_p_holm3", "rel_two_run_sb",
            "intercept_corr", "intercept_corr_lo", "intercept_corr_hi", "slope_corr", "slope_corr_lo",
            "slope_corr_hi", "intercept_loo_min", "intercept_loo_max", "max_leverage", "boot_n_floored"]
    print(T[T.x != "d_petco2"][["family"] + cols].round(4).to_string(index=False))
    print(T[T.x == "d_petco2"][["contrast", "region", "N", "mean_x", "petco2_dz", "petco2_p_onesample",
                                "rho_petco2_lr_rvt", "intercept", "intercept_p", "slope", "slope_p"]]
          .round(4).to_string(index=False))


# --------------------------------------------------------------------------- figure
def figure():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T = pd.read_csv(os.path.join(RD, PFX + "models.csv"))
    S = pd.read_csv(os.path.join(RD, PFX + "subjects.csv"))
    INK, MUTED, GRID, ACC, ACC2 = "#0b0b0b", "#52514e", "#dcdbd6", "#2a78d6", "#eb6834"
    names = {"caren_AInsL": "AInsL", "caren_AInsR": "AInsR", "caren_dACC": "dACC", GM: "Cortical grey matter"}
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.8), sharex="row")
    for i, con in enumerate(["C1", "C4"]):
        D = S[S.contrast == con]
        for j, r in enumerate(MASKS + [GM]):
            ax = axs[i, j]
            t = T[(T.contrast == con) & (T.x == "lr_rvt") & (T.region == r) & ~T.family.str.startswith("sens")].iloc[0]
            ax.axhline(0, color=GRID, lw=1, zorder=0); ax.axvline(0, color=GRID, lw=1, zorder=0)
            ax.scatter(D.lr_rvt, D[r], s=18, color=ACC, alpha=.75, lw=0)
            xx = np.linspace(min(D.lr_rvt.min(), 0) - .05, D.lr_rvt.max() + .05, 50)
            ax.plot(xx, t.intercept + t.slope * xx, color=ACC, lw=1.6, label="OLS")
            if t.corr_identifiable:
                ax.plot(xx, t.intercept_corr + t.slope_corr * xx, color=ACC2, lw=1.4, ls="--",
                        label="corrected for unreliability (very imprecise)")
            ax.errorbar([0], [t.intercept], yerr=[[t.intercept - t.intercept_lo], [t.intercept_hi - t.intercept]],
                        fmt="o", color=INK, ms=5, capsize=3, zorder=5)
            p = t.intercept_p_holm3 if r in MASKS else t.intercept_p
            ax.set_title(f"{con} {names[r]}\nintercept {t.intercept:+.3f} [{t.intercept_lo:+.3f}, "
                         f"{t.intercept_hi:+.3f}], p{'(Holm)' if r in MASKS else ''} {p:.2g}",
                         fontsize=8.5, color=INK, loc="left")
            ax.tick_params(labelsize=8, colors=MUTED)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if j == 0:
                ax.set_ylabel(f"{con} mask mean (PSC)", fontsize=9, color=INK)
            ax.set_xlabel("log RVT ratio (" + ("Match" if con == "C4" else "Intero") + " / Extero)",
                          fontsize=9, color=INK)
    axs[0, 0].legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.suptitle("Condition difference against the breathing-depth difference across participants "
                 "(black: intercept at equal depth, 95% HC3 CI)", fontsize=10, color=INK, x=.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, .95))
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "loc_depth_c1.png"), dpi=150)


if __name__ == "__main__":
    {"check": check, "run": run, "figure": figure}[sys.argv[1]]()
