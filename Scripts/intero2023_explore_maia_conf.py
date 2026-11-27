#!/usr/bin/env python
"""
intero2023_explore_maia_conf.py -- EXPLORATORY E9 (deviations.md, "After the update", 2026-11-26):
does the pacing-related ACC response that tracks MAIA (H6) also relate to confidence?

  python Scripts/intero2023_explore_maia_conf.py behav   # behavioural checks only, no imaging values
  python Scripts/intero2023_explore_maia_conf.py run     # the logged tests (imaging)
  python Scripts/intero2023_explore_maia_conf.py plot    # scatter PNG -> derivatives (not the repo)

Inputs (no new first-level fitting): <sub>_var-base_con-C1_effect (paced - free, all 47),
<sub>_var-base_con-C3_effect (Model B confidence modulator, N 44), <sub>_var-base_con-E8confjoint_effect
(E8 rating-period confidence modulator, N 44); Data/maia_scores.csv (maia_ieat_total, 46 of 47);
confidence of every paced epoch in the 90 registered runs (desc-model_events, the trials modelled).

conf_i   = mean confidence (1-7) over participant i's paced epochs (aware, unaware, wrong-direction,
           No Change; 7-14 epochs), all trials, no accuracy split.
C1_i     = mean of participant i's C1 effect over the registered ieat_ACCsphere (528 voxels).

PRIMARY  E9a     Spearman rho(C1_i, conf_i), N 47. Permutation of conf across participants,
                 10,000, seed 20261126, two-sided, p = (1 + #|rho*| >= |rho|) / 10,001.
                 Alongside: Pearson r and the JZS correlation BF10 (Wetzels & Wagenmakers 2012),
                 prior scale s = 1 (their default), s = 0.707 and 0.354 as robustness.
SECONDARY (Holm over the secondary p-values)
  E9b     rank partial correlation of C1_i and conf_i given MAIA (Pearson on residualised ranks), N 46;
          permutation of the conf residuals (Freedman-Lane), 10,000, two-sided. Descriptive chain:
          rho(MAIA, C1), rho(MAIA, conf), rho(C1, conf) in the same 46.
  E9c     MAIA (centred) as a covariate on C3 and on E8confjoint, in ieat_ACCsphere, voxel max-t FWE by
          permutation of the covariate (nilearn non_parametric_inference, as H6), two-sided, 10,000,
          seed 20261126; Bonferroni x2 over the two maps. N 43 each.
  E9a-LOO C1 in a leave-one-out H6 cluster: for each participant i, OLS t map of centred MAIA on C1
          over the group mask from the H6 sample minus i; threshold t > t_.001(df) one-tailed positive;
          6-connected components; take the one containing (3.5, 35.5, 3.5), else the one with the
          voxel nearest that point if within 10 mm, else the fold is EMPTY. sub-065 (no MAIA) uses
          the full H6 sample (it is never in it). C1_i(LOO) = mean of i's C1 in its fold cluster;
          then as E9a (Spearman, permutation, N 47). If any fold is empty, the LOO secondary is
          dropped and only the fold count is reported.
DESCRIPTIVE (no inference)
  type-2 AUC (Data/behaviour_subject.csv, N 31) vs C1_i: Spearman rho with a bootstrap 95% CI
  (5,000 resamples, seed 20261126). C1 sphere ICC(3,1) run 1 vs run 2 (two-run participants, N 43).
  Confidence split-half reliability (odd vs even paced epochs, Spearman-Brown), N 47.
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import integrate, ndimage, special, stats
import intero2023_glm_common as g
q = g.q
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
OD = os.path.join(g.OUT, "group", "EXPLORATORY_E9")
ROI_NAME = "ieat_ACCsphere"
PEAK = (3.5, 35.5, 3.5)          # H6 whole-brain positive cluster peak (pregenual ACC)
MAX_DIST_MM = 10.0
N_PERM, SEED = 10000, 20261126


# ------------------------------------------------------------------ behaviour (no imaging)
def behaviour():
    man = g.run_list()
    ev = pd.concat([pd.read_csv(f"{q.FL}/{s}_{e}_desc-model_events.tsv", sep="\t").assign(subject=s, ent=e)
                    for s, e in zip(man.subject, man.bids_entity)], ignore_index=True)
    p = ev[ev.condition.str.startswith("paced")].sort_values(["subject", "ent", "onset"])
    assert p.confidence.notna().all()
    d = p.groupby("subject").confidence.agg(n_paced="count", conf="mean", conf_sd="std")
    halves = p.assign(k=p.groupby("subject").cumcount() % 2).groupby(["subject", "k"]).confidence.mean().unstack()
    d["conf_odd"], d["conf_even"] = halves[0], halves[1]
    maia = pd.read_csv(os.path.join(g.REPO, "Data", "maia_scores.csv")).set_index("subject")
    d["maia"] = maia.maia_ieat_total.reindex(d.index)
    b = pd.read_csv(os.path.join(g.REPO, "Data", "behaviour_subject.csv")).set_index("subject")
    d["type2_auc"] = b.type2_auc.reindex(d.index)
    d["p_correct_dir"] = b.p_correct_dir.reindex(d.index)
    d["frac_conf_ge6"] = p.assign(h=p.confidence >= 6).groupby("subject").h.mean()
    return d, p


def behav():
    d, p = behaviour()
    print(f"N = {len(d)}; paced epochs {len(p)}; per person {d.n_paced.min()}-{d.n_paced.max()} "
          f"(median {d.n_paced.median():.0f})")
    print(f"mean confidence: range {d.conf.min():.2f}-{d.conf.max():.2f}, SD across people {d.conf.std():.2f}; "
          f"ratings >= 6: {100 * (p.confidence >= 6).mean():.0f}%")
    print("constant confidence:", d[d.conf_sd == 0].index.tolist())
    print("no MAIA:", d[d.maia.isna()].index.tolist())
    print(f"type-2 AUC available: {d.type2_auc.notna().sum()}")
    x = d.dropna(subset=["maia"])
    print(f"MAIA-conf  Spearman {stats.spearmanr(x.maia, x.conf).statistic:.3f}, "
          f"Pearson {stats.pearsonr(x.maia, x.conf).statistic:.3f}, n {len(x)}")
    print(f"MAIA-acc   Spearman {stats.spearmanr(x.maia, x.p_correct_dir).statistic:.3f}")
    for c in ("C1", "C3", "E8confjoint"):
        s = [x_ for x_ in d.index if os.path.exists(fl(x_, c))]
        print(f"{c}: {len(s)} maps, {sum(pd.notna(d.loc[s, 'maia']))} with MAIA")


# ------------------------------------------------------------------ statistics helpers
def fl(sub, con, run=None):
    pre = f"{sub}_{run}_" if run else f"{sub}_"
    return os.path.join(g.OUT, "firstlevel", sub, f"{pre}var-base_con-{con}_effect.nii.gz")


def rank(x):
    return stats.rankdata(x)


def perm_corr(x, y, rng):
    """Two-sided permutation p for Pearson r of x, y (pass ranks for Spearman)."""
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
    r = float((x * y).mean())
    P = np.array([rng.permutation(y) for _ in range(N_PERM)])
    rp = (P * x).mean(1)
    return r, float((1 + (np.abs(rp) >= abs(r) - 1e-12).sum()) / (N_PERM + 1))


def resid(y, z):
    Z = np.column_stack([np.ones_like(z), z])
    return y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]


def jzs_corr_bf10(r, n, s=1.0):
    """JZS Bayes factor for a Pearson correlation (Wetzels & Wagenmakers 2012, eq. 13), with the
    Zellner-Siow prior g ~ InvGamma(1/2, s^2 n / 2); s = 1 is their default."""
    a = s ** 2 * n / 2
    f = lambda gg: np.exp(0.5 * np.log(a) - special.gammaln(0.5) - 1.5 * np.log(gg) - a / gg
                          + (n - 2) / 2 * np.log1p(gg) - (n - 1) / 2 * np.log1p((1 - r ** 2) * gg))
    return float(integrate.quad(f, 0, np.inf, limit=200)[0])


def icc31(a, b):
    X = np.column_stack([a, b]); n, k = X.shape
    ms_r = k * ((X.mean(1) - X.mean()) ** 2).sum() / (n - 1)
    ms_e = (((X - X.mean(1, keepdims=True) - X.mean(0, keepdims=True) + X.mean()) ** 2).sum()) / ((n - 1) * (k - 1))
    return (ms_r - ms_e) / (ms_r + (k - 1) * ms_e)


def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); out = np.empty(m); run = 0.0
    for i, j in enumerate(o):
        run = max(run, min(1.0, (m - i) * ps[j])); out[j] = run
    return out


# ------------------------------------------------------------------ imaging
def load_masked(subs, con, mask):
    Y = np.vstack([nib.load(fl(s, con)).get_fdata(dtype=np.float32)[mask] for s in subs]).astype(np.float64)
    assert np.isfinite(Y).all(), f"non-finite values in {con}"
    return Y


def ols_t_vec(Y, x):
    X = np.column_stack([np.ones_like(x), x - x.mean()])
    XtXi = np.linalg.inv(X.T @ X); B = XtXi @ X.T @ Y; R = Y - X @ B
    s2 = (R ** 2).sum(0) / (len(x) - 2)
    return B[1] / np.sqrt(s2 * XtXi[1, 1])


def fold_cluster(tv, mask, df, affine, peak_ijk, peak_xyz):
    """Component of the one-tailed positive p < .001 map holding (or nearest within 10 mm of) the peak."""
    sup = np.zeros(mask.shape, bool); sup[mask] = tv > stats.t.isf(0.001, df)
    lab, n = ndimage.label(sup, structure=ndimage.generate_binary_structure(3, 1))
    if n == 0:
        return None, np.nan, False
    k = lab[peak_ijk]
    if k:
        return lab == k, 0.0, True
    ijk = np.argwhere(sup)
    dist = np.linalg.norm(nib.affines.apply_affine(affine, ijk) - np.array(peak_xyz), axis=1)
    j = int(np.argmin(dist))
    if dist[j] > MAX_DIST_MM:
        return None, float(dist[j]), False
    return lab == lab[tuple(ijk[j])], float(dist[j]), False


def run():
    import intero2023_group as G
    from nilearn.glm.second_level import non_parametric_inference
    os.makedirs(OD, exist_ok=True)
    d, _ = behaviour()
    gm_img = nib.load(g.GROUP_MASK); gmask = gm_img.get_fdata() > 0
    roi = nib.load(G.ROI(ROI_NAME)).get_fdata() > 0
    assert np.allclose(gm_img.affine, nib.load(G.ROI(ROI_NAME)).affine)
    subs = [s for s in d.index if os.path.exists(fl(s, "C1"))]
    assert len(subs) == 47
    d = d.loc[subs]
    d["C1_sphere"] = [nib.load(fl(s, "C1")).get_fdata(dtype=np.float32)[roi].mean() for s in subs]
    rows = []

    # ---------------- E9a primary
    rng = np.random.default_rng(SEED)
    rho, p = perm_corr(rank(d.C1_sphere.values), rank(d.conf.values), rng)
    r = stats.pearsonr(d.C1_sphere, d.conf).statistic
    row = dict(analysis="E9a", role="primary", x="C1_sphere", y="conf", covariate="", N=len(d),
               spearman_rho=rho, p_perm=p, pearson_r=r)
    for s_ in (1.0, 0.707, 0.354):
        row[f"BF10_s{s_}"] = jzs_corr_bf10(r, len(d), s_)
    rows.append(row)

    # ---------------- E9b partial (descriptive chain alongside)
    m = d.dropna(subset=["maia"])
    rng = np.random.default_rng(SEED)
    ex = resid(rank(m.C1_sphere.values), rank(m.maia.values)); ey = resid(rank(m.conf.values), rank(m.maia.values))
    pr, pp = perm_corr(ex, ey, rng)
    rows.append(dict(analysis="E9b", role="secondary", x="C1_sphere", y="conf", covariate="maia", N=len(m),
                     spearman_rho=pr, p_perm=pp))
    for a, b in (("maia", "C1_sphere"), ("maia", "conf"), ("C1_sphere", "conf")):
        rows.append(dict(analysis="E9b chain", role="descriptive", x=a, y=b, covariate="", N=len(m),
                         spearman_rho=stats.spearmanr(m[a], m[b]).statistic,
                         pearson_r=stats.pearsonr(m[a], m[b]).statistic))

    # ---------------- E9a-LOO
    h6 = [s for s in subs if pd.notna(d.loc[s, "maia"])]
    Y = load_masked(subs, "C1", gmask); idx = {s: i for i, s in enumerate(subs)}
    peak_ijk = tuple(np.round(nib.affines.apply_affine(np.linalg.inv(gm_img.affine), PEAK)).astype(int))
    Yh = Y[[idx[s] for s in h6]]; xh = d.loc[h6, "maia"].values
    full, _, _ = fold_cluster(ols_t_vec(Yh, xh), gmask, len(h6) - 2, gm_img.affine, peak_ijk, PEAK)
    full_size = int(full.sum()) if full is not None else 0
    print(f"full-sample check: cluster at the H6 peak = {full_size} voxels (H6 reported 687)")
    folds = []
    for s in subs:
        dset = [x_ for x_ in h6 if x_ != s]
        sel = [idx[x_] for x_ in dset]
        cl, dist, has_peak = fold_cluster(ols_t_vec(Y[sel], d.loc[dset, "maia"].values), gmask,
                                          len(dset) - 2, gm_img.affine, peak_ijk, PEAK)
        if cl is None:
            folds.append(dict(subject=s, n_def=len(dset), empty=True, size_vox=0, dist_mm=dist,
                              contains_peak=False, dice_full=np.nan, C1_loo=np.nan)); continue
        dice = 2 * (cl & full).sum() / (cl.sum() + full.sum()) if full is not None else np.nan
        val = nib.load(fl(s, "C1")).get_fdata(dtype=np.float32)[cl].mean()
        folds.append(dict(subject=s, n_def=len(dset), empty=False, size_vox=int(cl.sum()), dist_mm=dist,
                          contains_peak=has_peak, dice_full=float(dice), C1_loo=float(val)))
    F = pd.DataFrame(folds); F.to_csv(os.path.join(RD, "EXPLORATORY_E9_loo_folds.csv"),
                                       index=False, columns=[c for c in F.columns if c != "C1_loo"])
    n_empty = int(F["empty"].sum())
    if n_empty == 0:
        d["C1_loo"] = F.set_index("subject").C1_loo
        rng = np.random.default_rng(SEED)
        rl, pl = perm_corr(rank(d.C1_loo.values), rank(d.conf.values), rng)
        rows.append(dict(analysis="E9a-LOO", role="secondary", x="C1_loo", y="conf", covariate="", N=len(d),
                         spearman_rho=rl, p_perm=pl, pearson_r=stats.pearsonr(d.C1_loo, d.conf).statistic,
                         note=f"fold clusters {F.size_vox.min()}-{F.size_vox.max()} vox, "
                              f"median {F.size_vox.median():.0f}; Dice vs full {F.dice_full.min():.2f}-"
                              f"{F.dice_full.max():.2f}"))
        rows.append(dict(analysis="E9a-LOO chain", role="descriptive", x="maia", y="C1_loo", covariate="",
                         N=len(h6), spearman_rho=stats.spearmanr(d.loc[h6, "maia"], d.loc[h6, "C1_loo"]).statistic))
    else:
        rows.append(dict(analysis="E9a-LOO", role="secondary (dropped)", N=len(d),
                         note=f"{n_empty} of {len(F)} folds empty; dropped as logged"))

    # ---------------- E9c MAIA x confidence-modulator maps, sphere
    crow = []
    for con in ("C3", "E8confjoint"):
        cs = [s for s in subs if os.path.exists(fl(s, con)) and pd.notna(d.loc[s, "maia"])]
        imgs = [nib.load(fl(s, con)) for s in cs]
        x = d.loc[cs, "maia"].values
        dm = pd.DataFrame({"intercept": 1.0, "maia": x - x.mean()})
        lp = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="maia", mask=G.ROI(ROI_NAME),
                                      n_perm=N_PERM, two_sided_test=True, random_state=SEED, n_jobs=G.N_JOBS)
        tmap = G.ols_t(imgs, dm, "maia", G.ROI(ROI_NAME))
        lpv = lp.get_fdata()[roi]; tv = tmap.get_fdata()[roi]
        k = int(np.nanargmax(np.abs(tv))); ijk = np.argwhere(roi)[k]
        pmin = float(10 ** -np.nanmax(lpv))
        means = np.array([i.get_fdata(dtype=np.float32)[roi].mean() for i in imgs])
        crow.append(dict(analysis="E9c", role="secondary", map=con, roi=ROI_NAME, N=len(cs),
                         peak_abs_t_signed=round(float(tv[k]), 3),
                         peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(tmap.affine, ijk)],
                         p_fwe=pmin, p_fwe_bonf2=min(1.0, 2 * pmin), n_vox_p_lt_05=int((lpv > -np.log10(0.05)).sum()),
                         roi_mean_rho_maia=stats.spearmanr(x, means).statistic,
                         roi_mean=float(means.mean()), roi_mean_dz=float(means.mean() / means.std(ddof=1))))
        g.save_float32(tmap, os.path.join(OD, f"E9c_{con}_maia_t.nii.gz"))
        g.save_float32(lp, os.path.join(OD, f"E9c_{con}_maia_logpFWE.nii.gz"))
    C = pd.DataFrame(crow)

    # ---------------- Holm over the secondary family
    R = pd.DataFrame(rows)
    sec = [("E9b", R.loc[R.analysis == "E9b", "p_perm"].iloc[0])]
    sec += [(f"E9c {c}", p_) for c, p_ in zip(C["map"], C.p_fwe)]
    if n_empty == 0:
        sec.append(("E9a-LOO", R.loc[R.analysis == "E9a-LOO", "p_perm"].iloc[0]))
    hp = dict(zip([k for k, _ in sec], holm([p_ for _, p_ in sec])))
    R["p_holm_secondary"] = [hp.get(a, np.nan) for a in R.analysis]
    C["p_holm_secondary"] = [hp[f"E9c {c}"] for c in C["map"]]
    C["holm_m"] = R["holm_m"] = len(sec)

    # ---------------- descriptives
    a = d.dropna(subset=["type2_auc"])
    rng = np.random.default_rng(SEED)
    bs = []
    for _ in range(5000):
        i = rng.integers(0, len(a), len(a))
        bs.append(stats.spearmanr(a.C1_sphere.values[i], a.type2_auc.values[i]).statistic)
    R = pd.concat([R, pd.DataFrame([dict(analysis="type-2 AUC", role="descriptive", x="C1_sphere", y="type2_auc",
                                         N=len(a), spearman_rho=stats.spearmanr(a.C1_sphere, a.type2_auc).statistic,
                                         note=f"bootstrap 95% CI [{np.nanpercentile(bs, 2.5):.3f}, "
                                              f"{np.nanpercentile(bs, 97.5):.3f}]")])], ignore_index=True)
    man = g.run_list(); two = sorted(man.groupby("subject").size().pipe(lambda x: x[x == 2]).index)
    r1 = [nib.load(fl(s, "C1", "run-1")).get_fdata(dtype=np.float32)[roi].mean() for s in two]
    r2 = [nib.load(fl(s, "C1", "run-2")).get_fdata(dtype=np.float32)[roi].mean() for s in two]
    sh = stats.pearsonr(d.conf_odd, d.conf_even).statistic
    rel = pd.DataFrame([dict(measure="C1 ieat_ACCsphere mean, run 1 vs run 2", N=len(two),
                             icc31=icc31(np.array(r1), np.array(r2)), pearson_r=stats.pearsonr(r1, r2).statistic),
                        dict(measure="mean paced confidence, odd vs even epochs (Spearman-Brown)", N=len(d),
                             pearson_r=sh, spearman_brown=2 * sh / (1 + sh)),
                        dict(measure="mean paced confidence across participants", N=len(d),
                             mean=d.conf.mean(), sd=d.conf.std(), min=d.conf.min(), max=d.conf.max(),
                             frac_ratings_ge6=float((d.frac_conf_ge6 * d.n_paced).sum() / d.n_paced.sum()))])
    R.to_csv(os.path.join(RD, "EXPLORATORY_E9_correlations.csv"), index=False)
    C.to_csv(os.path.join(RD, "EXPLORATORY_E9c_maia_confmod.csv"), index=False)
    rel.to_csv(os.path.join(RD, "EXPLORATORY_E9_reliability.csv"), index=False)
    json.dump(dict(E9a=subs, E9b=h6, E9c={c: [s for s in subs if os.path.exists(fl(s, c)) and s in h6]
                                         for c in ("C3", "E8confjoint")}, type2_auc=a.index.tolist()),
              open(os.path.join(RD, "EXPLORATORY_E9_subjects.json"), "w"), indent=1)
    d.to_csv(os.path.join(OD, "E9_subject_values.csv"))          # per-participant values stay in derivatives
    pd.set_option("display.width", 250)
    print(R.round(4).to_string(index=False)); print(C.round(4).to_string(index=False))
    print(rel.round(3).to_string(index=False)); print(f"LOO folds empty: {n_empty}")


def plot():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    d = pd.read_csv(os.path.join(OD, "E9_subject_values.csv"), index_col=0)
    cmap = LinearSegmentedColormap.from_list("blue", ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    panels = [("C1_sphere", "C1 (paced − free), IEAT ACC sphere")]
    if "C1_loo" in d and d.C1_loo.notna().all():
        panels.append(("C1_loo", "C1 (paced − free), leave-one-out H6 cluster"))
    fig, axs = plt.subplots(1, len(panels), figsize=(5.2 * len(panels) + 0.8, 4.4), squeeze=False,
                            facecolor="#fcfcfb")
    for ax, (col, lab) in zip(axs[0], panels):
        ax.set_facecolor("#fcfcfb")
        hm = d.maia.notna()
        sc = ax.scatter(d.loc[hm, col], d.loc[hm, "conf"], c=d.loc[hm, "maia"], cmap=cmap, s=46,
                        edgecolor="#fcfcfb", linewidth=1.5, zorder=3)
        ax.scatter(d.loc[~hm, col], d.loc[~hm, "conf"], s=46, facecolor="none", edgecolor="#52514e",
                   linewidth=1.2, zorder=3, label="no MAIA")
        rho = stats.spearmanr(d[col], d.conf).statistic
        ax.set_title(f"{lab}\nSpearman ρ = {rho:.2f}, N = {len(d)}", fontsize=10, color="#0b0b0b", loc="left")
        ax.set_xlabel("ROI-mean effect (% signal)", color="#52514e"); ax.set_ylabel("mean confidence, paced epochs (1–7)",
                                                                                 color="#52514e")
        ax.axvline(0, color="#d9d8d4", lw=1, zorder=1); ax.grid(color="#ecebe8", lw=0.8); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#b5b4af")
        ax.tick_params(colors="#52514e")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    cb = fig.colorbar(sc, ax=axs[0].tolist(), fraction=0.04, pad=0.02); cb.set_label("MAIA (IEAT total)", color="#52514e")
    fig.suptitle("EXPLORATORY E9 — pacing-related ACC response and confidence", fontsize=11, x=0.01, y=1.06, ha="left")
    os.makedirs(OD, exist_ok=True)
    fig.savefig(os.path.join(OD, "E9_scatter.png"), dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote", os.path.join(OD, "E9_scatter.png"))


if __name__ == "__main__":
    {"behav": behav, "run": run, "plot": plot}[sys.argv[1]]()
