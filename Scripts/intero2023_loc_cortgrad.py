#!/usr/bin/env python
"""
intero2023_loc_cortgrad.py -- EXPLORATORY follow-ups on the localizer (OSF p5ztc), logged in
Preregistration/localizer_deviations.md (rows "Q1 cortical gradients, tier A" and "Q4 breath
tracking") before any alignment or brain-behaviour statistic was computed.

Not to be confused with intero2023_loc_gradient.py (stage C3: posterior -> anterior slope of C1
along the insula on MNI y). This script places whole-cortex condition maps on the Margulies et al.
(2016) principal functional gradients (G1-G3), and fits the reframed Q4 breath-tracking model.

  python Scripts/intero2023_loc_cortgrad.py check    # label-free: atlas grid, parcel coverage, label and
                                                     # centroid correspondence, gradient orientation, Q4 N + VIF
  python Scripts/intero2023_loc_cortgrad.py q1       # parcel maps, positive control (stop rule), spin and
                                                     # BrainSMASH alignment, paired sign-flip tests, profiles
  python Scripts/intero2023_loc_cortgrad.py q4       # FLAME 1 (C1, Intero vs Rest) with tracking PLV +
                                                     # covariates; ROI voxel FWE, whole-brain clusters, mask means
  python Scripts/intero2023_loc_cortgrad.py figure   # Figures/localizer/loc_cortgrad.png

Downloads (never committed; $OUT/group/EXPLORATORY_LOC_cortgrad/downloads/):
  Schaefer-400 7-network, MNI152NLin2009cAsym res-02 (templateflow S3) and fsLR 32k dlabel (CBIG GitHub);
  Margulies 2016 fcgradient01-03, fsLR 32k (neuromaps 0.0.7 / OSF 4mw3a). Spheres: templateflow tpl-fsLR 32k.
Tables: Analysis/results/EXPLORATORY_LOC_cortgrad_*.csv and EXPLORATORY_LOC_Q4_*.csv.
Imaging / parcel outputs: $OUT/group/EXPLORATORY_LOC_cortgrad/ and $OUT/group/LOCq4_<map>/.
"""
import os, sys, json, glob, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
import intero2023_group as G
g = L.g
warnings.simplefilter("ignore")

RD = os.path.join(g.REPO, "Analysis", "results")
QD = os.path.join(g.OUT, "group", "EXPLORATORY_LOC_cortgrad")
DL = os.path.join(QD, "downloads")
VOL_ATLAS = os.path.join(DL, "templateflow", "tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.nii.gz")
VOL_TSV = os.path.join(DL, "templateflow", "tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.tsv")
SURF_ATLAS = os.path.join(DL, "cbig", "Schaefer2018_400Parcels_7Networks_order.dlabel.nii")
GRAD = os.path.join(DL, "neuromaps", "annotations", "margulies2016", "fcgradient0{k}", "fsLR",
                    "source-margulies2016_desc-fcgradient0{k}_space-fsLR_den-32k_hemi-{h}_feature.func.gii")
TF = os.path.expanduser("~/.cache/templateflow/tpl-fsLR")
SPHERE = os.path.join(TF, "tpl-fsLR_hemi-{h}_den-32k_sphere.surf.gii")
MIDTH = os.path.join(DL, "neuromaps", "atlases", "fsLR", "tpl-fsLR_den-32k_hemi-{h}_midthickness.surf.gii")

N_PERM, SEED, N_BOOT = 10000, 20261128, 10000
MIN_FRAC = 0.5          # parcel kept if >= 50% of its atlas voxels lie in the group brain mask
GRADS = ["G1", "G2", "G3"]
PRIMARY_G = ["G1", "G2"]
# orientation rule (label-free, fixed before any condition map was parcellated):
#   G1: Default > mean(Vis, SomMot) (transmodal positive); G2: Vis > SomMot; G3: Cont > Default
ORIENT = {"G1": (["Default"], ["Vis", "SomMot"]), "G2": (["Vis"], ["SomMot"]), "G3": (["Cont"], ["Default"])}
# the nine maps: name -> (source, contrast, participant set)
MAPS = {
    "intero": ("loc", "intero", "intero"), "extero": ("loc", "extero", "all"), "match": ("loc", "match", "all"),
    "C6_self": ("loc", "C6", "all"), "C1": ("loc", "C1", "intero"), "C2": ("loc", "C2", "intero"),
    "C4": ("loc", "C4", "all"), "C5": ("loc", "C5", "intero"), "breath_paced_free": ("breath", "C1", "H1"),
}
# paired condition differences (Fisher z of the alignment): contrast -> (a, b, set)
PAIRED = {"C1": ("intero", "extero", "intero"), "C2": ("intero", "C6_self", "intero"),
          "C4": ("match", "extero", "all"), "C5": ("match", "intero", "intero")}
# Q4
Q4_PARCELS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere"]
Q4_MAPS = {"C1": "C1", "intero": "intero"}
Q4_COVS = ["tracking_plv", "match_plv", "period_cv", "protocol"]
VIF_MAX = 5.0


# =========================================================================== shared inputs
def sets():
    s, _, u = L.participant_sets()
    s = dict(s)
    s["H1"] = json.load(open(os.path.join(RD, "H1_var-base_set-all_subjects.json")))["subjects"]
    return s, u


def map_file(src, con, sub):
    if src == "loc":
        return L.fl_name(sub, "joint", 6, "base", con, "effect")
    return G.fl(sub, "base", con, "effect")          # swc7p breathing-run first level, 6 mm (H1: C1 = paced - free)


def vol_parcels():
    lab = np.asarray(nib.load(VOL_ATLAS).dataobj).astype(int)
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    n_all = np.bincount(lab.ravel(), minlength=401)[1:]
    n_in = np.bincount(lab[gm].ravel(), minlength=401)[1:]
    names = pd.read_csv(VOL_TSV, sep="\t").sort_values("index").name.tolist()
    keep = (n_in / n_all) >= MIN_FRAC
    return lab, gm, n_all, n_in, names, keep


def surf_parcel_gradients():
    """Parcel means of G1-G3 over fsLR 32k Schaefer-400 vertices (medial wall, label 0, excluded;
    vertices where the gradient file is exactly 0 also excluded). Returns (400 x 3) raw values."""
    c = nib.load(SURF_ATLAS)
    lab = np.asarray(c.dataobj)[0].astype(int)            # 64984 = L 32492 + R 32492
    out = np.zeros((400, 3)); nv = np.zeros((400, 3), int)
    for j, k in enumerate((1, 2, 3)):
        v = np.concatenate([nib.load(GRAD.format(k=k, h=h)).darrays[0].data for h in "LR"]).astype(float)
        ok = (lab > 0) & (v != 0)
        out[:, j] = np.bincount(lab[ok], weights=v[ok], minlength=401)[1:] / np.maximum(1, np.bincount(lab[ok], minlength=401)[1:])
        nv[:, j] = np.bincount(lab[ok], minlength=401)[1:]
    names = [c.header.get_axis(0).label[0][i][0] for i in range(1, 401)]
    return out, nv, lab, names


def network(name):
    return name.split("_")[2]


def oriented_gradients():
    raw, nv, _, names = surf_parcel_gradients()
    net = np.array([network(x) for x in names])
    Gm, flips = {}, {}
    for j, gname in enumerate(GRADS):
        hi, lo = ORIENT[gname]
        x = raw[:, j]
        d = x[np.isin(net, hi)].mean() - np.mean([x[net == n].mean() for n in lo])
        flips[gname] = bool(d < 0)
        Gm[gname] = -x if d < 0 else x
    return pd.DataFrame(Gm), flips, nv, names, net


def sphere_centroids():
    """neuromaps 'surface' centroid method on the fsLR 32k sphere: mean vertex coordinate of the
    parcel, snapped to the nearest vertex. Returns (400 x 3) coords and hemiid (0 L, 1 R)."""
    lab = np.asarray(nib.load(SURF_ATLAS).dataobj)[0].astype(int)
    cen, hem = np.zeros((400, 3)), np.zeros(400, int)
    for hi, h in enumerate("LR"):
        xyz = nib.load(SPHERE.format(h=h)).darrays[0].data.astype(float)
        lh = lab[hi * 32492:(hi + 1) * 32492]
        for k in np.unique(lh[lh > 0]):
            m = xyz[lh == k].mean(0)
            cen[k - 1] = xyz[np.argmin(((xyz - m) ** 2).sum(1))]
            hem[k - 1] = hi
    return cen, hem


# =========================================================================== label-free check
def vif(X):
    X = X - X.mean(0); out = []
    for j in range(X.shape[1]):
        y = X[:, j]; Z = np.c_[np.ones(len(y)), np.delete(X, j, 1)]
        b = np.linalg.lstsq(Z, y, rcond=None)[0]
        r2 = 1 - ((y - Z @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        out.append(1 / (1 - r2))
    return np.array(out)


def q4_covariates():
    s, u = sets()
    p = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C1_participants.csv")).set_index("subject")
    b = pd.read_csv(os.path.join(L.QCDIR, "breathing_by_condition.csv"))
    b = b[(b.level == "run") & (b.condition == "intero")].copy(); b["run"] = b.run.astype(int)
    b = b.merge(u[["subject", "run"]], on=["subject", "run"])          # registered (used) runs only
    pv = u.groupby("subject").protocol_version.first()
    D = pd.DataFrame({"tracking_plv": p.trk_tracking_plv, "match_plv": p.loc_match_plv,
                      "period_cv": b.groupby("subject").period_cv.mean()})
    D["protocol"] = (pv.reindex(D.index) == "v1").astype(float) - 0.5
    D["period_cv_n_runs"] = b.groupby("subject").size()
    D = D.reindex(s["intero"]).dropna(subset=Q4_COVS)
    return D


def check():
    rows = []
    add = lambda k, v, note="": rows.append(dict(check=k, value=v, note=note))
    va = nib.load(VOL_ATLAS); gm = nib.load(g.GROUP_MASK)
    add("vol_atlas_grid_equals_glm_grid", bool(va.shape == gm.shape and np.allclose(va.affine, gm.affine)), str(va.shape))
    lab, gmask, n_all, n_in, vnames, keep = vol_parcels()
    add("parcels_kept", int(keep.sum()), f"rule: >= {MIN_FRAC:.0%} of atlas voxels in group brain mask")
    for k in np.where(~keep)[0]:
        add("parcel_dropped", vnames[k], f"in-mask fraction {n_in[k] / n_all[k]:.2f}")
    add("kept_parcel_in_mask_voxels_min_median", f"{n_in[keep].min()} / {np.median(n_in[keep]):.0f}")
    raw, nv, slab, snames = surf_parcel_gradients()
    add("surface_vertices_per_parcel_min", int(nv.min()))
    same_net = all(a.split("_")[1:3] == b.split("_")[1:3] for a, b in zip(vnames, snames))
    add("label_index_hemi_network_match", bool(same_net),
        f"{sum(a != b for a, b in zip(vnames, snames))} names differ only by CBIG's 2020 renaming of sub-labels")
    # centroid correspondence: volumetric centroid (MNI 2009c) vs fsLR S1200 midthickness parcel centroid
    ijk = np.argwhere(lab > 0); xyz = nib.affines.apply_affine(va.affine, ijk); lv = lab[lab > 0]
    vc = np.stack([xyz[lv == k].mean(0) for k in range(1, 401)])
    sc = np.zeros((400, 3))
    for hi, h in enumerate("LR"):
        mxyz = nib.load(MIDTH.format(h=h)).darrays[0].data.astype(float)
        lh = slab[hi * 32492:(hi + 1) * 32492]
        for k in np.unique(lh[lh > 0]):
            sc[k - 1] = mxyz[lh == k].mean(0)
    dist = np.linalg.norm(vc - sc, axis=1)
    add("centroid_distance_mm_median_p95_max", f"{np.median(dist):.1f} / {np.percentile(dist, 95):.1f} / {dist.max():.1f}",
        "volumetric (MNI152NLin2009cAsym) vs fsLR S1200 midthickness parcel centroid; spaces differ slightly")
    nn = np.array([np.argmin(np.linalg.norm(sc - vc[k], axis=1)) == k for k in range(400)])
    add("centroid_nearest_is_same_label", f"{nn.sum()}/400")
    Gd, flips, _, _, net = oriented_gradients()
    for gname in GRADS:
        add(f"{gname}_sign_flipped_by_rule", flips[gname], f"rule {ORIENT[gname][0]} > {ORIENT[gname][1]}")
        for n in ["Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"]:
            add(f"{gname}_network_mean_{n}", round(float(Gd[gname][net == n].mean()), 3))
    rg = Gd[keep].corr(method="spearman")
    add("gradient_intercorrelation_spearman_G1G2_G1G3_G2G3",
        f"{rg.loc['G1', 'G2']:.3f} / {rg.loc['G1', 'G3']:.3f} / {rg.loc['G2', 'G3']:.3f}", "kept parcels")
    # Q4 covariates (no brain data)
    D = q4_covariates()
    vf = vif(D[Q4_COVS].values)
    add("q4_N", len(D), "Intero set ∩ tracking PLV ∩ Match PLV ∩ Intero period_cv")
    add("q4_n_protocol_v1", int((D.protocol > 0).sum()))
    add("q4_period_cv_two_run_participants", int((D.period_cv_n_runs == 2).sum()))
    for c, v in zip(Q4_COVS, vf):
        add(f"q4_vif_{c}", round(float(v), 3))
    add("q4_match_plv_dropped", bool(vf.max() > VIF_MAX), f"rule: drop Match PLV if any VIF > {VIF_MAX}")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_checks.csv"), index=False)
    print(R.to_string())


# =========================================================================== Q1
def parcel_matrix(src, con, subs, lab, gm, keep):
    li = lab[gm]
    cnt = np.bincount(li, minlength=401)[1:]
    def one(x):
        d = nib.load(map_file(src, con, x)).get_fdata(dtype=np.float32)[gm].astype(float)
        return (np.bincount(li, weights=d, minlength=401)[1:] / np.maximum(cnt, 1))[keep]
    with ThreadPoolExecutor(8) as ex:
        return np.stack(list(ex.map(one, subs)))


def ranks(a, axis=-1):
    return stats.rankdata(a, axis=axis)


def spearman_rows(Xr, gr):
    """Pearson of ranked rows of Xr with ranked vector gr."""
    Xc = Xr - Xr.mean(-1, keepdims=True); gc = gr - gr.mean()
    return (Xc @ gc) / (np.linalg.norm(Xc, axis=-1) * np.linalg.norm(gc))


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o]); out = np.empty(m); out[o] = np.minimum(1, adj)
    return out


def signflip_two(v, rng):
    v = np.asarray(v); obs = abs(v.mean())
    null = np.abs((rng.choice([-1.0, 1.0], size=(N_PERM, len(v))) * v).mean(1))
    return (1 + (null >= obs - 1e-15).sum()) / (N_PERM + 1)


def boot_ci(v, rng):
    b = v[rng.integers(0, len(v), size=(N_BOOT, len(v)))].mean(1)
    return np.quantile(b, .025), np.quantile(b, .975)


def q1():
    from neuromaps.nulls.spins import gen_spinsamples
    os.makedirs(QD, exist_ok=True)
    s, _ = sets()
    lab, gm, _, _, vnames, keep = vol_parcels()
    Gd, flips, _, snames, net = oriented_gradients()
    Gk = Gd[keep].reset_index(drop=True); netk = net[keep]
    # parcel maps (cached)
    P = {}
    for name, (src, con, st) in MAPS.items():
        f = os.path.join(QD, f"parcels_{name}.npz")
        subs = list(s[st])
        if os.path.exists(f) and list(np.load(f)["subjects"]) == subs:
            P[name] = (subs, np.load(f)["X"])
        else:
            X = parcel_matrix(src, con, subs, lab, gm, keep)
            np.savez_compressed(f, X=X, subjects=np.array(subs)); P[name] = (subs, X)
        print(name, P[name][1].shape, flush=True)
    grp = {k: X.mean(0) for k, (_, X) in P.items()}
    pd.DataFrame(dict(parcel=np.array(vnames)[keep], network=netk, **{k: v for k, v in grp.items()},
                      **{g_: Gk[g_].values for g_ in GRADS})).to_csv(os.path.join(QD, "group_parcel_maps.csv"), index=False)
    # spin permutations (Vasa), parcel level, fsLR 32k sphere, mirrored rotations
    sf = os.path.join(QD, f"spins_vasa_{N_PERM}_seed{SEED}.npy")
    if os.path.exists(sf):
        spins = np.load(sf)
    else:
        cen, hem = sphere_centroids()
        spins = gen_spinsamples(cen[keep], hem[keep], n_rotate=N_PERM, method="vasa", seed=SEED)
        np.save(sf, spins)
    gr = {g_: ranks(Gk[g_].values) for g_ in GRADS}
    # positive control first (stop rule): Self vs Rest (C6) toward the transmodal end of G1
    def align(name, g_):
        x = grp[name]; rx = ranks(x)
        r = spearman_rows(rx[None], gr[g_])[0]
        null = spearman_rows(rx[spins.T], gr[g_])
        return r, null
    r6, n6 = align("C6_self", "G1")
    p6 = (1 + (np.abs(n6) >= abs(r6) - 1e-12).sum()) / (N_PERM + 1)
    passed = bool(r6 > 0 and p6 < 0.05)
    ctrl = pd.DataFrame([dict(control="C6_self toward transmodal end of G1", N=len(P["C6_self"][0]), rho=r6, p_spin=p6,
                              rule="rho > 0 and two-sided spin p < .05", passed=passed)])
    ctrl.to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_control.csv"), index=False)
    print(ctrl.to_string(), flush=True)
    if not passed:
        print("POSITIVE CONTROL FAILED: stop rule applied; no further Q1 statistic computed.")
        return
    # group-map alignment: spin and BrainSMASH nulls
    from brainsmash.mapgen.base import Base
    ijk = np.argwhere((lab > 0) & gm); xyz = nib.affines.apply_affine(nib.load(VOL_ATLAS).affine, ijk)
    lv = lab[(lab > 0) & gm]
    vc = np.stack([xyz[lv == k].mean(0) for k in range(1, 401)])[keep]
    D = np.linalg.norm(vc[:, None] - vc[None], axis=-1)
    rows = []
    for name in MAPS:
        x = grp[name]
        sur = Base(x=x, D=D, seed=SEED)(n=N_PERM)                 # variogram-matched surrogates of the map
        rs = ranks(sur, axis=1)
        for g_ in GRADS:
            r, null = align(name, g_)
            nsm = spearman_rows(rs, gr[g_])
            rows.append(dict(map=name, gradient=g_, N=len(P[name][0]), n_parcels=int(keep.sum()), rho=r,
                             p_spin=(1 + (np.abs(null) >= abs(r) - 1e-12).sum()) / (N_PERM + 1),
                             spin_null_lo=np.quantile(null, .025), spin_null_hi=np.quantile(null, .975),
                             p_smash=(1 + (np.abs(nsm) >= abs(r) - 1e-12).sum()) / (N_PERM + 1),
                             smash_null_lo=np.quantile(nsm, .025), smash_null_hi=np.quantile(nsm, .975),
                             role="primary" if g_ in PRIMARY_G else "descriptive",
                             group_mean_psc=float(x.mean())))
        print("aligned", name, flush=True)
    A = pd.DataFrame(rows)
    for g_ in GRADS:
        m = A.gradient == g_
        A.loc[m, "p_spin_holm9"] = holm(A.loc[m, "p_spin"]); A.loc[m, "p_smash_holm9"] = holm(A.loc[m, "p_smash"])
    A.to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_alignment.csv"), index=False)
    # participant alignment (Fisher z of Spearman) and paired condition differences
    rng = np.random.default_rng(SEED)
    Z, srows = {}, []
    for name, (subs, X) in P.items():
        for g_ in GRADS:
            z = np.arctanh(spearman_rows(ranks(X, axis=1), gr[g_]))
            Z[(name, g_)] = pd.Series(z, index=subs)
            srows += [dict(subject=a, map=name, gradient=g_, fisher_z=b) for a, b in zip(subs, z)]
    pd.DataFrame(srows).to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_subjects.csv"), index=False)
    prow = []
    for con, (a, b, st) in PAIRED.items():
        subs = list(s[st])
        for g_ in GRADS:
            d = (Z[(a, g_)].reindex(subs) - Z[(b, g_)].reindex(subs)).values
            assert np.isfinite(d).all()
            lo, hi = boot_ci(d, rng); tt = stats.ttest_1samp(d, 0)
            prow.append(dict(contrast=con, a=a, b=b, gradient=g_, N=len(d), mean_z_a=Z[(a, g_)].reindex(subs).mean(),
                             mean_z_b=Z[(b, g_)].reindex(subs).mean(), diff_mean=d.mean(), ci_lo=lo, ci_hi=hi,
                             dz=d.mean() / d.std(ddof=1), t=tt.statistic, p_t=tt.pvalue, p_signflip=signflip_two(d, rng),
                             role="family" if g_ in PRIMARY_G else "descriptive"))
    Pd = pd.DataFrame(prow)
    fam = Pd.role == "family"
    Pd.loc[fam, "p_signflip_holm8"] = holm(Pd.loc[fam, "p_signflip"])
    Pd.to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_paired.csv"), index=False)
    # profiles: 10 gradient bins (deciles of parcels), mean PSC, participant-bootstrap CI
    pr = []
    for g_ in GRADS:
        order = np.argsort(Gk[g_].values); bins = np.array_split(order, 10)
        for name, (subs, X) in P.items():
            for bi, idx in enumerate(bins):
                v = X[:, idx].mean(1); lo, hi = boot_ci(v, rng)
                pr.append(dict(gradient=g_, bin=bi + 1, gradient_mean=Gk[g_].values[idx].mean(), n_parcels=len(idx),
                               map=name, N=len(v), psc_mean=v.mean(), ci_lo=lo, ci_hi=hi))
    pd.DataFrame(pr).to_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_profiles.csv"), index=False)
    print(A.round(4).to_string()); print(Pd.round(4).to_string())


# =========================================================================== Q4
def q4_dirs(name):
    return os.path.join(g.OUT, "group", f"LOCq4_{name}")


def q4_design():
    D = q4_covariates()
    covs = list(Q4_COVS)
    if vif(D[covs].values).max() > VIF_MAX:
        covs.remove("match_plv")
    X = pd.DataFrame({"intercept": 1.0}, index=D.index)
    for c in covs:
        X[c] = D[c] - D[c].mean()
    return X


def q4_flame(name):
    con = Q4_MAPS[name]; od = q4_dirs(name); os.makedirs(od, exist_ok=True)
    X = q4_design(); subs = list(X.index)
    ref = nib.load(g.GROUP_MASK)
    stack = lambda kind: nib.Nifti1Image(np.stack([nib.load(L.fl_name(x, "joint", 6, "base", con, kind)).get_fdata(dtype=np.float32)
                                                   for x in subs], -1), ref.affine)
    stack("effect").to_filename(f"{od}/cope.nii.gz"); stack("variance").to_filename(f"{od}/varcope.nii.gz")
    dofs = [json.load(open(os.path.join(L.fl_dir(x), f"{x}_task-localizer_firstlevel.json")))["fits"]["joint:sm6:base"]["dof"] for x in subs]
    m = ref.get_fdata() > 0
    nib.Nifti1Image(np.stack([m * d for d in dofs], -1).astype(np.float32), ref.affine).to_filename(f"{od}/dof.nii.gz")
    n, p = X.shape
    open(f"{od}/design.mat", "w").write(f"/NumWaves {p}\n/NumPoints {n}\n/Matrix\n" +
                                        "".join(" ".join(f"{v:.6f}" for v in row) + "\n" for row in X.values))
    c = " ".join("1" if col == "tracking_plv" else "0" for col in X.columns)
    open(f"{od}/design.con", "w").write(f"/NumWaves {p}\n/NumContrasts 1\n/Matrix\n{c}\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    sm = G.sh(f"smoothest -d {n - p} -r {od}/stats/res4d.nii.gz -m {g.GROUP_MASK}")
    meta = dict(label=f"LOCq4_{name}", con=con, subjects=subs, design=list(X.columns), test="tracking_plv", dofs=dofs,
                smoothest={l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()})
    os.remove(f"{od}/stats/res4d.nii.gz")
    json.dump(meta, open(f"{od}/loc_meta.json", "w"), indent=1)
    return name


def q4_clusters(name):
    from intero2023_loc_group import atlas_labels
    od = q4_dirs(name); st = f"{od}/stats"; smv = json.load(open(f"{od}/loc_meta.json"))["smoothest"]
    out = []
    for sgn, zf in ((+1, f"{st}/zstat1.nii.gz"), (-1, f"{st}/zneg.nii.gz")):
        if sgn < 0:
            G.sh(f"fslmaths {st}/zstat1.nii.gz -mul -1 {zf}")
        txt = G.sh(f"cluster -i {zf} -t 3.09 -p 1 -d {smv['DLH']} --volume={int(smv['VOLUME'])} --mm")
        lines = [l.split("\t") for l in txt.strip().splitlines()]
        for l in lines[1:]:
            d = dict(zip(lines[0], l))
            mni = [float(d["MAX X (mm)"]), float(d["MAX Y (mm)"]), float(d["MAX Z (mm)"])]
            out.append(dict(map=name, sign="+" if sgn > 0 else "-", size_vox=int(d["Voxels"]), peak_z=sgn * float(d["MAX"]),
                            peak_mni=mni, p_cluster=float(d["P"]), p_fwe_two_sided=min(1.0, 2 * float(d["P"]))))
    c = pd.DataFrame(out, columns=["map", "sign", "size_vox", "peak_z", "peak_mni", "p_cluster", "p_fwe_two_sided"])
    if len(c):
        c = c.sort_values("size_vox", ascending=False)
        c["labels"] = [atlas_labels(m) if (p < 0.05 or i < 5) else "" for i, (m, p) in enumerate(zip(c.peak_mni, c.p_fwe_two_sided))]
    return c


def q4_maskmean(y, X, rng):
    """OLS of mask mean on the design; tracking term: slope, t, p, partial r, Freedman-Lane permutation p,
    JZS BF10 = BF(full) / BF(reduced) (Zellner-Siow, r = .707 and .354)."""
    from intero2023_loc_group import jzs_reg_bf10
    Xf = X.values; j = list(X.columns).index("tracking_plv"); Xr = np.delete(Xf, j, 1)
    n, p = Xf.shape
    def fit(A, yy):
        b = np.linalg.lstsq(A, yy, rcond=None)[0]; e = yy - A @ b; return b, e
    b, e = fit(Xf, y); s2 = (e @ e) / (n - p); cov = s2 * np.linalg.inv(Xf.T @ Xf)
    t0 = b[j] / np.sqrt(cov[j, j])
    br, er = fit(Xr, y)
    xr = fit(Xr, Xf[:, j])[1]
    pr = stats.pearsonr(xr, er)[0]
    fl = []
    for _ in range(N_PERM):
        ys = Xr @ br + rng.permutation(er)
        bb, ee = fit(Xf, ys); ss = (ee @ ee) / (n - p)
        fl.append(bb[j] / np.sqrt(ss * np.linalg.inv(Xf.T @ Xf)[j, j]))
    fl = np.abs(np.array(fl))
    r2 = lambda ee: 1 - (ee @ ee) / ((y - y.mean()) @ (y - y.mean()))
    bf = {}
    for rs in (0.707, 0.354):
        bf[rs] = jzs_reg_bf10(r2(e), n, r=rs, p=p - 1) / jzs_reg_bf10(r2(er), n, r=rs, p=p - 2)
    return dict(slope=b[j], t=t0, p_param=2 * stats.t.sf(abs(t0), n - p), partial_r=pr,
                p_perm_freedman_lane=(1 + (fl >= abs(t0) - 1e-12).sum()) / (N_PERM + 1),
                bf10=bf[0.707], bf01=1 / bf[0.707], bf10_r354=bf[0.354],
                zero_order_spearman=stats.spearmanr(X["tracking_plv"], y)[0])


def q4():
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(q4_flame, Q4_MAPS))
    X = q4_design(); rng = np.random.default_rng(SEED)
    D = q4_covariates()
    roi_rows, wb, mm, subj = [], [], [], []
    for name in Q4_MAPS:
        od = q4_dirs(name); meta = json.load(open(f"{od}/loc_meta.json")); assert meta["subjects"] == list(X.index)
        z = nib.load(f"{od}/stats/zstat1.nii.gz").get_fdata(); aff = nib.load(f"{od}/stats/zstat1.nii.gz").affine
        cope = nib.load(f"{od}/cope.nii.gz").get_fdata(dtype=np.float32); smv = meta["smoothest"]
        for roi in Q4_PARCELS:
            rm = nib.load(G.ROI(roi)).get_fdata() > 0
            zz = np.where(rm, np.abs(z), -np.inf); ijk = np.unravel_index(np.argmax(zz), z.shape); zt = float(zz[ijk])
            resels = rm.sum() / smv["RESELS"]
            pv = min(1.0, 2 * G.grf_p(zt, resels)) if zt > 1 else 1.0
            zthr = float(G.sh(f"ptoz {0.05 / (5 * 2)} -g {resels}").strip())
            y = cope[rm].mean(0).astype(float)
            rec = dict(map=name, roi=roi, N=len(y), peak_z=float(z[ijk]),
                       peak_mni=[round(float(a), 1) for a in nib.affines.apply_affine(aff, ijk)], resels=resels,
                       p_fwe=pv, p_fwe_bonf5=min(1.0, 5 * pv), z_thresh_bonf5=zthr, n_vox_sig=int((np.abs(z[rm]) > zthr).sum()),
                       flame_ok=bool(np.isfinite(z[rm]).all() and (z[rm] != 0).mean() > 0.99))
            roi_rows.append(rec)
            mm.append(dict(map=name, roi=roi, N=len(y), mask_mean=y.mean(), **q4_maskmean(y, X, rng)))
            subj += [dict(map=name, roi=roi, subject=s_, mask_mean=v) for s_, v in zip(X.index, y)]
        wb.append(q4_clusters(name))
        print("q4", name, flush=True)
    R = pd.DataFrame(roi_rows); M = pd.DataFrame(mm)
    R["p_fwe_bonf10"] = np.minimum(1, 10 * R.p_fwe)
    R.to_csv(os.path.join(RD, "EXPLORATORY_LOC_Q4_roi.csv"), index=False)
    M.to_csv(os.path.join(RD, "EXPLORATORY_LOC_Q4_maskmeans.csv"), index=False)
    pd.concat(wb).to_csv(os.path.join(RD, "EXPLORATORY_LOC_Q4_wholebrain.csv"), index=False)
    cv = D[Q4_COVS].describe().T[["count", "mean", "std", "min", "50%", "max"]]
    cv["vif"] = vif(D[Q4_COVS].values); cv["in_model"] = [c in X.columns for c in Q4_COVS]
    cv.to_csv(os.path.join(RD, "EXPLORATORY_LOC_Q4_covariates.csv"))
    print(R.round(4).to_string()); print(M.round(4).to_string()); print(pd.concat(wb).head(20).to_string())


# =========================================================================== figure
def figure():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    COL = {"intero": "#2a78d6", "extero": "#eb6834", "match": "#1baf7a", "C6_self": "#eda100"}
    LAB = {"intero": "Intero", "extero": "Extero", "match": "Match", "C6_self": "Self", "C1": "C1 Intero − Extero",
           "C2": "C2 Intero − Self", "C4": "C4 Match − Extero", "C5": "C5 Match − Intero", "breath_paced_free": "Breathing paced − free"}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False})
    Pr = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_profiles.csv"))
    A = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_alignment.csv"))
    Pd = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_cortgrad_paired.csv"))
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.8), facecolor=SURF)
    for a in ax.ravel():
        a.set_facecolor(SURF); a.grid(axis="y", color="#e6e5e1", lw=0.6); a.set_axisbelow(True)
    pole = {"G1": ("unimodal", "transmodal"), "G2": ("somatomotor", "visual")}
    for k, g_ in enumerate(PRIMARY_G):
        a = ax[0, k]
        ends = {}
        for name, col in COL.items():
            d = Pr[(Pr.gradient == g_) & (Pr["map"] == name)]
            a.fill_between(d.bin, d.ci_lo, d.ci_hi, color=col, alpha=0.15, lw=0)
            a.plot(d.bin, d.psc_mean, color=col, lw=2, marker="o", ms=4)
            ends[name] = d.psc_mean.iloc[-1]
        ys = sorted(ends.items(), key=lambda kv: kv[1]); placed = []
        for name, y in ys:                      # direct labels, nudged apart by >= 0.022 PSC
            y = max(y, placed[-1] + 0.022) if placed else y; placed.append(y)
            a.annotate(LAB[name], (10.15, y), va="center", color=INK, fontsize=8.5)
        a.axhline(0, color=INK2, lw=0.8)
        a.set_xticks(range(1, 11)); a.set_xlim(0.6, 11.9)
        a.set_xlabel(f"{g_} decile  ({pole[g_][0]} → {pole[g_][1]})"); a.set_ylabel("Mean PSC vs Rest")
        a.set_title(f"{'AB'[k]}. Condition profile along {g_} (95% participant-bootstrap CI)", loc="left", fontsize=9.5)
    # alignment
    a = ax[1, 0]; order = list(MAPS)
    for k, g_ in enumerate(PRIMARY_G):
        d = A[A.gradient == g_].set_index("map").reindex(order); yy = np.arange(len(order)) + (k - 0.5) * 0.3
        a.hlines(yy, d.spin_null_lo, d.spin_null_hi, color="#c9c7c1", lw=4)
        a.plot(d.rho, yy, "o" if g_ == "G1" else "s", color=INK if g_ == "G1" else INK2, ms=6,
               mfc=INK if g_ == "G1" else SURF, label=f"{g_} (grey bar: 95% spin null)")
        for yi, r, p in zip(yy, d.rho, d.p_spin_holm9):
            if p < .05:
                a.annotate("*", (r, yi), xytext=(6, -3), textcoords="offset points", fontsize=10)
    a.axvline(0, color=INK2, lw=0.8); a.set_yticks(range(len(order))); a.set_yticklabels([LAB[o] for o in order])
    a.invert_yaxis(); a.set_xlabel("Spearman ρ, group parcel map × gradient")
    a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    a.set_title("C. Group-map alignment (* spin p Holm ×9 < .05)", loc="left", fontsize=9.5); a.grid(False)
    # paired
    a = ax[1, 1]; cons = list(PAIRED)
    for k, g_ in enumerate(PRIMARY_G):
        d = Pd[Pd.gradient == g_].set_index("contrast").reindex(cons); yy = np.arange(len(cons)) + (k - 0.5) * 0.3
        a.hlines(yy, d.ci_lo, d.ci_hi, color=INK if g_ == "G1" else INK2, lw=1.5)
        a.plot(d.diff_mean, yy, "o" if g_ == "G1" else "s", color=INK if g_ == "G1" else INK2, ms=6,
               mfc=INK if g_ == "G1" else SURF, label=g_)
        for yi, r, hi, p in zip(yy, d.diff_mean, d.ci_hi, d.p_signflip_holm8):
            if p < .05:
                a.annotate("*", (hi, yi), xytext=(4, -3), textcoords="offset points", fontsize=10)
    a.axvline(0, color=INK2, lw=0.8); a.set_yticks(range(len(cons)))
    a.set_yticklabels([f"{c}: {LAB[PAIRED[c][0]]} − {LAB[PAIRED[c][1]]}" for c in cons]); a.invert_yaxis()
    a.set_xlabel("Δ Fisher z of participant alignment (95% bootstrap CI)")
    a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    a.set_title("D. Paired differences in alignment (* sign-flip Holm ×8 < .05)", loc="left", fontsize=9.5); a.grid(False)
    fig.suptitle("Localizer maps on the Margulies (2016) cortical gradients (Schaefer-400, cortex only; exploratory)",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(g.REPO, "Figures", "localizer", "loc_cortgrad.png")
    fig.savefig(out, dpi=200, facecolor=SURF); print(out)


if __name__ == "__main__":
    {"check": check, "q1": q1, "q4": q4, "figure": figure}[sys.argv[1]]()
