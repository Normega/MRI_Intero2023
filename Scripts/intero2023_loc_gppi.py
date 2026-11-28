#!/usr/bin/env python
"""
intero2023_loc_gppi.py -- localizer stage C2 (EXPLORATORY; OSF p5ztc plan §6): generalised
psychophysiological interaction (gPPI; McLaren et al. 2012) from the gate-primary fROIs and
the a-priori masks, contrasts C1 (Intero - Extero) and C4 (Match - Extero), whole brain.
Choices: Analysis/localizer_stageC2_implementation_notes.md (committed before any statistic).

  python Scripts/intero2023_loc_gppi.py check sub-067   # design-level checks (no statistic)
  python Scripts/intero2023_loc_gppi.py subject sub-060 # first levels, one participant
  python Scripts/intero2023_loc_gppi.py flames          # FLAME 1 per seed x contrast
  python Scripts/intero2023_loc_gppi.py tables          # clusters, within-network summaries
  python Scripts/intero2023_loc_gppi.py subjectrvt sub-060  # +RVT refit, seeds with a surviving cluster
  python Scripts/intero2023_loc_gppi.py rvtcompare      # cluster means, base vs +RVT, same participants

First level = the stage-A joint two-run 6 mm base model (intero2023_loc_common.design,
build_joint_design; intero2023_firstlevel.fit, AR(1); PSC per run) plus, per seed:
  seed       the seed's mean 6 mm PSC signal, adjusted per run for the run's non-task
             columns (nuisance sheet without cosine_*, 200 s DCT, intercept); run-specific
  ppi_<c>    for every modelled condition c (intero, extero, match, self, instruction,
             intero_nc): deconvolved seed x 0/1 boxcar of c on a 16-bin microtime grid,
             convolved with the SPM canonical HRF, sampled at TR (k + 0.5); shared across runs
Deconvolution (p5ztc Transformations 13): neural signal = DCT basis (n_vol components) on
the microtime grid; ridge solution of y = H b, lambda by generalised cross-validation.
Outputs (outside the repo): $OUT/gppi/<sub>/, $OUT/group/LOC_gPPI_<seed>_<con>/.
"""
import os, sys, json, time, warnings, functools
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
from scipy.signal import fftconvolve
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
from nilearn.glm.first_level.hemodynamic_models import spm_hrf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
import intero2023_loc_group as LG
import intero2023_group as G
from intero2023_firstlevel import fit
import intero2023_loc_connectivity as LC
g = L.g
warnings.simplefilter("ignore")
NT = 16                                         # microtime bins per TR
PARCELS = LC.PARCELS
SEEDS = [f"froi{LC.SHORT[p]}" for p in PARCELS] + [f"mask{LC.SHORT[p]}" for p in PARCELS]
GCON = {"C1": {"ppi_intero": 1.0, "ppi_extero": -1.0}, "C4": {"ppi_match": 1.0, "ppi_extero": -1.0}}
PPI = [f"ppi_{c}" for c in L.TASK]
GPD = os.path.join(g.OUT, "gppi")
RD = LC.RD
NETS = LC.NETS
LAMBDA_GRID = 10.0 ** np.linspace(-6, 1, 71)    # x s_max^2


def gname(sub, seed, con, kind, var="base"):
    v = "" if var == "base" else f"_var-{var}"
    return os.path.join(GPD, sub, f"{sub}_task-localizer_gppi_seed-{seed}{v}_con-{con}_{kind}.nii.gz")


@functools.lru_cache(maxsize=None)
def seed_masks(sub):
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    out = {f"mask{LC.SHORT[p]}": nib.load(G.ROI(p)).get_fdata() > 0 for p in PARCELS}
    fr = json.load(open(os.path.join(g.OUT, "froi", "froi_index.json")))["subjects"]
    if sub in fr:
        out.update({f"froi{k}": v for k, v in LC.froi_masks(sub).items()})
    for k, v in out.items():
        assert (v <= gm).all() and v.sum() >= 50, (sub, k)
    return out


# --------------------------------------------------------------------------- deconvolution
def hrf_micro():
    return spm_hrf(g.TR, oversampling=NT)        # sampled every TR / NT s, sums to 1


def sample_idx(n):
    return np.arange(n) * NT + NT // 2           # microtime index of TR (k + 0.5)


def conv(x, h):
    return fftconvolve(x, h)[: len(x)]


def dct_basis(M, K):
    j = np.arange(M)[:, None]; k = np.arange(K)[None, :]
    B = np.cos(np.pi * (2 * j + 1) * k / (2 * M))
    return B / np.linalg.norm(B, axis=0)


def deconvolve(y, n):
    """Ridge deconvolution (Gitelman et al. 2003 style). Returns microtime neural signal and
    diagnostics (lambda ratio, r(reconvolved, y))."""
    M = n * NT; h = hrf_micro(); si = sample_idx(n)
    B = dct_basis(M, n)
    H = np.column_stack([conv(B[:, k], h)[si] for k in range(n)])
    U, s, Vt = np.linalg.svd(H, full_matrices=False)
    uy = U.T @ y
    best = None
    for lr in LAMBDA_GRID:
        lam = lr * s[0] ** 2
        f = s ** 2 / (s ** 2 + lam)
        yhat = U @ (f * uy)
        rss = float(((y - yhat) ** 2).sum())
        gcv = n * rss / (n - f.sum()) ** 2
        if best is None or gcv < best[0]:
            best = (gcv, lr, lam)
    _, lr, lam = best
    b = Vt.T @ ((s / (s ** 2 + lam)) * uy)
    x = B @ b
    yhat = H @ b
    return x, dict(lambda_ratio=float(lr), r_reconv=float(np.corrcoef(yhat, y)[0, 1]),
                   edge=bool(lr in (LAMBDA_GRID[0], LAMBDA_GRID[-1])))


def boxcar(ev, cond, n):
    t = np.arange(n * NT) * g.TR / NT
    b = np.zeros(n * NT)
    for e in ev[ev.trial_type == cond].itertuples():
        b[(t >= e.onset) & (t < e.onset + e.duration)] = 1.0
    return b


def run_ppi(y, dm, ev, n):
    """Seed adjustment, deconvolution and PPI columns for one run's design."""
    task = [c for c in L.TASK if c in dm.columns]
    Xn = dm.drop(columns=task).values.astype(float)
    yadj = y - Xn @ np.linalg.lstsq(Xn, y, rcond=None)[0]
    x, diag = deconvolve(yadj, n)
    h = hrf_micro(); si = sample_idx(n)
    cols = {"seed": yadj}
    psy_r = {}
    for c in task:
        bc = boxcar(ev, c, n)
        cols[f"ppi_{c}"] = conv(x * bc, h)[si]
        psy_r[c] = float(np.corrcoef(conv(bc, h)[si], dm[c].values)[0, 1])   # own HRF vs nilearn regressor
    return pd.DataFrame(cols, index=dm.index), diag, psy_r


def seed_designs(sub, runs, Y, smask_idx, physio=()):
    dms = L.run_designs(sub, runs, physio)
    out, diags = [], []
    for (k, dm), r in zip(dms, runs.itertuples()):
        assert k == int(r.run)
        ev, _ = L.read_events(sub, r.run, int(r.n_vol))
        y = Y[k][:, smask_idx].mean(axis=1).astype(np.float64)
        cols, diag, psy = run_ppi(y, dm, ev, int(r.n_vol))
        out.append(pd.concat([dm, cols], axis=1)); diags.append(dict(run=k, **diag, psych_r_min=min(psy.values())))
    return L.build_joint_design(out, task=L.TASK + PPI), diags


# --------------------------------------------------------------------------- first level
def load_data(sub, runs):
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    Y = {}
    for r in runs.itertuples():
        y = masker.transform(g.bold_path(sub, L.entity(r.run)))
        assert y.shape[0] == int(r.n_vol)
        Y[int(r.run)], _ = mean_scaling(y, axis=0)
    return masker, Y


def subject(sub, var="base", seeds=None):
    """var 'rvt': +rvt_rrf_z in the model (and in the seed adjustment), for `seeds` only
    (the pre-logged rule: seeds with a surviving cluster)."""
    t0 = time.time()
    runs = L.run_list(); runs = runs[runs.subject == sub]
    physio = tuple(L.PHYSIO_RVT) if var == "rvt" else ()
    if var == "rvt":
        assert runs.has_rvt.all(), sub
    os.makedirs(os.path.join(GPD, sub), exist_ok=True)
    masker, Y = load_data(sub, runs)
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    log = dict(subject=sub, seeds={})
    Ym = np.vstack([Y[int(r.run)] for r in runs.itertuples()])
    for seed, m in seed_masks(sub).items():
        if seeds is not None and seed not in seeds:
            continue
        si = np.flatnonzero(m[gm])
        dm, diags = seed_designs(sub, runs, Y, si, physio)
        specs = {c: s for c, s in GCON.items() if L.estimable(dm, s)}
        res, dof = fit(Ym, dm, specs)
        for c, (e, v) in res.items():
            if not (v > 0).all():
                raise ValueError(f"{sub} {seed} {c}: non-positive variance")
            for kind, arr in (("effect", e), ("variance", v), ("t", (e / np.sqrt(v)).astype(np.float32))):
                g.save_float32(masker.inverse_transform(arr), gname(sub, seed, c, kind, var))
        log["seeds"][seed] = dict(n_vox=int(m.sum()), dof=dof, n_reg=int(dm.shape[1]), contrasts=sorted(res),
                                  runs=diags)
    log["seconds"] = round(time.time() - t0, 1)
    js = f"{sub}_gppi.json" if var == "base" else f"{sub}_gppi_var-{var}.json"
    json.dump(log, open(os.path.join(GPD, sub, js), "w"), indent=1)
    print(f"{sub} gPPI done in {log['seconds']} s", flush=True)


def check(sub):
    """Design-level checks on one participant (no voxelwise statistic):
    (a) own microtime HRF convolution of each boxcar vs nilearn's regressor;
    (b) deconvolution diagnostics; (c) joint design rank and estimability;
    (d) recovery of a known neural signal from a simulated seed."""
    runs = L.run_list(); runs = runs[runs.subject == sub]
    masker, Y = load_data(sub, runs)
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    out = dict(subject=sub, seeds={})
    for seed, m in seed_masks(sub).items():
        dm, diags = seed_designs(sub, runs, Y, np.flatnonzero(m[gm]))
        X = dm.values.astype(float)
        out["seeds"][seed] = dict(runs=diags, rank=int(np.linalg.matrix_rank(X)), n_col=int(X.shape[1]),
                                  estimable={c: L.estimable(dm, s) for c, s in GCON.items()},
                                  max_abs_corr_ppi_vs_task=float(max(abs(np.corrcoef(dm[f"ppi_{c}"], dm[c])[0, 1])
                                                                     for c in ("intero", "extero", "match")
                                                                     if c in dm.columns)))
    # (d) simulation: known neural signal -> BOLD (+ noise) -> deconvolution
    rng = np.random.default_rng(20261127); n = 196
    from scipy.ndimage import gaussian_filter1d
    xt = gaussian_filter1d(rng.standard_normal(n * NT), sigma=NT * 1.5)
    y = conv(xt, hrf_micro())[sample_idx(n)]
    sim = {}
    for snr in (np.inf, 4.0, 1.0):
        yy = y + (0 if np.isinf(snr) else rng.standard_normal(n) * y.std() / snr)
        xh, dg = deconvolve(yy - yy.mean(), n)
        lowpass = lambda a: gaussian_filter1d(a, NT)            # compare at ~TR resolution
        sim[str(snr)] = dict(r_neural=float(np.corrcoef(lowpass(xh), lowpass(xt))[0, 1]), **dg)
    out["simulation"] = sim
    json.dump(out, open(os.path.join(g.REPO, "Analysis", "localizer_gppi_design_check.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


# --------------------------------------------------------------------------- group
def gsets():
    s, _, _ = L.participant_sets()
    fr = json.load(open(os.path.join(g.OUT, "froi", "froi_index.json")))["subjects"]
    return {("mask", "C1"): s["intero"], ("mask", "C4"): s["all"], ("froi", "C1"): fr, ("froi", "C4"): fr}


def labels():
    return [(seed, con) for seed in SEEDS for con in GCON]


def fit_label(seed, con):
    return f"gPPI_{seed}_{con}"


def run_flame(sc):
    seed, con = sc
    subs = gsets()[(seed[:4], con)]
    od = LG.fit_dir(fit_label(seed, con)); os.makedirs(od, exist_ok=True)
    ref = nib.load(g.GROUP_MASK); m = ref.get_fdata() > 0
    dof = [json.load(open(os.path.join(GPD, x, f"{x}_gppi.json")))["seeds"][seed]["dof"] for x in subs]
    stack = lambda kind: nib.Nifti1Image(np.stack([nib.load(gname(x, seed, con, kind)).get_fdata(dtype=np.float32)
                                                   for x in subs], -1), ref.affine)
    stack("effect").to_filename(f"{od}/cope.nii.gz"); stack("variance").to_filename(f"{od}/varcope.nii.gz")
    nib.Nifti1Image(np.stack([m * d for d in dof], -1).astype(np.float32), ref.affine).to_filename(f"{od}/dof.nii.gz")
    n = len(subs)
    open(f"{od}/design.mat", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1.000000\n" * n)
    open(f"{od}/design.con", "w").write("/NumWaves 1\n/NumContrasts 1\n/Matrix\n1\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    sm = G.sh(f"smoothest -d {n - 1} -r {od}/stats/res4d.nii.gz -m {g.GROUP_MASK}")
    meta = dict(label=fit_label(seed, con), seed=seed, con=con, subjects=subs, dofs=dof,
                smoothest={l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()})
    os.remove(f"{od}/stats/res4d.nii.gz")
    json.dump(meta, open(f"{od}/loc_meta.json", "w"), indent=1)
    return fit_label(seed, con)


def flames():
    with ThreadPoolExecutor(int(os.environ.get("FLAME_JOBS", "10"))) as ex:
        for lab in ex.map(run_flame, labels()):
            print(f"{lab}: fitted", flush=True)


def tables():
    cl, nw = [], []
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    netm = {n: nib.load(os.path.join(g.REPO, "ROIs", f"net_{n}.nii")).get_fdata() > 0 for n in NETS}
    for seed, con in labels():
        lab = fit_label(seed, con)
        c = LG.clusters(lab, f"gPPI {con} seed {seed}")
        F = json.load(open(os.path.join(LG.fit_dir(lab), "loc_meta.json")))
        c = c.assign(seed=seed, gppi_contrast=con, N=len(F["subjects"]))
        cl.append(c)
        # within-network summaries: per participant, mean effect over the network minus the seed voxels
        cope = nib.load(os.path.join(LG.fit_dir(lab), "cope.nii.gz")).get_fdata(dtype=np.float32)
        for i, x in enumerate(F["subjects"]):
            sm = seed_masks(x)[seed]
            for n in NETS:
                mm = netm[n] & ~sm
                nw.append(dict(seed=seed, gppi_contrast=con, subject=x, network=n, mean_effect=float(cope[..., i][mm].mean())))
    C = pd.concat(cl, ignore_index=True)
    C.to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_gppi_clusters.csv"), index=False)
    W = pd.DataFrame(nw); rows = []
    for (seed, con, n), d in W.groupby(["seed", "gppi_contrast", "network"], sort=False):
        y = d.mean_effect.values; t, p = stats.ttest_1samp(y, 0)
        rows.append(dict(seed=seed, gppi_contrast=con, network=n, N=len(y), mean=float(y.mean()),
                         dz=float(y.mean() / y.std(ddof=1)), t=float(t), p=float(p)))
    R = pd.DataFrame(rows)
    R["p_fdr50"] = np.nan
    for con in GCON:
        k = R.gppi_contrast == con
        R.loc[k, "p_fdr50"] = LC.fdr(R.loc[k, "p"].values)
    R.to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_gppi_networks.csv"), index=False)
    # first-level diagnostics (deconvolution), all participants and seeds
    q = []
    for x in gsets()[("mask", "C4")]:
        js = json.load(open(os.path.join(GPD, x, f"{x}_gppi.json")))
        for seed, v in js["seeds"].items():
            for r in v["runs"]:
                q.append(dict(subject=x, seed=seed, dof=v["dof"], **r))
    Q = pd.DataFrame(q)
    Q.groupby("seed").agg(n_runs=("run", "size"), r_reconv_median=("r_reconv", "median"), r_reconv_min=("r_reconv", "min"),
                          lambda_ratio_median=("lambda_ratio", "median"), n_lambda_at_grid_edge=("edge", "sum"),
                          psych_r_min=("psych_r_min", "min"), dof_median=("dof", "median")).reset_index() \
        .to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_gppi_firstlevel_qc.csv"), index=False)
    print("gPPI tables written", flush=True)


def surviving():
    C = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_gppi_clusters.csv"))
    return C[C.p_fwe_two_sided < 0.05].reset_index(drop=True)


def rvt_seeds():
    return sorted(surviving().seed.unique())


def rvt_compare():
    """Pre-logged physiology rule: for each surviving cluster, the per-participant mean gPPI
    effect inside the cluster, base vs +RVT, in participants with RVT on every used run."""
    s, _, _ = L.participant_sets()
    rows = []
    for _, c in surviving().iterrows():
        lab = fit_label(c.seed, c.gppi_contrast); od = LG.fit_dir(lab)
        mask, _, _ = cluster_mask(lab, c.sign, c.peak_mni, int(c.size_vox))
        subs = [x for x in json.load(open(os.path.join(od, "loc_meta.json")))["subjects"] if x in s["rvt"]]
        b = np.array([nib.load(gname(x, c.seed, c.gppi_contrast, "effect")).get_fdata(dtype=np.float32)[mask].mean() for x in subs])
        r = np.array([nib.load(gname(x, c.seed, c.gppi_contrast, "effect", "rvt")).get_fdata(dtype=np.float32)[mask].mean() for x in subs])
        tb, pb = stats.ttest_1samp(b, 0); tr, pr = stats.ttest_1samp(r, 0)
        rows.append(dict(seed=c.seed, gppi_contrast=c.gppi_contrast, sign=c.sign, size_vox=int(c.size_vox),
                         peak_mni=c.peak_mni, N_full=int(c.N), N_rvt=len(subs), base_mean=float(b.mean()),
                         base_dz=float(b.mean() / b.std(ddof=1)), base_p=float(pb), rvt_mean=float(r.mean()),
                         rvt_dz=float(r.mean() / r.std(ddof=1)), rvt_p=float(pr),
                         ratio_rvt_over_base=float(r.mean() / b.mean()), r_base_rvt=float(np.corrcoef(b, r)[0, 1])))
    pd.DataFrame(rows).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_gppi_rvt.csv"), index=False)
    print(pd.DataFrame(rows).round(3).to_string())


def cluster_mask(lab, sign, peak_mni, size):
    from scipy import ndimage
    zimg = nib.load(os.path.join(LG.fit_dir(lab), "stats", "zstat1.nii.gz")); z = zimg.get_fdata()
    sg = 1 if sign == "+" else -1
    L_, _ = ndimage.label(sg * z > 3.09, structure=np.ones((3, 3, 3)))
    pk = json.loads(peak_mni) if isinstance(peak_mni, str) else peak_mni
    ijk = np.round(nib.affines.apply_affine(np.linalg.inv(zimg.affine), pk)).astype(int)
    m = L_ == L_[tuple(ijk)]
    assert L_[tuple(ijk)] > 0 and m.sum() == size, (lab, m.sum(), size)
    return m, z, zimg.affine


def short_label(txt):
    t = str(txt).split(";")[0].split(",")[0]
    t = t.split("% ", 1)[-1] if "% " in t else t
    t = t.replace("Inferior Frontal Gyrus", "IFG").replace("Middle Frontal Gyrus", "MFG").replace(" Cortex", "") \
        .replace("Lateral Occipital", "Lat. occipital").replace("Supramarginal Gyrus", "Supramarginal") \
        .replace("posterior division", "post.").replace("superior division", "sup.").replace("Temporal Fusiform", "Temporal fusiform")
    return t


def figure():
    """Maps with a cluster that also survives x20 (p_FWE two-sided < .0025). One axial slice per
    surviving cluster of those maps, at the cluster's peak; only surviving clusters are coloured."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from nilearn import datasets, image
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
    S = surviving()
    keep = S[S.p_fwe_two_sided < 0.05 / 20][["seed", "gppi_contrast"]].drop_duplicates()
    order = [("froiAInsL", "C1"), ("froiAInsR", "C1"), ("maskAInsR", "C1"), ("froidACC", "C4"), ("froiPIns", "C4"),
             ("maskPIns", "C4")]
    maps = [k for k in order if ((keep.seed == k[0]) & (keep.gppi_contrast == k[1])).any()]
    assert len(maps) == len(keep), (maps, keep)
    cname = {"C1": "C1 Intero − Extero", "C4": "C4 Match − Extero"}
    ref = nib.load(g.GROUP_MASK)
    bg = image.resample_to_img(datasets.load_mni152_template(resolution=2), ref, interpolation="continuous",
                               force_resample=True, copy_header=True).get_fdata()
    gm = ref.get_fdata() > 0
    ncol = max(len(S[(S.seed == a) & (S.gppi_contrast == b)]) for a, b in maps)
    norm = Normalize(-5.5, 5.5)
    xi = np.flatnonzero(gm.any(axis=(1, 2))); yi = np.flatnonzero(gm.any(axis=(0, 2)))
    x0, x1, y0, y1 = xi[0] - 1, xi[-1] + 2, yi[0] - 1, yi[-1] + 2
    nr = len(maps)
    fig = plt.figure(figsize=(2.45 * ncol + 0.4, 3.15 * nr + 1.0), facecolor=SURF)
    gs = fig.add_gridspec(2 * nr + 1, ncol, height_ratios=[0.16, 1] * nr + [0.06], hspace=0.28, wspace=0.06,
                          left=0.02, right=0.98, top=0.935, bottom=0.04)
    vmin, vmax = np.percentile(bg[gm], 1), np.percentile(bg[gm], 99.5)
    for i, (seed, con) in enumerate(maps):
        lab = fit_label(seed, con)
        cl = S[(S.seed == seed) & (S.gppi_contrast == con)].sort_values("p_fwe_two_sided")
        m = np.zeros(ref.shape, bool)
        for _, c in cl.iterrows():
            mm, z, aff = cluster_mask(lab, c.sign, c.peak_mni, int(c.size_vox)); m |= mm
        zz = np.where(m, z, np.nan)
        kind = "fROI" if seed.startswith("froi") else "mask"
        top = fig.add_subplot(gs[2 * i, :]); top.axis("off")
        top.text(0, 0.0, f"{kind} {seed[4:]} seed · gPPI {cname[con]} · N {int(cl.N.iloc[0])}",
                 fontsize=10.5, color=INK, transform=top.transAxes, va="bottom")
        for j in range(ncol):
            ax = fig.add_subplot(gs[2 * i + 1, j]); ax.axis("off"); ax.set_facecolor(SURF)
            if j >= len(cl):
                continue
            c = cl.iloc[j]
            pk = json.loads(c.peak_mni) if isinstance(c.peak_mni, str) else c.peak_mni
            k = int(round(nib.affines.apply_affine(np.linalg.inv(aff), pk)[2]))
            sl = np.where(gm[:, :, k], bg[:, :, k], np.nan)[x0:x1, y0:y1].T
            ax.imshow(sl, cmap="gray", origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax)
            ax.imshow(zz[x0:x1, y0:y1, k].T, cmap=cmap, norm=norm, origin="lower", interpolation="nearest")
            star = "" if c.p_fwe_two_sided < 0.05 / 20 else " (×20 n.s.)"
            ax.set_title(f"{short_label(c.labels)}, z = {int(round(pk[2]))}\n{int(c.size_vox)} vox, p {c.p_fwe_two_sided:.1g}{star}",
                         fontsize=7.8, color=INK2, pad=3)
    cax = fig.add_subplot(gs[-1, max(0, ncol // 2 - 1): ncol // 2 + 2])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label("FLAME z; red = stronger seed coupling in the first condition", color=INK2, fontsize=8.5)
    cb.outline.set_visible(False)
    fig.suptitle("Localizer gPPI (exploratory): maps with a cluster surviving ×20 across the 20 maps.\n"
                 "One axial slice per surviving cluster (whole-brain cluster FWE, z > 3.09, two-sided), at its peak; "
                 "image left = participant left.", x=0.02, ha="left", fontsize=10.5, color=INK)
    fig.savefig(os.path.join(LC.FIG, "loc_C2_gppi_maps.png"), dpi=150, facecolor=SURF)
    print("wrote Figures/localizer/loc_C2_gppi_maps.png")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "subject":
        subject(sys.argv[2])
    elif cmd == "subjectrvt":
        subject(sys.argv[2], "rvt", set(rvt_seeds()))
    elif cmd == "check":
        check(sys.argv[2])
    else:
        {"flames": flames, "tables": tables, "rvtcompare": rvt_compare, "figure": figure}[cmd]()
