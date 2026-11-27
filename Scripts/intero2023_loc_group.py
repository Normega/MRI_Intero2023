#!/usr/bin/env python
"""
intero2023_loc_group.py -- registered localizer group tests, stage A (OSF p5ztc):
L1, L2, L3, the L1 breathing-robustness co-report, the registered sensitivity
analyses for L1-L3, and the descriptive whole-brain maps for C1-C6.
Choices: Analysis/localizer_glm_implementation_notes.md (committed before any statistic).

  python Scripts/intero2023_loc_group.py flames    # every FLAME 1 fit, in parallel (FSL flameo)
  python Scripts/intero2023_loc_group.py tables    # mask tests, permutation checks, robustness,
                                                   # sensitivity, cluster tables -> Analysis/results/LOC_*

Reuses intero2023_group (sh, grf_p, reml_meta, ols_t, peak_mni, N_PERM, SEED, ROI) and
intero2023_explore_awareness (jzs_bf10, bf_plus0), unchanged.
Maps: $OUT/group/LOC_<fit>/ (outside the repo).
"""
import os, sys, json, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats, integrate, ndimage
from scipy.special import gammaln
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
import intero2023_group as G
from intero2023_explore_awareness import jzs_bf10, bf_plus0
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
GD = os.path.join(g.OUT, "group")
N_JOBS = int(os.environ.get("GROUP_NJOBS", "12"))
R_BF = 0.707
TESTS = {   # registered confirmatory tests (plan §9, p5ztc Inference criteria)
    "L1": dict(con="C1", set="intero", rois=["caren_AInsL", "caren_AInsR", "caren_dACC"], tail="two", k=3, maia=False),
    "L2": dict(con="C5", set="intero", rois=["caren_PIns", "caren_COper", "caren_SM"], tail="neg", k=3, maia=False),
    "L3": dict(con="C1", set="L3", rois=["ieat_ACCsphere"], tail="two", k=1, maia=True),
}
DESCRIPTIVE = {"C1": "intero", "C2": "intero", "C5": "intero", "intero": "intero",
               "C4": "all", "C6": "all", "extero": "all", "match": "all"}


# --------------------------------------------------------------------------- inputs
def covariates(u):
    """ΔRVT and Extero entrainment per participant, from the step-1 run-level rows,
    restricted to the runs that enter the imaging model; protocol effect-coded."""
    b = pd.read_csv(os.path.join(L.QCDIR, "breathing_by_condition.csv"))
    b = b[b.level == "run"].copy(); b["run"] = b.run.astype(int)
    b = b.merge(u[["subject", "run"]], on=["subject", "run"])
    ps = b.groupby(["subject", "condition"])[["rvt_z_mean", "plv"]].mean().reset_index()
    rvt = ps.pivot(index="subject", columns="condition", values="rvt_z_mean")
    plv = ps.pivot(index="subject", columns="condition", values="plv")
    null = ps[ps.condition.isin(["rest", "self", "intero"])].groupby("subject").plv.mean()
    c = pd.DataFrame({"drvt": rvt["intero"] - rvt["extero"], "entrain": plv["extero"] - null})
    pv = u.groupby("subject").protocol_version.first()
    c = c.reindex(sorted(u.subject.unique()))
    c["protocol"] = np.where(pv.reindex(c.index) == "v1", 0.5, -0.5)
    return c


def fit_spec(label):
    """label -> (test, con, var, subjects, covariate names). Sensitivity labels are '<T>_<s>'."""
    s, maia, u = L.participant_sets()
    cov = covariates(u)
    if label in DESCRIPTIVE:
        return dict(test=None, con=label, var="base", subs=s[DESCRIPTIVE[label]], covs=[], maia=False)
    T, _, sens = label.partition("_")
    t = TESTS[T]; subs = list(s[t["set"]]); var = "base"; covs = []
    if sens in ("hp128", "hp320"):
        var = sens
    elif sens == "nov1":
        subs = [x for x in subs if x not in s["v1"]]
    elif sens in ("covprotocol", "covdrvt", "coventrain"):
        covs = [sens[3:]]
        subs = [x for x in subs if pd.notna(cov.loc[x, covs[0]])]
    elif sens in ("co2base", "co2"):
        subs = [x for x in subs if x in s["co2"]]; var = "co2" if sens == "co2" else "base"
    elif sens in ("rvtbase", "rvt"):
        subs = [x for x in subs if x in s["rvt_tworun"]]; var = "rvt" if sens == "rvt" else "base"
    elif sens:
        raise KeyError(label)
    return dict(test=T, con=t["con"], var=var, subs=subs, covs=covs, maia=t["maia"])


SENS = ["hp128", "hp320", "nov1", "covprotocol", "covdrvt", "coventrain", "co2base", "co2", "rvtbase", "rvt"]
ALL_FITS = ["L1", "L2", "L3"] + [f"{T}_{x}" for T in TESTS for x in SENS] + \
    [d for d in DESCRIPTIVE if d not in ("C1", "C5")]


def fit_dir(label):
    same = {"C1": "L1", "C5": "L2"}   # the C1 and C5 whole-brain maps are the L1 / L2 fits
    return os.path.join(GD, "LOC_" + same.get(label, label))


# --------------------------------------------------------------------------- FLAME
def design_matrix(sp):
    s, maia, u = L.participant_sets()
    cov = covariates(u)
    cols = {"intercept": np.ones(len(sp["subs"]))}
    if sp["maia"]:
        x = maia.loc[sp["subs"], "maia_ieat_total"].astype(float).values
        cols["maia"] = x - x.mean()
    for c in sp["covs"]:
        x = cov.loc[sp["subs"], c].astype(float).values
        cols[c] = x - x.mean()
    X = pd.DataFrame(cols, index=sp["subs"])
    test_col = "maia" if sp["maia"] else "intercept"
    return X, test_col


def dofs(subs, var):
    return [json.load(open(os.path.join(L.fl_dir(x), f"{x}_task-localizer_firstlevel.json")))
            ["fits"][f"joint:sm6:{var}"]["dof"] for x in subs]


def run_flame(label):
    sp = fit_spec(label); od = fit_dir(label)
    os.makedirs(od, exist_ok=True)
    X, test_col = design_matrix(sp)
    meta = dict(label=label, con=sp["con"], var=sp["var"], subjects=sp["subs"], design=list(X.columns))
    mf = os.path.join(od, "loc_meta.json")
    if os.path.exists(mf) and os.path.exists(os.path.join(od, "stats", "zstat1.nii.gz")):
        old = json.load(open(mf))
        if all(old.get(k) == meta[k] for k in meta):
            return label, "cached"
    eff = [L.fl_name(x, "joint", 6, sp["var"], sp["con"], "effect") for x in sp["subs"]]
    var = [L.fl_name(x, "joint", 6, sp["var"], sp["con"], "variance") for x in sp["subs"]]
    ref = nib.load(g.GROUP_MASK)
    stack = lambda files: nib.Nifti1Image(np.stack([nib.load(f).get_fdata(dtype=np.float32) for f in files], -1),
                                          ref.affine)
    stack(eff).to_filename(f"{od}/cope.nii.gz"); stack(var).to_filename(f"{od}/varcope.nii.gz")
    m = ref.get_fdata() > 0
    nib.Nifti1Image(np.stack([m * d for d in dofs(sp["subs"], sp["var"])], -1).astype(np.float32),
                    ref.affine).to_filename(f"{od}/dof.nii.gz")
    n, p = X.shape
    open(f"{od}/design.mat", "w").write(f"/NumWaves {p}\n/NumPoints {n}\n/Matrix\n" +
                                        "".join(" ".join(f"{v:.6f}" for v in row) + "\n" for row in X.values))
    c = " ".join("1" if col == test_col else "0" for col in X.columns)
    open(f"{od}/design.con", "w").write(f"/NumWaves {p}\n/NumContrasts 1\n/Matrix\n{c}\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    sm = G.sh(f"smoothest -d {n - p} -r {od}/stats/res4d.nii.gz -m {g.GROUP_MASK}")
    meta["smoothest"] = {l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()}
    meta["dofs"] = dofs(sp["subs"], sp["var"])
    os.remove(f"{od}/stats/res4d.nii.gz")          # large; smoothness is recorded
    json.dump(meta, open(mf, "w"), indent=1)
    return label, "fitted"


def flames():
    with ThreadPoolExecutor(int(os.environ.get("FLAME_JOBS", "10"))) as ex:
        for label, st in ex.map(run_flame, ALL_FITS):
            print(f"{label}: {st}", flush=True)


# --------------------------------------------------------------------------- statistics
def jzs_reg_bf10(r2, n, r=R_BF, p=1):
    """JZS (Zellner-Siow) regression BF10 for p covariates from R^2 (Liang et al. 2008;
    Rouder & Morey 2012), prior g ~ InvGamma(1/2, r^2 n / 2), as BayesFactor linearReg.R2stat."""
    a, b = 0.5, r ** 2 * n / 2
    def f(gg):
        lg = ((n - 1 - p) / 2) * np.log1p(gg) - ((n - 1) / 2) * np.log1p(gg * (1 - r2)) \
            + a * np.log(b) - gammaln(a) - (a + 1) * np.log(gg) - b / gg
        return np.exp(lg)
    return float(integrate.quad(f, 0, np.inf, limit=200)[0])


def load_fit(label):
    od = fit_dir(label); meta = json.load(open(os.path.join(od, "loc_meta.json")))
    st = os.path.join(od, "stats")
    return dict(od=od, meta=meta, z=nib.load(f"{st}/zstat1.nii.gz").get_fdata(),
                cope=nib.load(f"{od}/cope.nii.gz").get_fdata(dtype=np.float32),
                varcope=nib.load(f"{od}/varcope.nii.gz").get_fdata(dtype=np.float32),
                revar=nib.load(f"{st}/mean_random_effects_var1.nii.gz").get_fdata(),
                affine=nib.load(f"{st}/zstat1.nii.gz").affine)


def flame_roi(label, T):
    """FLAME + GRF mask tests for test T on fit `label`. Returns (roi rows, subject rows)."""
    t = TESTS[T]; F = load_fit(label); sp = fit_spec(label)
    subs = F["meta"]["subjects"]; assert subs == sp["subs"]
    X, test_col = design_matrix(sp)
    smv = F["meta"]["smoothest"]; z = F["z"]
    tails = 2 if t["tail"] == "two" else 1
    rows, srows = [], []
    for roi in t["rois"]:
        rm = nib.load(G.ROI(roi)).get_fdata() > 0
        zr = z[rm]
        flame_ok = bool(np.isfinite(zr).all() and (zr != 0).mean() > 0.99)
        resels = rm.sum() / smv["RESELS"]
        if t["tail"] == "two":
            zz = np.where(rm, np.abs(z), -np.inf)
        elif t["tail"] == "neg":
            zz = np.where(rm, -z, -np.inf)
        else:
            zz = np.where(rm, z, -np.inf)
        ijk = np.unravel_index(np.argmax(zz), z.shape); zt = float(zz[ijk])
        p = min(1.0, tails * G.grf_p(zt, resels)) if zt > 1 else 1.0
        zthr = float(G.sh(f"ptoz {0.05 / (t['k'] * tails)} -g {resels}").strip())
        nsig = int((np.abs(zr) > zthr).sum() if t["tail"] == "two" else ((zr if t["tail"] == "pos" else -zr) > zthr).sum())
        y = F["cope"][rm].mean(axis=0); v = F["varcope"][rm].mean(axis=0); re = float(F["revar"][rm].mean())
        w = 1 / (v + re); share = w / w.sum()
        rec = dict(test=T, label=f"LOC_{label}", engine="FLAME1+GRF", contrast=sp["con"], variant=sp["var"],
                   covariates=";".join(sp["covs"]), roi=roi, N=len(subs), direction={"two": "two-sided", "neg": "negative", "pos": "positive"}[t["tail"]],
                   peak_z=round(float(z[ijk]), 3), peak_mni=[round(float(a), 1) for a in nib.affines.apply_affine(F["affine"], ijk)],
                   resels=round(float(resels), 2), p_fwe=p, p_fwe_bonf=min(1.0, t["k"] * p), bonferroni_k=t["k"],
                   z_thresh_bonf=round(zthr, 3), n_vox_sig=nsig, flame_ok=flame_ok,
                   roi_mean=float(y.mean()), max_weight_share=float(share.max()))
        if t["maia"]:
            x = X["maia"].values
            r_, p_r = stats.pearsonr(x, y)
            ols = stats.linregress(x, y)
            rec.update(maia_r=float(r_), maia_slope=float(ols.slope), maia_slope_t=float(ols.slope / ols.stderr),
                       maia_p_two=float(p_r), bf10=jzs_reg_bf10(r_ ** 2, len(y)), bf10_r354=jzs_reg_bf10(r_ ** 2, len(y), r=0.354))
        else:
            tt = stats.ttest_1samp(y, 0)
            meta_ = G.reml_meta(y, v)
            rec.update(roi_dz=float(y.mean() / y.std(ddof=1)), roi_t=float(tt.statistic), roi_p_two=float(tt.pvalue),
                       bf10=jzs_bf10(float(tt.statistic), len(y), R_BF),
                       reml_mu=meta_["mu"], reml_z=meta_["z"], reml_tau2=meta_["tau2"])
            if t["tail"] == "neg":
                rec["bf_dir"] = bf_plus0(-y, R_BF)[1]       # BF-0 (effect < 0)
        rows.append(rec)
        for s_, yy, vv, sh_ in zip(subs, y, v, share):
            srows.append(dict(label=f"LOC_{label}", roi=roi, subject=s_, mask_mean=float(yy), mask_mean_varcope=float(vv),
                              flame_weight_share=float(sh_), **{c: float(X.loc[s_, c]) for c in X.columns if c != "intercept"}))
    return rows, srows


def perm_roi(label, T):
    """Permutation voxel max-t check (10,000 sign-flips / MAIA permutations, seed 20261126)."""
    from nilearn.glm.second_level import non_parametric_inference
    from nilearn.image import math_img
    t = TESTS[T]; sp = fit_spec(label); X, test_col = design_matrix(sp)
    imgs = [nib.load(L.fl_name(x, "joint", 6, sp["var"], sp["con"], "effect")) for x in sp["subs"]]
    dm = X.reset_index(drop=True)
    two = t["tail"] == "two"
    use = imgs if t["tail"] != "neg" else [math_img("-i", i=i) for i in imgs]
    rows = []
    for roi in t["rois"]:
        rm = nib.load(G.ROI(roi)).get_fdata() > 0
        lp = non_parametric_inference(use, design_matrix=dm, second_level_contrast=test_col, mask=G.ROI(roi),
                                      n_perm=G.N_PERM, two_sided_test=two, random_state=G.SEED, n_jobs=N_JOBS)
        tmap = G.ols_t(use, dm, test_col, G.ROI(roi))
        lpv = lp.get_fdata()[rm]; pmin = float(10 ** -np.nanmax(lpv))
        tv = tmap.get_fdata()
        tt = np.where(rm, np.abs(tv) if two else tv, -np.inf)
        ijk = np.unravel_index(np.argmax(tt), tv.shape)
        t_signed = float(tv[ijk]) * (-1 if t["tail"] == "neg" else 1)
        rows.append(dict(test=T, label=f"LOC_{label}", engine="permutation max-t (OLS)", contrast=sp["con"], variant=sp["var"],
                         covariates="", roi=roi, N=len(sp["subs"]),
                         direction={"two": "two-sided", "neg": "negative", "pos": "positive"}[t["tail"]],
                         peak_t=round(t_signed, 3), peak_mni=[round(float(a), 1) for a in nib.affines.apply_affine(tmap.affine, ijk)],
                         p_fwe=pmin, p_fwe_bonf=min(1.0, t["k"] * pmin), bonferroni_k=t["k"],
                         n_vox_sig=int((lpv > -np.log10(0.05 / t["k"])).sum())))
        g.save_float32(tmap, os.path.join(fit_dir(label), f"perm_{roi}_t.nii.gz"))
        g.save_float32(lp, os.path.join(fit_dir(label), f"perm_{roi}_logpFWE.nii.gz"))
    return rows


def atlas_labels(mni):
    out = []
    for atlas in ("Harvard-Oxford Cortical Structural Atlas", "Harvard-Oxford Subcortical Structural Atlas",
                  "Cerebellar Atlas in MNI152 space after normalization with FNIRT"):
        try:
            txt = G.sh(f'atlasquery -a "{atlas}" -c {mni[0]},{mni[1]},{mni[2]}').strip()
        except Exception as ex:
            txt = f"<error {ex}>"
        lab = txt.split("<br>", 1)[-1] if "<br>" in txt else txt
        if lab and "No label found" not in lab:
            out.append(lab.strip())
    return "; ".join(out)


def clusters(label, contrast_name):
    """Whole-brain cluster FWE at z > 3.09, both tails, GRF cluster p x 2 (as H3)."""
    od = fit_dir(label); st = os.path.join(od, "stats"); smv = json.load(open(f"{od}/loc_meta.json"))["smoothest"]
    out = []
    for sgn, zf in ((+1, f"{st}/zstat1.nii.gz"), (-1, f"{st}/zneg.nii.gz")):
        if sgn < 0:
            G.sh(f"fslmaths {st}/zstat1.nii.gz -mul -1 {zf}")
        txt = G.sh(f"cluster -i {zf} -t 3.09 -p 1 -d {smv['DLH']} --volume={int(smv['VOLUME'])} --mm")
        lines = [l.split("\t") for l in txt.strip().splitlines()]
        for l in lines[1:]:
            d = dict(zip(lines[0], l))
            mni = [float(d["MAX X (mm)"]), float(d["MAX Y (mm)"]), float(d["MAX Z (mm)"])]
            out.append(dict(contrast=contrast_name, sign="+" if sgn > 0 else "-", size_vox=int(d["Voxels"]),
                            peak_z=sgn * float(d["MAX"]), peak_mni=mni, p_cluster=float(d["P"]),
                            p_fwe_two_sided=min(1.0, 2 * float(d["P"]))))
    c = pd.DataFrame(out)
    if len(c):
        c = c.sort_values("size_vox", ascending=False)
        c["labels"] = [atlas_labels(m) if p < 0.05 or i < 10 else "" for i, (m, p) in
                       enumerate(zip(c.peak_mni, c.p_fwe_two_sided))]
    return c


def conjunction():
    """Descriptive C3 group map: conjunction-null min of the C1 and C2 FLAME z maps (same 48)."""
    z1 = load_fit("L1")["z"]; f2 = load_fit("C2"); z2 = f2["z"]
    assert load_fit("L1")["meta"]["subjects"] == f2["meta"]["subjects"]
    m = nib.load(g.GROUP_MASK).get_fdata() > 0
    out = []
    for sgn in (+1, -1):
        zc = np.where(m, np.minimum(sgn * z1, sgn * z2), 0)
        lab, n = ndimage.label(zc > 3.09, structure=np.ones((3, 3, 3)))
        for k in range(1, n + 1):
            vox = lab == k
            ijk = np.unravel_index(np.argmax(np.where(vox, zc, -np.inf)), zc.shape)
            mni = [round(float(a), 1) for a in nib.affines.apply_affine(f2["affine"], ijk)]
            out.append(dict(contrast="C3 conjunction min(z_C1, z_C2)", sign="+" if sgn > 0 else "-",
                            size_vox=int(vox.sum()), peak_z=round(sgn * float(zc[ijk]), 3), peak_mni=mni,
                            p_cluster=np.nan, p_fwe_two_sided=np.nan))
        g.save_float32(nib.Nifti1Image(zc.astype(np.float32), f2["affine"]),
                       os.path.join(GD, f"LOC_C3conj_{'pos' if sgn > 0 else 'neg'}_minz.nii.gz"))
    c = pd.DataFrame(out)
    if len(c):
        c = c.sort_values("size_vox", ascending=False)
        c["labels"] = [atlas_labels(m) if s >= 20 else "" for m, s in zip(c.peak_mni, c.size_vox)]
    return c


def write(label, roi=None, wb=None, subj=None):
    os.makedirs(RD, exist_ok=True)
    for kind, d in (("roi", roi), ("wholebrain", wb), ("subjects", subj)):
        if d is not None:
            pd.DataFrame(d).to_csv(os.path.join(RD, f"LOC_{label}_{kind}.csv"), index=False)


def tables():
    s, maia, u = L.participant_sets()
    sens_rows = []
    # --- confirmatory tests, with the permutation check and whole-brain maps
    for T in TESTS:
        r, sr = flame_roi(T, T)
        pr = perm_roi(T, T)
        wb = clusters(T, f"{TESTS[T]['con']}" + (" x MAIA" if TESTS[T]["maia"] else ""))
        write(T, r + pr, wb, sr)
        sens_rows += [dict(x, sensitivity="primary") for x in r]
        print(f"{T}: done", flush=True)
    # --- sensitivity analyses
    for T in TESTS:
        for x in SENS:
            r, sr = flame_roi(f"{T}_{x}", T)
            write(f"{T}_{x}", r, None, sr)
            sens_rows += [dict(y, sensitivity=x) for y in r]
    S = pd.DataFrame(sens_rows)
    S.to_csv(os.path.join(RD, "LOC_sensitivity.csv"), index=False)
    # --- L1 breathing-robustness rule (same participants, base vs +RVT)
    rb = []
    b = pd.read_csv(os.path.join(RD, "LOC_L1_rvtbase_subjects.csv"))
    rv = pd.read_csv(os.path.join(RD, "LOC_L1_rvt_subjects.csv"))
    for roi in TESTS["L1"]["rois"]:
        y0 = b[b.roi == roi].set_index("subject").mask_mean
        y1 = rv[rv.roi == roi].set_index("subject").mask_mean.reindex(y0.index)
        assert y1.notna().all()
        m0, m1 = float(y0.mean()), float(y1.mean())
        ratio = m1 / m0 if m0 != 0 else np.nan
        robust = bool(np.sign(m0) == np.sign(m1) and ratio >= 0.5)
        f0 = S[(S.sensitivity == "rvtbase") & (S.test == "L1") & (S.roi == roi)].iloc[0]
        f1 = S[(S.sensitivity == "rvt") & (S.test == "L1") & (S.roi == roi)].iloc[0]
        rb.append(dict(roi=roi, N=len(y0), base_mask_mean=m0, rvt_mask_mean=m1, ratio_rvt_over_base=ratio,
                       same_sign=bool(np.sign(m0) == np.sign(m1)), verdict="breathing-robust" if robust else "breathing-dependent",
                       base_dz=float(y0.mean() / y0.std(ddof=1)), rvt_dz=float(y1.mean() / y1.std(ddof=1)),
                       r_base_rvt=float(np.corrcoef(y0, y1)[0, 1]),
                       base_peak_z=f0.peak_z, base_p_fwe_bonf=f0.p_fwe_bonf, rvt_peak_z=f1.peak_z, rvt_p_fwe_bonf=f1.p_fwe_bonf,
                       subjects=";".join(y0.index)))
    pd.DataFrame(rb).to_csv(os.path.join(RD, "LOC_L1_breathing_robustness.csv"), index=False)
    # --- descriptive whole-brain maps
    for con in DESCRIPTIVE:
        lab = {"C1": "L1", "C5": "L2"}.get(con, con)
        F = load_fit(lab)
        wb = clusters(lab, con) if con not in ("C1", "C5") else pd.read_csv(os.path.join(RD, f"LOC_{lab}_wholebrain.csv"))
        if con in ("C1", "C5"):
            wb["contrast"] = con
        write(con, None, wb, pd.DataFrame(dict(subject=F["meta"]["subjects"])))
    write("C3conj", None, conjunction(), pd.DataFrame(dict(subject=load_fit("L1")["meta"]["subjects"])))
    print("tables written", flush=True)


if __name__ == "__main__":
    {"flames": flames, "tables": tables}[sys.argv[1]]()
