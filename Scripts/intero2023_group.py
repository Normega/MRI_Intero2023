#!/usr/bin/env python
"""
intero2023_group.py -- group-level tests for H1, H2, H4, H6 (nilearn permutation)
and H3 (FSL FLAME 1 + GRF), prereg v0.8 §8; choices in Analysis/glm_implementation_notes.md.

Usage:
  python Scripts/intero2023_group.py <analysis> [<variant>] [<subset>]
    analysis : H1 | H2 | H4 | H6 | H3 | H3perrun | S7ols
    variant  : base (default) | rvt | co2 | both
    subset   : all (default) | S1 | S2 | S2b | S4 | S6 | S7 | both
Writes maps to $OUT/group/<label>/ and tables to Analysis/results/<label>_*.csv.
"""
import os, sys, json, subprocess, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats, optimize
from nilearn.glm.second_level import non_parametric_inference
from nilearn.maskers import NiftiMasker
from nilearn.image import math_img, new_img_like
import intero2023_glm_common as g
q = g.q
warnings.simplefilter("ignore")
N_PERM, SEED, N_JOBS = 10000, 20261126, int(os.environ.get("GROUP_NJOBS", "12"))
ATLAS = os.environ.get("GROUP_ATLAS", "caren")        # caren | ho (atlas sensitivity, never inferential)
WB = os.environ.get("GROUP_WB", "1") == "1"           # whole-brain test (primary analyses)
COV = os.environ.get("GROUP_COV", "")                 # "" | protocol (S5)
FSL = os.environ.get("FSLDIR", "/home/normega/fsl")
ROI = lambda n: os.path.join(g.REPO, "ROIs", f"{n}.nii")
SETS = dict(H1=["caren_PIns", "caren_COper", "caren_SM"],
            H2=["caren_AInsL", "caren_AInsR", "caren_dACC"],
            H3=["caren_AInsL", "caren_AInsR", "caren_dACC"],
            H4=["caren_AInsL", "caren_AInsR", "caren_dACC"],
            H6=["ieat_ACCsphere"])
DESCRIPTIVE = dict(H6=["caren_dACC"])
SETS["S3"] = SETS["H2"]
if ATLAS == "ho":
    for k in ("H1", "H2", "H3", "H4", "S3"):
        SETS[k] = [r.replace("caren_", "ho_") for r in SETS[k]]
CON = dict(H1="C1", H2="C3", H4="C5", H6="C1", H3="C2joint", H3perrun="C2", S7ols="C2joint", S3="C3")
SIGN = dict(H1=-1, H2=+1, H3=+1, H4=+1)          # predicted direction (H6 two-sided)


def fl(sub, var, con, kind):
    return os.path.join(g.OUT, "firstlevel", sub, f"{sub}_var-{var}_con-{con}_{kind}.nii.gz")


def subject_sets():
    man = g.run_list()
    for c in ("has_rvt", "has_petco2"):
        man[c] = man[c].astype(str).str.upper() == "TRUE"
    gb = man.groupby("subject")
    allsub = sorted(man.subject.unique())
    s = dict(all=allsub)
    s["S1"] = sorted(gb.has_petco2.all().pipe(lambda x: x[x]).index)
    s["S2"] = sorted(gb.rvt_source.apply(lambda x: (x == "biopac").all()).pipe(lambda x: x[x]).index)
    s["S2b"] = sorted(gb.apply(lambda d: (d.has_rvt & d.rvt_source.isin(["biopac", "pmu"])).all())
                      .pipe(lambda x: x[x]).index)
    s["both"] = sorted(set(s["S1"]) & set(s["S2"]))
    s["S6"] = sorted(gb.size().pipe(lambda x: x[x == 2]).index)
    # S4: verified pacer adherence on at least one retained run (reproduces the registered 26)
    adh = pd.read_csv(os.path.join(g.REPO, "Data", "resp_epochs.csv")).merge(
        man[["subject", "bids_entity"]], on=["subject", "bids_entity"])
    s["S4"] = sorted(adh[adh.adherence_verified.astype(str).str.upper() == "TRUE"].subject.unique())
    h3 = pd.read_csv(os.path.join(g.REPO, "Data", "h3_sets.csv"))
    s["H3"] = sorted(h3[h3.in_H3].subject); s["S7"] = sorted(h3[h3.in_S7].subject)
    # H4: >= 3 high and >= 3 low salience change trials pooled (prereg §4)
    h4 = []
    for sub, d in man.groupby("subject"):
        ev = pd.concat([pd.read_csv(f"{q.FL}/{sub}_{e}_desc-model_events.tsv", sep="\t")
                        for e in d.bids_entity])
        ch = ev[ev.condition.isin(q.CHANGE_ONLY)]
        if (ch.salience == "high").sum() >= 3 and (ch.salience == "low").sum() >= 3:
            h4.append(sub)
    s["H4"] = h4
    maia = pd.read_csv(os.path.join(g.REPO, "Data", "maia_scores.csv")).set_index("subject")
    s["H6"] = [x for x in allsub if x in maia.index and pd.notna(maia.loc[x, "maia_ieat_total"])]
    return s, maia, man


def subjects_for(analysis, variant, subset, sets):
    base = dict(H1=sets["all"], H2=sets["all"], H4=sets["H4"], H6=sets["H6"], H3=sets["H3"],
                H3perrun=sets["H3"], S7ols=sets["S7"], S3=sets["all"])[analysis]
    if analysis in ("H3", "H3perrun") and subset == "S7":
        base = sets["S7"]
    keep = [x for x in base if subset in ("all", "S7") or x in sets[subset]]
    con = CON[analysis]
    keep = [x for x in keep if os.path.exists(fl(x, variant, con, "effect"))]
    return keep


def peak_mni(img, mask_img, sign):
    d = img.get_fdata() * sign; m = nib.load(mask_img).get_fdata() > 0
    d = np.where(m, d, -np.inf); ijk = np.unravel_index(np.argmax(d), d.shape)
    return [round(float(x), 1) for x in nib.affines.apply_affine(img.affine, ijk)], float(d[ijk] * sign)


def roi_mean(imgs, mask):
    mk = NiftiMasker(mask_img=mask).fit()
    return np.array([mk.transform(i).mean() for i in imgs])


def nilearn_tests(analysis, subs, variant, label, maia):
    od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
    con = CON[analysis]
    imgs = [nib.load(fl(s, variant, con, "effect")) for s in subs]
    rows = []
    if analysis == "H6":
        x = maia.loc[subs, "maia_ieat_total"].astype(float).values
        dm = pd.DataFrame({"intercept": 1.0, "maia": x - x.mean()})
        tests = [("maia", +1, True, SETS["H6"] + DESCRIPTIVE["H6"])]
        contrast = "maia"
    elif analysis == "S3":   # cohort check: CO2-available (S1) vs the rest, H2 contrast, two-sided
        s1 = set(subject_sets()[0]["S1"])
        x = np.array([1.0 if s_ in s1 else 0.0 for s_ in subs])
        dm = pd.DataFrame({"intercept": 1.0, "co2group": x - x.mean()})
        tests = [("co2group", +1, True, SETS["S3"])]
        contrast = "co2group"
    else:
        dm = pd.DataFrame({"intercept": np.ones(len(subs))})
        contrast = "intercept"
        sgn = SIGN[analysis]
        tests = [("predicted", sgn, False, SETS[analysis]), ("opposite", -sgn, False, SETS[analysis])]
    if COV == "protocol":
        pv = protocol_of(subs); dm["protocol"] = pv - pv.mean()
    for tag, sgn, two, rois in tests:
        use = imgs if (two or sgn > 0) else [math_img("-i", i=i) for i in imgs]
        for roi in rois:
            descr = roi in DESCRIPTIVE.get(analysis, [])
            k_bonf = 1 if (descr or analysis == "H6") else 3
            rmask = nib.load(ROI(roi)).get_fdata() > 0
            res = non_parametric_inference(use, design_matrix=dm, second_level_contrast=contrast,
                                           mask=ROI(roi), n_perm=N_PERM, two_sided_test=two,
                                           random_state=SEED, n_jobs=N_JOBS)
            tmap, lp = res["t"], res["logp_max_t"]
            lpv = lp.get_fdata()[rmask]
            pmin = float(10 ** -np.nanmax(lpv))
            if two:   # peak of |t|, reported with its sign
                coord, _ = peak_mni(math_img("abs(i)", i=tmap), ROI(roi), +1)
                ijk = tuple(np.round(nib.affines.apply_affine(np.linalg.inv(tmap.affine), coord)).astype(int))
                t_signed = float(tmap.get_fdata()[ijk])
            else:     # tmap is in the tested direction; report the signed original t
                coord, t_tested = peak_mni(tmap, ROI(roi), +1)
                t_signed = sgn * t_tested
            means = roi_mean(imgs, ROI(roi))
            rows.append(dict(analysis=analysis, label=label, variant=variant, test=tag,
                             direction=("two-sided" if two else ("positive" if sgn > 0 else "negative")),
                             roi=roi, N=len(subs), peak_t_signed=round(t_signed, 3), peak_mni=coord,
                             p_fwe_min=pmin, p_fwe_bonf=min(1.0, pmin * k_bonf), bonferroni_k=k_bonf,
                             n_vox_sig_bonf=int((lpv > -np.log10(0.05 / k_bonf)).sum()), descriptive=descr,
                             roi_mean=np.nan if analysis in ("H6", "S3") else float(means.mean()),
                             roi_dz=np.nan if analysis in ("H6", "S3") else float(means.mean() / means.std(ddof=1))))
            tmap.to_filename(os.path.join(od, f"{tag}_{roi}_t.nii.gz"))
            lp.to_filename(os.path.join(od, f"{tag}_{roi}_logpFWE.nii.gz"))
    if not WB:
        return pd.DataFrame(rows), pd.DataFrame()
    # whole brain: cluster-forming p < .001 per tail, cluster-size FWE two-sided
    wb = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast=contrast,
                                  mask=g.GROUP_MASK, n_perm=N_PERM, two_sided_test=True,
                                  threshold=0.002, random_state=SEED, n_jobs=N_JOBS)
    for k in ("t", "logp_max_size"):
        wb[k].to_filename(os.path.join(od, f"wholebrain_{k}.nii.gz"))
    clusters = cluster_table(wb["t"], wb["logp_max_size"], len(subs) - dm.shape[1])
    return pd.DataFrame(rows), clusters


def protocol_of(subs):
    b = pd.read_csv(os.path.join(g.REPO, "Data", "breathing_trials.csv"))
    pv = b.groupby("subject").protocol_version.first()
    return np.array([1.0 if pv[s_] == "v1" else 0.0 for s_ in subs])


def cluster_table(tmap, logp, df):
    from scipy import ndimage
    thr = stats.t.isf(0.001, df)
    t = tmap.get_fdata(); lp = logp.get_fdata(); out = []
    for sgn in (+1, -1):
        lab, n = ndimage.label(sgn * t > thr, structure=ndimage.generate_binary_structure(3, 1))
        for k in range(1, n + 1):
            vox = lab == k
            ijk = np.unravel_index(np.argmax(np.where(vox, sgn * t, -np.inf)), t.shape)
            out.append(dict(sign="+" if sgn > 0 else "-", size_vox=int(vox.sum()),
                            peak_t=round(float(t[ijk]), 3),
                            peak_mni=[round(float(x), 1) for x in nib.affines.apply_affine(tmap.affine, ijk)],
                            p_fwe=float(10 ** -lp[ijk])))
    d = pd.DataFrame(out)
    return d.sort_values("size_vox", ascending=False) if len(d) else d


# ------------------------------------------------------------------ H3 / FLAME
def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                       env=dict(os.environ, FSLDIR=FSL, FSLOUTPUTTYPE="NIFTI_GZ",
                                PATH=f"{FSL}/share/fsl/bin:{FSL}/bin:" + os.environ["PATH"]))
    if r.returncode:
        raise RuntimeError(f"{cmd}\n{r.stderr}")
    return r.stdout


def grf_p(z, resels):
    """3-D EC-density GRF peak p (as FSL ptoz -g)."""
    return float(min(1.0, resels * (4 * np.log(2)) ** 1.5 * (2 * np.pi) ** -2 * (z ** 2 - 1) * np.exp(-z ** 2 / 2)))


def reml_meta(y, v):
    def nll(t2):
        w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
        return 0.5 * (np.log(v + t2).sum() + np.log(w.sum()) + (w * (y - mu) ** 2).sum())
    t2 = optimize.minimize_scalar(nll, bounds=(0, max(1e-8, 10 * y.var())), method="bounded").x
    w = 1 / (v + t2); mu = (w * y).sum() / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(mu=float(mu), se=float(se), z=float(mu / se), p_one=float(stats.norm.sf(mu / se)), tau2=float(t2))


def flame_tests(analysis, subs, variant, label):
    od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
    con = CON[analysis]
    eff = [fl(s, variant, con, "effect") for s in subs]
    var = [fl(s, variant, con, "variance") for s in subs]
    dofs = []
    for s in subs:
        js = json.load(open(os.path.join(g.OUT, "firstlevel", s, f"{s}_firstlevel.json")))
        if con == "C2joint":
            dofs.append(js["joint"][variant]["dof"])
        else:   # per-run FFX: sum of the dofs of the runs that formed C2
            dofs.append(sum(r["variants"][variant]["dof_A"] for r in js["runs"]
                            if variant in r["variants"] and "C2" in r["variants"][variant]["contrasts"]))
    ref = nib.load(g.GROUP_MASK)
    stack = lambda files: nib.Nifti1Image(np.stack([nib.load(f).get_fdata(dtype=np.float32) for f in files], -1),
                                          ref.affine)
    stack(eff).to_filename(f"{od}/cope.nii.gz"); stack(var).to_filename(f"{od}/varcope.nii.gz")
    m = ref.get_fdata() > 0
    nib.Nifti1Image(np.stack([m * d for d in dofs], -1).astype(np.float32), ref.affine).to_filename(f"{od}/dof.nii.gz")
    n = len(subs)
    if COV == "protocol":
        pv = protocol_of(subs); pv = pv - pv.mean()
        open(f"{od}/design.mat", "w").write(f"/NumWaves 2\n/NumPoints {n}\n/Matrix\n" +
                                            "".join(f"1 {x:.6f}\n" for x in pv))
        open(f"{od}/design.con", "w").write("/NumWaves 2\n/NumContrasts 1\n/Matrix\n1 0\n")
    else:
        open(f"{od}/design.mat", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
        open(f"{od}/design.con", "w").write("/NumWaves 1\n/NumContrasts 1\n/Matrix\n1\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
       f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
       f"--cs={od}/design.grp --runmode=flame1")
    st = f"{od}/stats"
    sm = sh(f"smoothest -d {n - (2 if COV else 1)} -r {st}/res4d.nii.gz -m {g.GROUP_MASK}")
    smv = {l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()}
    z = nib.load(f"{st}/zstat1.nii.gz").get_fdata()
    revar = nib.load(f"{st}/mean_random_effects_var1.nii.gz").get_fdata()
    varimg = nib.load(f"{od}/varcope.nii.gz").get_fdata()
    copeimg = nib.load(f"{od}/cope.nii.gz").get_fdata()
    rows, wrows = [], []
    for roi in SETS["H3"]:
        rm = nib.load(ROI(roi)).get_fdata() > 0
        resels = rm.sum() / smv["RESELS"]
        zr = np.where(rm, z, -np.inf); ijk = np.unravel_index(np.argmax(zr), z.shape)
        zmax = float(z[ijk]); p = grf_p(zmax, resels) if zmax > 1 else 1.0
        zthr = float(sh(f"ptoz {0.05 / 3} -g {resels}").strip())
        vmean = varimg[rm].mean(axis=0); re = float(revar[rm].mean())
        w = 1 / (vmean + re); share = w / w.sum()
        y = copeimg[rm].mean(axis=0)
        meta = reml_meta(y, vmean)
        rows.append(dict(analysis=analysis, label=label, variant=variant, roi=roi, N=n,
                         peak_z=round(zmax, 3), peak_mni=[round(float(x), 1) for x in nib.affines.apply_affine(ref.affine, ijk)],
                         resels=round(float(resels), 2), p_fwe=p, p_fwe_bonf=min(1.0, 3 * p),
                         z_thresh_bonf=round(zthr, 3), n_vox_sig=int((z[rm] > zthr).sum()),
                         roi_mean=float(y.mean()), roi_dz_unweighted=float(y.mean() / y.std(ddof=1)),
                         reml_mu=meta["mu"], reml_z=meta["z"], reml_p_one=meta["p_one"], reml_tau2=meta["tau2"],
                         max_weight_share=float(share.max())))
        for s_, sh_ in zip(subs, share):
            wrows.append(dict(roi=roi, subject=s_, weight_share=float(sh_)))
    # whole brain, two tails, cluster FWE by GRF, Bonferroni over tails
    cl = []
    for sgn, zf in ((+1, f"{st}/zstat1.nii.gz"), (-1, f"{st}/zneg.nii.gz")):
        if sgn < 0:
            sh(f"fslmaths {st}/zstat1.nii.gz -mul -1 {zf}")
        txt = sh(f"cluster -i {zf} -t 3.09 -p 1 -d {smv['DLH']} --volume={int(smv['VOLUME'])} --mm")
        lines = [l.split("\t") for l in txt.strip().splitlines()]
        hdr = lines[0]
        for l in lines[1:]:
            d = dict(zip(hdr, l))
            cl.append(dict(sign="+" if sgn > 0 else "-", size_vox=int(d["Voxels"]), p_cluster=float(d["P"]),
                           p_fwe_two_sided=min(1.0, 2 * float(d["P"])), peak_z=float(d["MAX"]),
                           peak_mni=[float(d["MAX X (mm)"]), float(d["MAX Y (mm)"]), float(d["MAX Z (mm)"])]))
    json.dump(dict(smoothest=smv, dofs=dofs, subjects=subs), open(f"{od}/flame_meta.json", "w"), indent=1)
    return pd.DataFrame(rows), pd.DataFrame(cl), pd.DataFrame(wrows)


def ols_s7(subs, variant, label):
    """S7's second fit: the v0.5 OLS one-sample model, same ROI/permutation machinery as H2/H4."""
    od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
    imgs = [nib.load(fl(s, variant, "C2joint", "effect")) for s in subs]
    dm = pd.DataFrame({"intercept": np.ones(len(subs))}); rows = []
    for roi in SETS["H3"]:
        res = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="intercept",
                                       mask=ROI(roi), n_perm=N_PERM, two_sided_test=False,
                                       random_state=SEED, n_jobs=N_JOBS)
        lp = res["logp_max_t"].get_fdata()[nib.load(ROI(roi)).get_fdata() > 0]
        coord, tpk = peak_mni(res["t"], ROI(roi), +1)
        pmin = float(10 ** -np.nanmax(lp))
        rows.append(dict(analysis="S7ols", label=label, variant=variant, roi=roi, N=len(subs),
                         peak_t=round(tpk, 3), peak_mni=coord, p_fwe_min=pmin, p_fwe_bonf=min(1, 3 * pmin)))
    return pd.DataFrame(rows)


def main(analysis, variant="base", subset="all"):
    sets, maia, _ = subject_sets()
    subs = subjects_for(analysis, variant, subset, sets)
    label = f"{analysis}_var-{variant}_set-{subset}" + (f"_atlas-{ATLAS}" if ATLAS != "caren" else "") \
        + (f"_cov-{COV}" if COV else "")
    rd = os.path.join(g.REPO, "Analysis", "results"); os.makedirs(rd, exist_ok=True)
    print(f"{label}: N = {len(subs)}")
    if analysis in ("H3", "H3perrun"):
        roi, cl, w = flame_tests(analysis, subs, variant, label)
        w.to_csv(f"{rd}/{label}_weights.csv", index=False)
    elif analysis == "S7ols":
        roi, cl = ols_s7(subs, variant, label), pd.DataFrame()
    else:
        roi, cl = nilearn_tests(analysis, subs, variant, label, maia)
    roi.to_csv(f"{rd}/{label}_roi.csv", index=False)
    cl.to_csv(f"{rd}/{label}_wholebrain.csv", index=False)
    json.dump(dict(label=label, N=len(subs), subjects=subs), open(f"{rd}/{label}_subjects.json", "w"), indent=1)
    print(f"{label}: wrote Analysis/results/{label}_*.csv")


if __name__ == "__main__":
    main(*sys.argv[1:])
