#!/usr/bin/env python
"""
intero2023_explore_awareness.py -- EXPLORATORY E3, E4, E5 (deviations.md, 2026-11-26).

  python Scripts/intero2023_explore_awareness.py e3subject sub-057   # E3 joint Model A, change - No Change
  python Scripts/intero2023_explore_awareness.py e3group             # E3 group map + awareness in its mask
  python Scripts/intero2023_explore_awareness.py e4                  # TFCE aware vs unaware
  python Scripts/intero2023_explore_awareness.py e5                  # Bayes factors in the H3 masks
  python Scripts/intero2023_explore_awareness.py e6                  # change - No Change in the H3 masks

E3a: joint two-run Model A (as H3: shared task columns, run-specific nuisance and intercept,
     PSC per run, AR(1)); contrast = unweighted mean of the change conditions present
     (paced_aware, paced_unaware, paced_wrongdir) minus paced_nochange.
E3b: mask = voxels of E3a clusters with p_FWE < .05 (either sign); the registered H3 FLAME z
     map tested in it, both tails, GRF voxel FWE (mask resels, H3 smoothness), Bonferroni x2.
E4:  TFCE, sign-flip, two-sided, 5,000 permutations, equal-weighted, on the H3 C2joint effects.
E5:  JZS Bayes factors (Rouder et al. 2009) on ROI-mean effects; r = 0.707 (0.5, 1.0 robustness).
E6:  E3a contrast in the three H3 masks: two-sided voxel max-t FWE (10,000 sign-flips), x3; mask-mean dz and BF10.
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import integrate, stats
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
from nilearn.glm.second_level import non_parametric_inference
import intero2023_glm_common as g
from intero2023_firstlevel import fit
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
CHANGE = ["paced_aware", "paced_unaware", "paced_wrongdir"]


def e3subject(sub):
    man = g.run_list(); runs = man[man.subject == sub]
    od = os.path.join(g.OUT, "firstlevel", sub)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    Ys, dms = [], []
    for r in runs.itertuples():
        Ys.append(mean_scaling(masker.transform(g.bold_path(sub, r.bids_entity)), axis=0)[0])
        ev, cf, _ = g.read_run_physio(sub, r.bids_entity, int(r.n_vol), [])
        dms.append(q.model_a(ev, cf, int(r.n_vol)))
    J = hj.build_joint_design(dms)
    present = [c for c in CHANGE if c in J.columns]
    log = dict(present=present)
    if "paced_nochange" in J.columns and present:
        spec = {c: 1.0 / len(present) for c in present}; spec["paced_nochange"] = -1.0
        out, dof = fit(np.vstack(Ys), J, {"E3": spec})
        e, v = out["E3"]
        g.save_float32(masker.inverse_transform(e), os.path.join(od, f"{sub}_var-base_con-E3changejoint_effect.nii.gz"))
        g.save_float32(masker.inverse_transform(v), os.path.join(od, f"{sub}_var-base_con-E3changejoint_variance.nii.gz"))
        log["dof"] = dof
    else:
        log["dof"] = None
    json.dump(log, open(os.path.join(od, f"{sub}_e3.json"), "w"))
    print(sub, log)


def e3group():
    import intero2023_group as G
    man = g.run_list(); subs = sorted(man.subject.unique())
    files = [G.fl(s, "base", "E3changejoint", "effect") for s in subs]
    subs = [s for s, f in zip(subs, files) if os.path.exists(f)]; files = [f for f in files if os.path.exists(f)]
    label = "EXPLORATORY_E3a_change_vs_nochange"
    od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
    dm = pd.DataFrame({"intercept": np.ones(len(files))})
    wb = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept",
                                  mask=g.GROUP_MASK, n_perm=G.N_PERM, two_sided_test=True,
                                  threshold=0.002, random_state=G.SEED, n_jobs=G.N_JOBS)
    for k in ("t", "logp_max_size"):
        g.save_float32(wb[k], os.path.join(od, f"wholebrain_{k}.nii.gz"))
    cl = G.cluster_table(wb["t"], wb["logp_max_size"], len(files) - 1)
    cl.to_csv(os.path.join(RD, f"{label}_wholebrain.csv"), index=False)
    json.dump(dict(N=len(files), subjects=subs), open(os.path.join(RD, f"{label}_subjects.json"), "w"))
    # E3b: awareness inside the change-responsive mask
    from scipy import stats as st
    t = wb["t"].get_fdata(); lp = wb["logp_max_size"].get_fdata()
    thr = st.t.isf(0.001, len(files) - 1)
    mask = (np.abs(t) > thr) & (lp > -np.log10(0.05))
    ref = nib.load(g.GROUP_MASK)
    mimg = nib.Nifti1Image(mask.astype(np.uint8), ref.affine); mimg.to_filename(os.path.join(od, "E3_change_mask.nii.gz"))
    h3 = os.path.join(g.OUT, "group", "H3_var-base_set-all")
    meta = json.load(open(os.path.join(h3, "flame_meta.json")))
    z = nib.load(os.path.join(h3, "stats", "zstat1.nii.gz")).get_fdata()
    rows = []
    if mask.sum():
        R = mask.sum() / meta["smoothest"]["RESELS"]
        for tail, sgn in (("aware>unaware", 1), ("aware<unaware", -1)):
            zz = np.where(mask, sgn * z, -np.inf); ijk = np.unravel_index(np.argmax(zz), z.shape)
            zmax = float(zz[ijk]); p = G.grf_p(zmax, R) if zmax > 1 else 1.0
            rows.append(dict(analysis="EXPLORATORY E3b", mask_vox=int(mask.sum()), resels=round(float(R), 2),
                             tail=tail, peak_z=round(zmax, 3),
                             peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(ref.affine, ijk)],
                             p_fwe_tail=p, p_fwe_bonf2=min(1.0, 2 * p),
                             n_vox_z_gt_thr=int((sgn * z[mask] > float(G.sh(f"ptoz {0.05 / 2} -g {R}").strip())).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(RD, "EXPLORATORY_E3b_awareness_in_change_mask.csv"), index=False)
    print(label, "N", len(files), "clusters p<.05:", int((cl.p_fwe < 0.05).sum()) if len(cl) else 0,
          "mask voxels", int(mask.sum()))


def e4():
    import intero2023_group as G
    h3 = json.load(open(os.path.join(RD, "H3_var-base_set-all_subjects.json")))["subjects"]
    files = [G.fl(s, "base", "C2joint", "effect") for s in h3]
    label = "EXPLORATORY_E4_tfce_aware_vs_unaware"
    od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
    dm = pd.DataFrame({"intercept": np.ones(len(files))})
    res = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept",
                                   mask=g.GROUP_MASK, n_perm=5000, two_sided_test=True, tfce=True,
                                   random_state=G.SEED, n_jobs=G.N_JOBS)
    for k in ("t", "tfce", "logp_max_tfce", "logp_max_t"):
        g.save_float32(res[k], os.path.join(od, f"{k}.nii.gz"))
    m = nib.load(g.GROUP_MASK).get_fdata() > 0
    lpt = res["logp_max_tfce"].get_fdata(); tt = res["t"].get_fdata(); ref = nib.load(g.GROUP_MASK)
    rows = []
    for sgn, tail in ((1, "aware>unaware"), (-1, "aware<unaware")):
        sel = m & (sgn * tt > 0)
        ijk = np.unravel_index(np.argmax(np.where(sel, lpt, -np.inf)), lpt.shape)
        rows.append(dict(analysis="EXPLORATORY E4", N=len(files), tail=tail,
                         min_p_fwe_tfce=float(10 ** -lpt[ijk]), t_at_min=round(float(tt[ijk]), 3),
                         mni=[round(float(c), 1) for c in nib.affines.apply_affine(ref.affine, ijk)],
                         n_vox_p_lt_05=int((sel & (lpt > -np.log10(0.05))).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(RD, f"{label}.csv"), index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def jzs_bf10(t, n, r):
    """One-sample JZS Bayes factor (Rouder et al., 2009, eq. 1), two-sided."""
    v = n - 1
    num = (1 + t ** 2 / v) ** (-(v + 1) / 2)
    f = lambda gg: (1 + n * gg * r ** 2) ** -0.5 * (1 + t ** 2 / ((1 + n * gg * r ** 2) * v)) ** (-(v + 1) / 2) \
        * (2 * np.pi) ** -0.5 * gg ** -1.5 * np.exp(-1 / (2 * gg))
    return integrate.quad(f, 0, np.inf)[0] / num


def bf_plus0(y, r):
    """Directional BF+0 = BF10 x 2 x posterior mass of delta > 0 (approximated by sampling the posterior)."""
    n = len(y); t = stats.ttest_1samp(y, 0).statistic
    bf10 = jzs_bf10(t, n, r)
    # posterior P(delta > 0): integrate the unnormalised posterior over delta with a Cauchy(0, r) prior
    d = np.linspace(-3, 3, 6001)
    like = stats.nct.pdf(t, n - 1, d * np.sqrt(n)); prior = stats.cauchy.pdf(d, 0, r)
    post = like * prior; post_pos = post[d > 0].sum() / post.sum()
    prior_pos = 0.5
    return bf10, bf10 * post_pos / prior_pos


def e5():
    import intero2023_group as G
    rows = []
    for target, con, sfile in (("H3", "C2joint", "H3_var-base_set-all_subjects.json"),
                               ("H2", "C3", "H2_var-base_set-all_subjects.json")):
        subs = json.load(open(os.path.join(RD, sfile)))["subjects"]
        for roi in G.SETS["H3"]:
            mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
            y = np.array([mk.transform(G.fl(s, "base", con, "effect")).mean() for s in subs])
            t = stats.ttest_1samp(y, 0).statistic
            for r in (0.5, 0.707, 1.0):
                b10, bp = bf_plus0(y, r)
                rows.append(dict(analysis="EXPLORATORY E5", target=target, roi=roi, N=len(y), dz=float(y.mean() / y.std(ddof=1)),
                                 t=float(t), prior_r=r, BF10=b10, BF01=1 / b10, BFplus0=bp, BF0plus=1 / bp))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(RD, "EXPLORATORY_E5_bayes_factors.csv"), index=False)
    print(d[d.prior_r == 0.707].round(3).to_string(index=False))

def e6():
    import intero2023_group as G
    subs = json.load(open(os.path.join(RD, "EXPLORATORY_E3a_change_vs_nochange_subjects.json")))["subjects"]
    files = [G.fl(s, "base", "E3changejoint", "effect") for s in subs]
    imgs = [nib.load(f) for f in files]
    dm = pd.DataFrame({"intercept": np.ones(len(files))})
    rows = []
    for roi in G.SETS["H3"]:
        rmask = nib.load(G.ROI(roi)).get_fdata() > 0
        lp = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="intercept",
                                      mask=G.ROI(roi), n_perm=G.N_PERM, two_sided_test=True,
                                      random_state=G.SEED, n_jobs=G.N_JOBS)
        tmap = G.ols_t(imgs, dm, "intercept", G.ROI(roi))
        lpv = lp.get_fdata()[rmask]; tv = tmap.get_fdata()[rmask]
        k = int(np.nanargmax(np.abs(tv))); ijk = np.argwhere(rmask)[k]
        pmin = float(10 ** -np.nanmax(lpv))
        mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
        y = np.array([mk.transform(f).mean() for f in files])
        row = dict(analysis="EXPLORATORY E6", roi=roi, N=len(files), peak_abs_t_signed=round(float(tv[k]), 3),
                   peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(tmap.affine, ijk)],
                   p_fwe=pmin, p_fwe_bonf3=min(1.0, 3 * pmin), n_vox_p_lt_05=int((lpv > -np.log10(0.05)).sum()),
                   mask_mean=float(y.mean()), mask_mean_dz=float(y.mean() / y.std(ddof=1)),
                   mask_mean_t=float(stats.ttest_1samp(y, 0).statistic))
        for r in (0.5, 0.707, 1.0):
            row[f"BF10_r{r}"] = jzs_bf10(row["mask_mean_t"], len(y), r)
        rows.append(row)
    d = pd.DataFrame(rows); d.to_csv(os.path.join(RD, "EXPLORATORY_E6_change_vs_nochange_h3masks.csv"), index=False)
    print(d.round(4).to_string(index=False))


if __name__ == "__main__":
    a = sys.argv[1]
    {"e3subject": lambda: e3subject(sys.argv[2]), "e3group": e3group, "e4": e4, "e5": e5, "e6": e6}[a]()
