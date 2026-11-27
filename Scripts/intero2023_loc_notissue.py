#!/usr/bin/env python
"""
intero2023_loc_notissue.py -- EXPLORATORY sensitivity (logged in Preregistration/localizer_deviations.md).
Localizer 6 mm joint base model refitted WITHOUT the fMRIPrep `csf` and `white_matter`
regressors (everything else identical), then FLAME 1 on the eight descriptive maps,
compared with the registered base maps.

  python Scripts/intero2023_loc_notissue.py first sub-060     # one participant
  python Scripts/intero2023_loc_notissue.py flames
  python Scripts/intero2023_loc_notissue.py summary
First levels: <sub>_task-localizer_loc-joint_sm-6_var-notissue_con-<c>_<kind>.nii.gz
Group:        $OUT/group/LOCnotissue_<con>/
"""
import os, sys, json, subprocess, tempfile, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
import intero2023_group as G
g = L.g
warnings.simplefilter("ignore")
GD = os.path.join(g.OUT, "group")
RD = os.path.join(g.REPO, "Analysis", "results")
VAR = "notissue"
CONS = {"intero": "intero", "extero": "all", "match": "all", "C6": "all",
        "C1": "intero", "C2": "intero", "C4": "all", "C5": "intero"}
BASE_DIR = {"C1": "LOC_L1", "C5": "LOC_L2"}
ROIS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere", "caren_SM"]
HO = os.path.join(os.environ.get("FSLDIR", "/home/normega/fsl"), "data", "atlases", "HarvardOxford",
                  "HarvardOxford-sub-maxprob-thr25-2mm.nii.gz")
HO_LAB = {"lat_ventricles": [3, 14], "cerebral_wm": [1, 12], "cerebral_cortex": [2, 13]}


def first(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    runs = L.run_list(); runs = runs[runs.subject == sub]
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    Y, dms = {}, []
    for r in runs.itertuples():
        y = masker.transform(g.bold_path(sub, L.entity(r.run)))
        Y[int(r.run)], _ = mean_scaling(y, axis=0)
        ev, _ = L.read_events(sub, r.run, int(r.n_vol))
        cf = L.read_confounds(sub, r.run, int(r.n_vol)).drop(columns=["csf", "white_matter"])
        dms.append((int(r.run), L.design(ev, cf, int(r.n_vol))))
    dm = L.build_joint_design([d for _, d in dms])
    Ym = np.vstack([Y[k] for k, _ in dms])
    specs = {c: s for c, s in L.CONTRASTS.items() if L.estimable(dm, s)}
    out, dof = fit(Ym, dm, specs)
    for c, (e, v) in out.items():
        assert (v > 0).all()
        g.save_float32(masker.inverse_transform(e), L.fl_name(sub, "joint", 6, VAR, c, "effect"))
        g.save_float32(masker.inverse_transform(v), L.fl_name(sub, "joint", 6, VAR, c, "variance"))
    json.dump(dict(subject=sub, dof=dof, n_reg=int(dm.shape[1]), contrasts=sorted(out)),
              open(os.path.join(L.fl_dir(sub), f"{sub}_task-localizer_notissue.json"), "w"), indent=1)
    print(sub, "done", flush=True)


def subjects(con):
    s, _, _ = L.participant_sets()
    return list(s[CONS[con]])


def run(con):
    subs = subjects(con); od = os.path.join(GD, f"LOCnotissue_{con}"); os.makedirs(od, exist_ok=True)
    base_meta = json.load(open(os.path.join(GD, BASE_DIR.get(con, f"LOC_{con}"), "loc_meta.json")))
    assert base_meta["subjects"] == subs, con
    ref = nib.load(g.GROUP_MASK); m = ref.get_fdata() > 0
    stack = lambda kind: nib.Nifti1Image(np.stack([nib.load(L.fl_name(x, "joint", 6, VAR, con, kind))
                                                   .get_fdata(dtype=np.float32) for x in subs], -1), ref.affine)
    stack("effect").to_filename(f"{od}/cope.nii.gz"); stack("variance").to_filename(f"{od}/varcope.nii.gz")
    dofs = [json.load(open(os.path.join(L.fl_dir(x), f"{x}_task-localizer_notissue.json")))["dof"] for x in subs]
    nib.Nifti1Image(np.stack([m * d for d in dofs], -1).astype(np.float32), ref.affine).to_filename(f"{od}/dof.nii.gz")
    n = len(subs)
    open(f"{od}/design.mat", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1.000000\n" * n)
    open(f"{od}/design.con", "w").write("/NumWaves 1\n/NumContrasts 1\n/Matrix\n1\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    sm = G.sh(f"smoothest -d {n - 1} -r {od}/stats/res4d.nii.gz -m {g.GROUP_MASK}")
    meta = dict(con=con, var=VAR, subjects=subs, dofs=dofs,
                smoothest={l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()})
    os.remove(f"{od}/stats/res4d.nii.gz")
    json.dump(meta, open(f"{od}/loc_meta.json", "w"), indent=1)
    return con


def flames():
    with ThreadPoolExecutor(8) as ex:
        for c in ex.map(run, CONS):
            print(c, "fitted", flush=True)


def cluster_masked(od):
    meta = json.load(open(f"{od}/loc_meta.json")); sm = meta["smoothest"]
    zimg = nib.load(f"{od}/stats/zstat1.nii.gz"); z = zimg.get_fdata()
    keep = np.zeros(z.shape); out = {}
    with tempfile.TemporaryDirectory() as td:
        for sgn, tail in ((1, "pos"), (-1, "neg")):
            zf, idx = f"{td}/z.nii.gz", f"{td}/i.nii.gz"
            nib.save(nib.Nifti1Image((sgn * z).astype(np.float32), zimg.affine), zf)
            subprocess.run(["cluster", "-i", zf, "-t", "3.09", "-p", "0.025", "-d", str(sm["DLH"]),
                            f"--volume={int(sm['VOLUME'])}", f"--oindex={idx}", "--no_table"],
                           check=True, capture_output=True)
            k = nib.load(idx).get_fdata() > 0
            out[tail] = int(k.sum()); keep += sgn * k
    return keep, out


def summary():
    from nilearn.maskers import NiftiMasker
    from nilearn import image
    ref = nib.load(g.GROUP_MASK); gm = ref.get_fdata() > 0
    ho = image.resample_to_img(HO, ref, interpolation="nearest", force_resample=True,
                               copy_header=True).get_fdata()
    masks = {"whole_brain": gm}
    masks.update({k: gm & np.isin(ho, v) for k, v in HO_LAB.items()})
    for r in ROIS:
        masks[r] = gm & (image.resample_to_img(os.path.join(g.REPO, "ROIs", f"{r}.nii"), ref,
                         interpolation="nearest", force_resample=True, copy_header=True).get_fdata() > 0)
    rows, maps = [], {}
    for con in CONS:
        bd = os.path.join(GD, BASE_DIR.get(con, f"LOC_{con}")); nd = os.path.join(GD, f"LOCnotissue_{con}")
        kb, cb = cluster_masked(bd); kn, cn = cluster_masked(nd)
        pos_base = kb > 0                          # base-model positive clusters
        b = nib.load(f"{bd}/cope.nii.gz").get_fdata(dtype=np.float32)
        n_ = nib.load(f"{nd}/cope.nii.gz").get_fdata(dtype=np.float32)
        N = b.shape[-1]
        for name, mk in list(masks.items()) + [("base_positive_clusters", pos_base)]:
            if mk.sum() == 0:
                continue
            bv, nv = b[mk].mean(0), n_[mk].mean(0)
            d = nv - bv
            rows.append(dict(con=con, N=N, region=name, n_vox=int(mk.sum()),
                             base_mean=bv.mean(), base_dz=bv.mean() / bv.std(ddof=1),
                             notissue_mean=nv.mean(), notissue_dz=nv.mean() / nv.std(ddof=1),
                             change_mean=d.mean(), change_dz=d.mean() / d.std(ddof=1),
                             base_vox_pos=cb["pos"], base_vox_neg=cb["neg"],
                             notissue_vox_pos=cn["pos"], notissue_vox_neg=cn["neg"],
                             base_pos_vox_still_pos=int(((kn > 0) & pos_base).sum())))
        print(con, N, "base", cb, "notissue", cn, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RD, "EXPLORATORY_LOC_notissue_comparison.csv"), index=False)


def figure():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from nilearn import plotting
    DIV = LinearSegmentedColormap.from_list(
        "bgr", ["#0d366b", "#2a78d6", "#86b6ef", "#f0efec", "#f3a58a", "#e34948", "#8a1f1e"])
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    names = {"intero": "Intero vs Rest", "extero": "Extero vs Rest", "match": "Match vs Rest",
             "C6": "Self vs Rest", "C1": "Intero − Extero", "C2": "Intero − Self",
             "C4": "Match − Extero", "C5": "Match − Intero"}
    fig = plt.figure(figsize=(14, 2.05 * len(CONS) + 1.0), facecolor=SURF)
    gs = fig.add_gridspec(len(CONS), 2, left=0.14, right=0.995, top=0.94, bottom=0.04, hspace=0.12, wspace=0.03)
    for i, con in enumerate(CONS):
        for j, (d, lab) in enumerate([(os.path.join(GD, BASE_DIR.get(con, f"LOC_{con}")), "Registered model (with CSF + WM regressors)"),
                                      (os.path.join(GD, f"LOCnotissue_{con}"), "Without CSF + WM regressors")]):
            keep, _ = cluster_masked(d)
            zimg = nib.load(f"{d}/stats/zstat1.nii.gz")
            img = nib.Nifti1Image(np.where(keep != 0, zimg.get_fdata(), 0).astype(np.float32), zimg.affine)
            ax = fig.add_subplot(gs[i, j])
            plotting.plot_stat_map(img, display_mode="z", cut_coords=[-5, 10, 22, 30, 50], cmap=DIV, vmax=8,
                                   threshold=3.09, colorbar=False, annotate=False, black_bg=False,
                                   axes=ax, draw_cross=False)
            if i == 0:
                ax.set_title(lab, fontsize=12, color=INK, pad=6)
        n = len(subjects(con))
        y = 0.94 - (i + 0.5) * (0.90 / len(CONS))
        fig.text(0.006, y + 0.008, names[con], fontsize=11.5, fontweight="bold", color=INK, va="center")
        fig.text(0.006, y - 0.014, f"N {n}", fontsize=9, color=INK2, va="center")
    fig.suptitle("Localizer maps with and without the CSF and white-matter nuisance regressors "
                 "(FLAME 1 z, cluster-forming z > 3.09, cluster p < .05 two-sided)",
                 fontsize=12.5, color=INK, x=0.006, ha="left", y=0.99)
    fig.text(0.5, 0.008, "Red = above Rest / first > second; blue = below. Axial z = −5, 10, 22, 30, 50 mm. "
             "Colour scale ±8. Exploratory sensitivity.", ha="center", fontsize=9, color=INK2)
    out = os.path.join(g.REPO, "Figures", "localizer", "loc_notissue_comparison.png")
    fig.savefig(out, dpi=120, facecolor=SURF); plt.close(fig); print(out)


if __name__ == "__main__":
    {"first": lambda: first(sys.argv[2]), "flames": flames, "summary": summary, "figure": figure}[sys.argv[1]]()
