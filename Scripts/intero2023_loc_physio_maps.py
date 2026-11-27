#!/usr/bin/env python
"""
intero2023_loc_physio_maps.py -- EXPLORATORY (logged in Preregistration/localizer_deviations.md).
Physiology-corrected versions of the descriptive localizer group maps:
for each map, FLAME 1 on the base and on the +RVT (or +PetCO2) joint 6 mm first levels,
in the same participants. Then region means (base vs corrected) and a figure.

  python Scripts/intero2023_loc_physio_maps.py flames
  python Scripts/intero2023_loc_physio_maps.py summary
Maps: $OUT/group/LOCphys_<con>_<physio>_<var>/ (outside the repo).
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
CONS = {"intero": "intero", "extero": "all", "match": "all", "C6": "all",
        "C1": "intero", "C2": "intero", "C4": "all", "C5": "intero"}
PHYS = {"rvt": "rvt", "co2": "co2"}
ROIS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere", "caren_SM"]


def subjects(con, phys):
    s, _, _ = L.participant_sets()
    return [x for x in s[CONS[con]] if x in s[PHYS[phys]]
            and os.path.exists(L.fl_name(x, "joint", 6, phys, con, "effect"))]


def dofs(subs, var):
    return [json.load(open(os.path.join(L.fl_dir(x), f"{x}_task-localizer_firstlevel.json")))
            ["fits"][f"joint:sm6:{var}"]["dof"] for x in subs]


def odir(con, phys, var):
    return os.path.join(GD, f"LOCphys_{con}_{phys}_{var}")


def run(job):
    con, phys, var = job
    subs = subjects(con, phys); od = odir(con, phys, var); os.makedirs(od, exist_ok=True)
    meta = dict(con=con, physio=phys, var=var, subjects=subs)
    mf = os.path.join(od, "loc_meta.json")
    if os.path.exists(mf) and os.path.exists(f"{od}/stats/zstat1.nii.gz"):
        old = json.load(open(mf))
        if all(old.get(k) == meta[k] for k in meta):
            return job, "cached"
    ref = nib.load(g.GROUP_MASK); m = ref.get_fdata() > 0
    stack = lambda kind: nib.Nifti1Image(np.stack([nib.load(L.fl_name(x, "joint", 6, var, con, kind))
                                                   .get_fdata(dtype=np.float32) for x in subs], -1), ref.affine)
    stack("effect").to_filename(f"{od}/cope.nii.gz"); stack("variance").to_filename(f"{od}/varcope.nii.gz")
    nib.Nifti1Image(np.stack([m * d for d in dofs(subs, var)], -1).astype(np.float32),
                    ref.affine).to_filename(f"{od}/dof.nii.gz")
    n = len(subs)
    open(f"{od}/design.mat", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1.000000\n" * n)
    open(f"{od}/design.con", "w").write("/NumWaves 1\n/NumContrasts 1\n/Matrix\n1\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={g.GROUP_MASK} --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    sm = G.sh(f"smoothest -d {n - 1} -r {od}/stats/res4d.nii.gz -m {g.GROUP_MASK}")
    meta["smoothest"] = {l.split()[0]: float(l.split()[1]) for l in sm.strip().splitlines()}
    os.remove(f"{od}/stats/res4d.nii.gz")
    json.dump(meta, open(mf, "w"), indent=1)
    return job, "fitted"


def flames():
    jobs = [(c, p, v) for c in CONS for p in PHYS for v in ("base", p)]
    with ThreadPoolExecutor(int(os.environ.get("FLAME_JOBS", "10"))) as ex:
        for job, st in ex.map(run, jobs):
            print(job, st, flush=True)


def cluster_counts(od):
    meta = json.load(open(f"{od}/loc_meta.json")); sm = meta["smoothest"]
    zimg = nib.load(f"{od}/stats/zstat1.nii.gz"); z = zimg.get_fdata()
    keep = np.zeros(z.shape); out = {}
    with tempfile.TemporaryDirectory() as td:
        for sgn, tail in ((1, "pos"), (-1, "neg")):
            zf = f"{td}/z.nii.gz"; idx = f"{td}/i.nii.gz"
            nib.save(nib.Nifti1Image((sgn * z).astype(np.float32), zimg.affine), zf)
            subprocess.run(["cluster", "-i", zf, "-t", "3.09", "-p", "0.025", "-d", str(sm["DLH"]),
                            f"--volume={int(sm['VOLUME'])}", f"--oindex={idx}", "--no_table"],
                           check=True, capture_output=True)
            k = nib.load(idx).get_fdata() > 0
            out[tail] = int(k.sum()); keep += sgn * k
    nib.save(nib.Nifti1Image(np.where(keep != 0, z, 0).astype(np.float32), zimg.affine),
             f"{od}/zstat1_clusterFWE.nii.gz")
    return out


def summary():
    from nilearn.maskers import NiftiMasker
    gm = NiftiMasker(g.GROUP_MASK).fit()
    maskers = {r: NiftiMasker(os.path.join(g.REPO, "ROIs", f"{r}.nii")).fit() for r in ROIS}
    rows = []
    for con in CONS:
        for phys in PHYS:
            per = {}
            for var in ("base", phys):
                od = odir(con, phys, var)
                cc = cluster_counts(od)
                cope = nib.load(f"{od}/cope.nii.gz")
                reg = {"whole_brain": gm.transform(cope).mean(1)}
                reg.update({r: mk.transform(cope).mean(1) for r, mk in maskers.items()})
                per[var] = (cc, reg)
            n = len(json.load(open(f"{odir(con, phys, 'base')}/loc_meta.json"))["subjects"])
            for r in ["whole_brain"] + ROIS:
                b, c = per["base"][1][r], per[phys][1][r]
                d = c - b
                rows.append(dict(con=con, physio=phys, N=n, region=r,
                                 base_mean=b.mean(), base_dz=b.mean() / b.std(ddof=1),
                                 corr_mean=c.mean(), corr_dz=c.mean() / c.std(ddof=1),
                                 retained=c.mean() / b.mean() if b.mean() != 0 else np.nan,
                                 change_dz=d.mean() / d.std(ddof=1),
                                 base_vox_pos=per["base"][0]["pos"], base_vox_neg=per["base"][0]["neg"],
                                 corr_vox_pos=per[phys][0]["pos"], corr_vox_neg=per[phys][0]["neg"]))
            print(con, phys, n, per["base"][0], per[phys][0], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RD, "EXPLORATORY_LOC_physio_corrected_maps.csv"), index=False)


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
    cols = [("rvt", "base", "Base"), ("rvt", "rvt", "+RVT"), ("co2", "base", "Base"), ("co2", "co2", "+PetCO2")]
    fig = plt.figure(figsize=(18, 2.05 * len(CONS) + 1.0), facecolor=SURF)
    gs = fig.add_gridspec(len(CONS), 4, left=0.11, right=0.995, top=0.94, bottom=0.04,
                          hspace=0.12, wspace=0.03)
    for i, con in enumerate(CONS):
        for j, (phys, var, lab) in enumerate(cols):
            od = odir(con, phys, var)
            meta = json.load(open(f"{od}/loc_meta.json"))
            ax = fig.add_subplot(gs[i, j])
            plotting.plot_stat_map(f"{od}/zstat1_clusterFWE.nii.gz", display_mode="z",
                                   cut_coords=[-5, 10, 30, 50], cmap=DIV, vmax=8, threshold=3.09,
                                   colorbar=False, annotate=False, black_bg=False, axes=ax,
                                   draw_cross=False)
            if i == 0:
                ax.set_title(f"{lab}  ({'RVT' if phys == 'rvt' else 'PetCO2'} subset)",
                             fontsize=12, color=INK, pad=6)
        n_r = len(json.load(open(f"{odir(con, 'rvt', 'base')}/loc_meta.json"))["subjects"])
        n_c = len(json.load(open(f"{odir(con, 'co2', 'base')}/loc_meta.json"))["subjects"])
        y = 0.94 - (i + 0.5) * (0.90 / len(CONS))
        fig.text(0.006, y + 0.008, names[con], fontsize=11.5, fontweight="bold", color=INK, va="center")
        fig.text(0.006, y - 0.014, f"N {n_r} (RVT) / {n_c} (CO2)", fontsize=9, color=INK2, va="center")
    fig.suptitle("Localizer maps before and after breathing correction, same participants "
                 "(FLAME 1 z, cluster-forming z > 3.09, cluster p < .05 two-sided)",
                 fontsize=13, color=INK, x=0.006, ha="left", y=0.99)
    fig.text(0.5, 0.008, "Red = above Rest / first > second; blue = below Rest / first < second. "
             "Axial z = −5, 10, 30, 50 mm. Colour scale ±8. Exploratory.",
             ha="center", fontsize=9, color=INK2)
    out = os.path.join(g.REPO, "Figures", "localizer", "loc_physio_corrected.png")
    fig.savefig(out, dpi=120, facecolor=SURF); plt.close(fig)
    print(out)


if __name__ == "__main__":
    {"flames": flames, "summary": summary, "figure": figure}[sys.argv[1]]()
