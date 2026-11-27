"""Descriptive figures for the localizer conditions (OSF p5ztc, stage A group maps).

Figures (written to Figures/localizer/):
  loc_conditions_vs_rest.png  each task condition vs Rest (FLAME z, cluster FWE)
  loc_contrasts.png           pairwise condition contrasts + C3 conjunction
  loc_roi_conditions.png      ROI-mean percent signal change per condition vs Rest

Thresholding matches the registered whole-brain tables: cluster-forming z > 3.09,
GRF cluster p < .05 two-sided (p < .025 per tail; FSL `cluster` with the map's
smoothest DLH/volume). The C3 conjunction (min z of C1 and C2) has no GRF
cluster p, so it is shown at |z| > 3.09 with an extent of 20 voxels, descriptive.
These figures are descriptive and change no registered result.

Run: ~/miniforge3/envs/neuro/bin/python Scripts/intero2023_loc_figures.py
"""
import json
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import nibabel as nib
import numpy as np
from nilearn import datasets, image, plotting, surface
from nilearn.maskers import NiftiMasker

REPO = "/media/normega/FB_Cluster_Data/intero2024/MRI_Intero2023"
GLM = "/media/normega/FB_Cluster_Data/intero2024/derivatives/intero2023_glm"
GRP = os.path.join(GLM, "group")
FL = os.path.join(GLM, "firstlevel")
OUT = os.path.join(REPO, "Figures", "localizer")
os.makedirs(OUT, exist_ok=True)

INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
# Diverging blue <-> red with a neutral gray midpoint (reference palette).
DIV = LinearSegmentedColormap.from_list(
    "bgr", ["#0d366b", "#2a78d6", "#86b6ef", "#f0efec", "#f3a58a", "#e34948", "#8a1f1e"])
# Categorical slots 1-4 in fixed order.
COND_COL = {"Intero": "#2a78d6", "Extero": "#eb6834", "Match": "#1baf7a", "Self": "#eda100"}
VMAX = 8.0

FSAV = datasets.load_fsaverage("fsaverage5")
FSAV_DATA = datasets.load_fsaverage_data(mesh="fsaverage5", data_type="sulcal")


def cluster_masked(label, extent=None):
    """Return the FLAME z map masked to clusters surviving GRF p < .025 per tail."""
    d = os.path.join(GRP, label)
    meta = json.load(open(os.path.join(d, "loc_meta.json")))
    sm = meta["smoothest"]
    zimg = nib.load(os.path.join(d, "stats", "zstat1.nii.gz"))
    z = zimg.get_fdata()
    keep = np.zeros(z.shape, bool)
    n = {}
    with tempfile.TemporaryDirectory() as td:
        for sgn, tail in ((1, "pos"), (-1, "neg")):
            zf = os.path.join(td, f"z{tail}.nii.gz")
            nib.save(nib.Nifti1Image((sgn * z).astype(np.float32), zimg.affine), zf)
            idx = os.path.join(td, f"idx{tail}.nii.gz")
            subprocess.run(
                ["cluster", "-i", zf, "-t", "3.09", "-p", "0.025", "-d", str(sm["DLH"]),
                 f"--volume={int(sm['VOLUME'])}", f"--oindex={idx}", "--no_table"],
                check=True, capture_output=True)
            m = nib.load(idx).get_fdata() > 0
            keep |= m
            n[tail] = int(m.sum())
    out = np.where(keep, z, 0).astype(np.float32)
    return nib.Nifti1Image(out, zimg.affine), len(meta["subjects"]), n


def conjunction_masked():
    f = os.path.join(GRP, "LOC_C3conj_pos_minz.nii.gz")
    g = os.path.join(GRP, "LOC_C3conj_neg_minz.nii.gz")
    a, b = nib.load(f), nib.load(g)
    pos, neg = a.get_fdata(), b.get_fdata()
    # pos = min(z1, z2); neg = min(-z1, -z2) (see intero2023_loc_group.conjunction)
    z = np.where(pos > 3.09, pos, 0) - np.where(neg > 3.09, neg, 0)
    img = image.threshold_img(nib.Nifti1Image(z.astype(np.float32), a.affine),
                              threshold=3.09, cluster_threshold=20, two_sided=True,
                              copy_header=True)
    zz = img.get_fdata()
    return img, 48, {"pos": int((zz > 0).sum()), "neg": int((zz < 0).sum())}


def surf_row(fig, gs_row, img, title, ncol_offset=0):
    """Four surface views (L lat, L med, R med, R lat) + an axial strip."""
    tex = {h: surface.vol_to_surf(img, FSAV["pial"].parts[h], radius=3.0,
                                  interpolation="linear") for h in ("left", "right")}
    views = [("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]
    for j, (h, v) in enumerate(views):
        ax = fig.add_subplot(gs_row[j], projection="3d")
        plotting.plot_surf_stat_map(
            FSAV["inflated"].parts[h], tex[h], hemi=h, view=v,
            bg_map=FSAV_DATA.data.parts[h], bg_on_data=True,
            cmap=DIV, vmax=VMAX, threshold=3.09, colorbar=False, axes=ax, figure=fig)
        ax.set_facecolor(SURF)
    return tex


def map_figure(rows, fname, suptitle):
    nr = len(rows)
    fig = plt.figure(figsize=(16, 2.9 * nr + 0.9), facecolor=SURF)
    outer = fig.add_gridspec(nr, 2, width_ratios=[4.2, 3.2], hspace=0.35, wspace=0.02,
                             left=0.01, right=0.99, top=1 - 0.9 / (2.9 * nr + 0.9), bottom=0.05)
    for i, (title, img, n, counts) in enumerate(rows):
        g = outer[i, 0].subgridspec(1, 4, wspace=0.0)
        surf_row(fig, g, img, title)
        ax = fig.add_subplot(outer[i, 1])
        plotting.plot_stat_map(img, display_mode="z", cut_coords=[-10, 0, 10, 25, 40, 55],
                               cmap=DIV, vmax=VMAX, threshold=3.09, colorbar=(i == nr - 1),
                               annotate=False, black_bg=False, axes=ax, draw_cross=False)
        sub = (f"N = {n}   ·   red = first > second, blue = first < second"
               if "−" in title else f"N = {n}   ·   red = above Rest, blue = below Rest")
        fig.text(0.012, 1 - (i * 2.9 + 0.95) / (2.9 * nr + 0.9) + 0.005,
                 title, fontsize=13, fontweight="bold", color=INK, va="bottom")
        fig.text(0.012, 1 - (i * 2.9 + 0.95) / (2.9 * nr + 0.9) - 0.012,
                 sub, fontsize=9.5, color=INK2, va="top")
    fig.text(0.5, 0.012,
             "FLAME 1 z (6 mm). Cluster-forming z > 3.09, GRF cluster p < .05 two-sided. "
             "Axial slices z = −10, 0, 10, 25, 40, 55 mm (MNI). Colour scale ±8.",
             ha="center", fontsize=9, color=INK2)
    fig.suptitle(suptitle, fontsize=15, color=INK, x=0.012, ha="left", y=0.998)
    fig.savefig(os.path.join(OUT, fname), dpi=130, facecolor=SURF)
    plt.close(fig)


def roi_figure():
    """Per-participant mask means of each condition vs Rest (6 mm joint, Intero set)."""
    subs = json.load(open(os.path.join(GRP, "LOC_intero", "loc_meta.json")))["subjects"]
    rois = [("caren_AInsL", "Anterior insula L"), ("caren_AInsR", "Anterior insula R"),
            ("caren_dACC", "dACC"), ("caren_PIns", "Posterior insula"),
            ("ieat_ACCsphere", "IEAT ACC sphere"), ("caren_SM", "Sensorimotor")]
    cons = [("intero", "Intero"), ("extero", "Extero"), ("match", "Match"), ("C6", "Self")]
    vals = {}
    for r, _ in rois:
        mk = NiftiMasker(os.path.join(REPO, "ROIs", f"{r}.nii")).fit()
        for c, cn in cons:
            fs = [os.path.join(FL, s, f"{s}_task-localizer_loc-joint_sm-6_var-base_con-{c}_effect.nii.gz")
                  for s in subs]
            vals[(r, cn)] = mk.transform(fs).mean(axis=1)
    rows = []
    fig, axes = plt.subplots(1, len(rois), figsize=(16, 4.4), sharey=True, facecolor=SURF)
    rng = np.random.default_rng(1)
    for ax, (r, rn) in zip(axes, rois):
        ax.set_facecolor(SURF)
        ax.axhline(0, color="#b5b4ad", lw=1, zorder=0)
        for k, (_, cn) in enumerate(cons):
            v = vals[(r, cn)]
            m, se = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
            ax.scatter(k + rng.uniform(-0.18, 0.18, len(v)), v, s=9, color=COND_COL[cn],
                       alpha=0.35, lw=0, zorder=1)
            ax.errorbar(k, m, yerr=1.96 * se, fmt="o", ms=8, color=COND_COL[cn],
                        mec=SURF, mew=2, elinewidth=2, capsize=0, zorder=3)
            rows.append(dict(roi=r, condition=cn, N=len(v), mean_psc=m, ci95=1.96 * se,
                             dz=m / v.std(ddof=1)))
        ax.set_xticks(range(len(cons)), [c for _, c in cons], fontsize=10, color=INK)
        ax.set_title(rn, fontsize=11.5, color=INK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#b5b4ad")
        ax.tick_params(colors=INK2)
    axes[0].set_ylabel("% signal change vs Rest", fontsize=11, color=INK)
    allv = np.concatenate(list(vals.values()))
    pad = 0.05 * (allv.max() - allv.min())
    axes[0].set_ylim(allv.min() - pad, allv.max() + pad)
    fig.suptitle(f"Region means for each condition vs Rest (N = {len(subs)}, Intero set; "
                 "large dot = mean ± 95% CI, small dots = participants)",
                 fontsize=12.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUT, "loc_roi_conditions.png"), dpi=130, facecolor=SURF)
    plt.close(fig)
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(REPO, "Analysis", "results",
                                           "LOC_figure_roi_condition_means.csv"), index=False)


def main():
    log = {}
    cond = []
    for lab, title in (("LOC_intero", "Intero (track the breath) vs Rest"),
                       ("LOC_extero", "Extero (track the circle) vs Rest"),
                       ("LOC_match", "Match (breathe with the circle) vs Rest"),
                       ("LOC_C6", "Self (trait-word judgements) vs Rest")):
        img, n, c = cluster_masked(lab)
        cond.append((title, img, n, c)); log[lab] = c
    map_figure(cond, "loc_conditions_vs_rest.png", "Localizer conditions vs Rest")

    con = []
    for lab, title in (("LOC_L1", "Intero − Extero (C1; L1)"),
                       ("LOC_C2", "Intero − Self (C2)"),
                       ("LOC_L2", "Match − Intero (C5; L2)"),
                       ("LOC_C4", "Match − Extero (C4)")):
        img, n, c = cluster_masked(lab)
        con.append((title, img, n, c)); log[lab] = c
    img, n, c = conjunction_masked()
    con.append(("Conjunction: Intero − Extero AND Intero − Self (C3, min z, descriptive)", img, n, c))
    log["C3conj"] = c
    map_figure(con, "loc_contrasts.png", "Localizer condition contrasts")
    roi_figure()
    json.dump(log, open(os.path.join(OUT, "figure_voxel_counts.json"), "w"), indent=1)
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
