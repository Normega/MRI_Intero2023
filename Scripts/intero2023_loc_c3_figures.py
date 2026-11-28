#!/usr/bin/env python
"""
intero2023_loc_c3_figures.py -- localizer stage C3 (EXPLORATORY) figures, from the stage-C3 tables only.

  python Scripts/intero2023_loc_c3_figures.py

Figures/localizer/loc_C3_t4_decoding.png   T4 index by outcome; decoding balanced accuracy vs chance
Figures/localizer/loc_C3_rsa_gradient.png  crossnobis RDMs; insula C1 profile along y and slopes
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
RD = os.path.join(g.REPO, "Analysis", "results")
FIG = os.path.join(g.REPO, "Figures", "localizer")
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
RED, MID = "#e34948", "#f0efec"
INK, INK2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})


def strip(ax, groups, colors, labels, ylabel, rng):
    for k, (vals, col) in enumerate(zip(groups, colors)):
        x = k + rng.uniform(-0.17, 0.17, len(vals))
        ax.scatter(x, vals, s=14, color=col, alpha=.55, edgecolor="none", zorder=2)
        m = np.mean(vals); se = np.std(vals, ddof=1) / np.sqrt(len(vals))
        ax.errorbar(k + .32, m, yerr=1.96 * se, fmt="o", ms=6, color=INK, mfc=col, mec=INK, lw=1.4, capsize=0, zorder=3)
    ax.axhline(0, color=INK2, lw=.8, ls=(0, (3, 3)), zorder=1)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel); ax.set_xlim(-.6, len(groups) - .3)
    ax.yaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)


def fig_t4_decoding():
    E = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_epochs.csv"))
    M = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_mixed.csv"))
    D = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_decoding.csv"))
    rng = np.random.default_rng(1)
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.2), gridspec_kw=dict(hspace=.42, wspace=.28))
    order = [("aware", BLUE, "Aware"), ("unaware", ORANGE, "Unaware"), ("nochange", AQUA, "No Change"), ("free", YELLOW, "Free")]
    for ax, col, title, ylab in ((axs[0, 0], "dot_own", "A  T4 index (dot product, own C1 map)", "C1 expression (PSC)"),
                                 (axs[0, 1], "corr_own", "B  T4 index (spatial correlation, own C1 map)", "r(epoch pattern, C1 map)")):
        pm = E.groupby(["subject", "outcome"])[col].mean().unstack()
        groups = [pm[o].dropna().values for o, _, _ in order]
        strip(ax, groups, [c for _, c, _ in order], [f"{l}\n(n {len(gv)})" for (_, _, l), gv in zip(order, groups)], ylab, rng)
        ax.set_title(title, loc="left")
        r = M[(M["index"] == col) & (M.model == "registered")].set_index("term")
        f = r.loc["F: outcome"]; au = r.loc["contrast: aware - unaware"]; pf = M[(M["index"] == col) & (M.model == "paced_vs_free") & (M.term == "epoch_typepaced")].iloc[0]
        txt = (f"Outcome F({f.num_df:.0f}, {f.df:.0f}) = {f.estimate:.2f}, p = {f.p:.2f}\n"
               f"aware − unaware {au.estimate:+.3f} [{au.ci_lo:+.3f}, {au.ci_hi:+.3f}]\n"
               f"paced − free {pf.estimate:+.3f} [{pf.ci_lo:+.3f}, {pf.ci_hi:+.3f}]")
        ax.text(.02, .98, txt, transform=ax.transAxes, va="top", fontsize=7.6, color=INK2,
                bbox=dict(boxstyle="round,pad=.3", fc=SURF, ec=GRID))
        lo, hi = min(v.min() for v in groups), max(v.max() for v in groups)
        ax.set_ylim(lo - .08 * (hi - lo), hi + .42 * (hi - lo))
    masks = ["whole_brain", "SAL", "SMN", "DMN", "insula"]; mlab = ["Whole brain", "SAL", "SMN", "DMN", "Insula"]
    ax = axs[1, 0]
    for j, (feat, colr, lab) in enumerate((("raw", BLUE, "Raw PSC maps"), ("pattern", ORANGE, "Pattern only (map mean and scale removed)"))):
        d = D[D.features == feat].set_index("mask").loc[masks]
        x = np.arange(len(masks)) + (j - .5) * .28
        ax.errorbar(x, d.bal_acc, yerr=[d.bal_acc - d.ci_lo, d.ci_hi - d.bal_acc], fmt="o", ms=7, color=colr, lw=1.6,
                    capsize=0, label=lab, zorder=3)
        ax.scatter(x, d.null_q95, marker="_", s=160, color=INK2, lw=1.6, zorder=2,
                   label="Permutation null, 95th percentile" if j == 0 else None)
    ax.axhline(.5, color=INK2, lw=.8, ls=(0, (3, 3)))
    ax.set_xticks(range(len(masks))); ax.set_xticklabels(mlab)
    ax.set_ylabel("Balanced accuracy"); ax.set_ylim(.3, 1.02)
    ax.set_title("C  Intero vs Extero, leave-one-participant-out", loc="left")
    ax.yaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
    ax = axs[1, 1]
    d = D[D.features == "raw"].set_index("mask").loc[masks]
    shares = [("Match", d.match_called_intero, BLUE), ("Extero", 1 - d.extero_called_extero, ORANGE),
              ("Self", d.self_called_intero, AQUA)]
    for j, (lab, v, colr) in enumerate(shares):
        x = np.arange(len(masks)) + (j - 1) * .24
        ax.scatter(x, v, s=46, color=colr, edgecolor=SURF, lw=1.5, zorder=3, label=f"Held-out {lab} maps")
    for k, mk in enumerate(masks):
        ax.text(k + .24, d.self_called_intero.iloc[k] - .045, "S", ha="center", va="top", fontsize=7, color=INK2)
        ax.text(k - .24, d.match_called_intero.iloc[k] + .04, "M", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.axhline(.5, color=INK2, lw=.8, ls=(0, (3, 3)))
    ax.set_xticks(range(len(masks))); ax.set_xticklabels(mlab); ax.set_ylim(0, 1.13); ax.set_yticks(np.arange(0, 1.01, .2))
    ax.set_ylabel("Share classified as Intero")
    ax.set_title("D  Intero-vs-Extero classifier applied to other conditions (raw)", loc="left")
    ax.yaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
    ax.legend(loc="upper center", fontsize=7.5, frameon=False, ncol=3, handletextpad=.2, columnspacing=1.0)
    axs[1, 0].legend(loc="lower left", fontsize=7.5, frameon=False)
    fig.text(.01, .005, "Exploratory (OSF p5ztc, stage C3). A–B: points = participant means; black = mean ± 95% CI; text = "
             "registered mixed model (index ~ confidence + outcome + (1|participant)). C: 95% participant-bootstrap CI; chance = 0.5. "
             "D: classifier trained on the other participants' Intero vs Extero maps.",
             fontsize=7, color=INK2)
    out = os.path.join(FIG, "loc_C3_t4_decoding.png")
    fig.savefig(out, dpi=170, bbox_inches="tight"); plt.close(fig); print(out)


def fig_rsa_gradient():
    D = pd.read_csv(os.path.join(g.OUT, "loc_c3", "rsa", "distances_long.csv"))
    meta = json.load(open(os.path.join(RD, "EXPLORATORY_LOC_C3_rsa_models.json")))
    G = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient.csv"))
    S = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient_subjects.csv"))
    P = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_gradient_profile.csv"))
    rois = [("caren_AInsL", "AInsL"), ("caren_dACC", "dACC"), ("net_SAL", "SAL"), ("net_SMN", "SMN"), ("net_VIS", "VIS")]
    lab = ["I", "E", "M", "S", "R"]
    fig = plt.figure(figsize=(12.5, 6.6))
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 1.15], hspace=.32, wspace=.38)
    mean = D.groupby(["roi", "pair"]).d.mean()
    for k, (roi, nm) in enumerate(rois):
        ax = fig.add_subplot(gs[0, k])
        R = np.zeros((5, 5))
        for pr in meta["pairs"]:
            a, b = pr.split("-"); i, j = lab.index(a), lab.index(b)
            R[i, j] = R[j, i] = 100 * mean.loc[(roi, pr)]
        vmax = np.abs(R).max()
        im = ax.imshow(R, cmap=DIV, norm=TwoSlopeNorm(0, -vmax, vmax))
        for i in range(5):
            for j in range(5):
                if i != j:
                    ax.text(j, i, f"{R[i, j]:.1f}", ha="center", va="center", fontsize=6.3,
                            color=INK if abs(R[i, j]) < .6 * vmax else SURF)
        ax.set_xticks(range(5)); ax.set_yticks(range(5)); ax.set_xticklabels(lab); ax.set_yticklabels(lab)
        ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
        ax.set_title(("A  " if k == 0 else "") + f"{nm}", loc="left")
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.03); cb.ax.tick_params(labelsize=6.5); cb.outline.set_visible(False)
    ax = fig.add_subplot(gs[1, :3])
    ax.fill_between(P.y_bin_center, P.C1_mean - 1.96 * P.C1_se, P.C1_mean + 1.96 * P.C1_se, color=BLUE, alpha=.18, lw=0)
    ax.plot(P.y_bin_center, P.C1_mean, color=BLUE, lw=2, marker="o", ms=4)
    ax.axhline(0, color=INK2, lw=.8, ls=(0, (3, 3)))
    ax.set_xlabel("MNI y (mm), posterior → anterior"); ax.set_ylabel("C1 Intero − Extero (PSC)")
    ax.set_title("B  C1 along the insula (CAREN AInsL ∪ AInsR ∪ PIns; 4 mm bins; mean ± 95% CI; N 48)", loc="left", fontsize=9)
    ax.yaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[1, 3:])
    specs = [("primary", "C1", BLUE), ("intero_vs_rest", "Intero vs Rest", ORANGE), ("extero_vs_rest", "Extero vs Rest", AQUA)]
    rng = np.random.default_rng(2)
    groups = [S[S.spec == s].slope_per10mm.values for s, _, _ in specs]
    strip(ax, groups, [c for _, _, c in specs], [l for _, l, _ in specs], "Slope on y (PSC per 10 mm)", rng)
    for k, (s, _, _) in enumerate(specs):
        r = G.set_index("spec").loc[s]
        ax.text(k, ax.get_ylim()[1] if False else np.max(groups[k]) + .02, f"dz {r.dz:+.2f}\np {r.p_signflip:.2g}",
                ha="center", va="bottom", fontsize=7.2, color=INK2)
    ax.margins(y=.18)
    ax.set_title("C  Slope along y per participant", loc="left", fontsize=9)
    ax.set_xlabel("positive = higher anteriorly", color=INK2, fontsize=8)
    fig.text(.01, .01, "Exploratory (OSF p5ztc, stage C3). A: group mean cross-run crossnobis distance × 100 (Ledoit–Wolf noise-normalised, "
             "per voxel), N 44; I Intero, E Extero, M Match, S Self, R Rest (zero pattern). Colour scale per panel.",
             fontsize=7, color=INK2)
    out = os.path.join(FIG, "loc_C3_rsa_gradient.png")
    fig.savefig(out, dpi=170, bbox_inches="tight"); plt.close(fig); print(out)


if __name__ == "__main__":
    fig_t4_decoding()
    fig_rsa_gradient()
