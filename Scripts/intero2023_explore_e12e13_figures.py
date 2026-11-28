#!/usr/bin/env python
"""
intero2023_explore_e12e13_figures.py -- figures for EXPLORATORY E12 (adherence gate) and E13 (awareness decoder).

  python Scripts/intero2023_explore_e12e13_figures.py

Reads Analysis/results/EXPLORATORY_E12_epochs.csv, EXPLORATORY_E13_participants.csv, EXPLORATORY_E13_decoder.csv and
the E13 null distributions in $OUT/group/EXPLORATORY_E13/. Writes Figures/EXPLORATORY_E12_gate.png and
Figures/EXPLORATORY_E13_decoder.png.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
import intero2023_explore_adherence as e12
RD = os.path.join(g.REPO, "Analysis", "results")
FIG = os.path.join(g.REPO, "Figures")
OD = os.path.join(g.OUT, "group", "EXPLORATORY_E13")
BLUE, ORANGE, GREEN, GOLD = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, SURF, MUTED, GRID = "#0b0b0b", "#fcfcfb", "#5c5c5a", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})


def fig_e12():
    P = pd.read_csv(os.path.join(RD, "EXPLORATORY_E12_epochs.csv"))
    D = P[P.plv.notna() & P.mae.notna()].copy()
    pr = D.groupby(["subject", "run"])
    D["plv_c"] = D.plv - pr.plv.transform("mean"); D["mae_c"] = D.mae - pr.mae.transform("mean")
    D["plv_r"] = e12.resid_within(D, "plv"); D["mae_r"] = e12.resid_within(D, "mae")
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=False)
    for ax, (x, y, ttl) in zip(axs, (("plv_c", "mae_c", "Centred within participant × run"),
                                     ("plv_r", "mae_r", "Also net of change type and direction"))):
        r = np.corrcoef(D[x], D[y])[0, 1]
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.scatter(D[x], D[y], s=14, color=BLUE, alpha=0.55, edgecolor=SURF, lw=0.5)
        b = np.polyfit(D[x], D[y], 1); xx = np.linspace(D[x].min(), D[x].max(), 2)
        ax.plot(xx, np.polyval(b, xx), color=ORANGE, lw=2)
        ax.set_title(ttl, fontsize=9, loc="left")
        ax.text(0.97, 0.95, f"r = {r:.2f}\ngate: r ≤ −0.40 → {'pass' if r <= -0.40 else 'fail'}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, color=INK)
        ax.set_xlabel("Paced PLV (deviation)")
    axs[0].set_ylabel("Period error, MAE s (deviation)")
    fig.suptitle(f"E12 Step 0: epoch PLV vs period error ({len(D)} epochs, {D.subject.nunique()} participants)",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "EXPLORATORY_E12_gate.png"), dpi=200); plt.close(fig)


def fig_e13():
    R = pd.read_csv(os.path.join(RD, "EXPLORATORY_E13_decoder.csv")).set_index("role")
    S = pd.read_csv(os.path.join(RD, "EXPLORATORY_E13_participants.csv"))
    nulls = [("primary", "null_primary.npy", "Brain (primary)", BLUE),
             ("control", "null_confreg.npy", "Brain, salience + direction removed", GREEN),
             ("control (benchmark)", "null_physio.npy", "Breathing only (benchmark)", ORANGE)]
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(width_ratios=[1.3, 1]))
    ax = axs[0]
    for k, (role, f, lab, col) in enumerate(nulls):
        nl = np.load(os.path.join(OD, f)); obs = R.loc[role, "mean_within_auc"]; p = R.loc[role, "p_one_sided"]
        lo, hi = np.percentile(nl, [2.5, 97.5])
        ax.plot([lo, hi], [k, k], color=GRID, lw=6, solid_capstyle="round")
        ax.plot([R.loc[role, "ci95_lo"], R.loc[role, "ci95_hi"]], [k, k], color=col, lw=2)
        ax.scatter([obs], [k], s=50, color=col, edgecolor=SURF, lw=1.5, zorder=3)
        ax.text(0.315, k + 0.28, f"{lab}: AUC {obs:.3f}, p {p:.3f}, N {int(R.loc[role, 'n_participants_scored'])}",
                fontsize=8, color=INK)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ax.set_yticks([]); ax.set_ylim(-0.5, 2.8); ax.set_xlim(0.30, 0.74)
    ax.set_xlabel("Mean within-participant AUC")
    ax.set_title("Observed (dot, 95% bootstrap CI)\nvs permutation null (grey, 95% range)", fontsize=9, loc="left")
    ax = axs[1]
    v = S.primary_auc.dropna().sort_values().values
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.scatter(v, np.arange(len(v)), s=14, color=BLUE, edgecolor=SURF, lw=0.5)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ax.set_yticks([]); ax.set_xlabel("Within-participant AUC (primary)")
    ax.set_title(f"Per participant (N {len(v)}), sorted", fontsize=9, loc="left")
    fig.suptitle("E13: decoding aware vs unaware paced epochs (SAL ∪ SMN LSS betas, LOSO)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "EXPLORATORY_E13_decoder.png"), dpi=200); plt.close(fig)


if __name__ == "__main__":
    fig_e12(); fig_e13()
