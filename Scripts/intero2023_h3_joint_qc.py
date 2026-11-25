#!/usr/bin/env python
"""
intero2023_h3_joint_qc.py -- design QC for the joint two-run H3 model (prereg v0.7).

H3 (aware > unaware) is fitted in ONE first-level model per participant that
spans both gate-passing runs:

  * task regressors (all Model A conditions + temporal derivatives) are SHARED
    across runs -- one aware column, one unaware column, ... stacked over runs;
  * nuisance regressors (base set: motion 6+6, FD, csf, white_matter, cosines,
    spikes, non-steady-state) and the intercept are RUN-SPECIFIC (block-diagonal).

Why: a run with aware trials but no unaware trial cannot yield aware - unaware
on its own, so the per-run fixed-effects form (v0.6) discards those aware
trials. In the joint model the shared free / rating / subrating regressors link
the runs, so every trial enters. Assumption: each condition's response is the
same in both runs (behaviourally supported; each run is scaled to percent signal
change before concatenation).

This script builds the joint design and the per-run designs from the event and
confound sheets only (no BOLD), and reports, per H3 participant, the variance
multiplier c'(X'X)^-1 c of aware - unaware under (a) the joint model and (b) the
per-run fixed-effects combination. `build_joint_design()` is the reference
construction for the analysis.

Output: Analysis/h3_joint_design_qc.csv
Usage:  python Scripts/intero2023_h3_joint_qc.py      (from the repo root)
Script version: 1.0 (2026-11-25)
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import intero2023_twomodel_qc as q      # read_run() drops physio; model_a() = Model A

warnings.simplefilter("ignore")
TASK_PREFIX = ("paced_", "free", "rating", "subrating")
C2 = {"paced_aware": 1.0, "paced_unaware": -1.0}


def cvar(dm, spec):
    """c'(X'X)^-1 c, column-standardised for numerical stability; nan if not estimable."""
    if not set(spec) <= set(dm.columns):
        return np.nan
    X = dm.values.astype(float)
    sd = X.std(axis=0); sd[sd < 1e-12] = 1.0
    Xs = X / sd
    c = np.array([spec.get(k, 0.0) for k in dm.columns]) / sd   # beta_s = sd * beta
    u, sv, vt = np.linalg.svd(Xs, full_matrices=False)
    rank = int((sv > max(Xs.shape) * np.finfo(float).eps * sv.max()).sum())
    if rank < Xs.shape[1]:
        Vr = vt[:rank]
        if np.linalg.norm(c - Vr.T @ (Vr @ c)) > 1e-8 * max(1.0, np.linalg.norm(c)):
            return np.nan
    return float(c @ np.linalg.pinv(Xs.T @ Xs, rcond=1e-12) @ c)


def build_joint_design(run_designs):
    """Stack per-run Model A designs: shared task columns, run-specific nuisance + intercept."""
    task = sorted(set().union(*[[c for c in dm.columns if c.startswith(TASK_PREFIX)]
                                for dm in run_designs]))
    parts = []
    for k, dm in enumerate(run_designs):
        nuis = dm.drop(columns=[c for c in dm.columns if c in task]).add_prefix(f"run{k + 1}_")
        parts.append((dm.reindex(columns=task, fill_value=0.0).reset_index(drop=True),
                      nuis.reset_index(drop=True)))
    ncols = [c for _, n in parts for c in n.columns]
    return pd.concat([pd.concat([t, n.reindex(columns=ncols, fill_value=0.0)], axis=1)
                      for t, n in parts], ignore_index=True)


def main():
    man = pd.concat([pd.read_csv(f"{q.FL}/firstlevel_manifest.csv"),
                     pd.read_csv(f"{q.FL}/firstlevel_manifest_candidates.csv")])
    man = man[(man.motion_gate_fail.astype(str).str.upper() == "FALSE") & (man.subject != "sub-056")]
    man = man.drop_duplicates(["subject", "bids_entity"])
    sets = pd.read_csv("Data/h3_sets.csv")
    rows = []
    for sb in sets.subject[sets.in_H3]:
        dms = []
        for _, r in man[man.subject == sb].sort_values("bids_entity").iterrows():
            ev, cf, _ = q.read_run(sb, r.bids_entity, int(r.n_vol))
            dms.append(q.model_a(ev, cf, int(r.n_vol)))
        per_run = [cvar(dm, C2) for dm in dms]
        ok = [v for v in per_run if np.isfinite(v)]
        v_ffx = 1.0 / sum(1.0 / v for v in ok) if ok else np.nan
        J = build_joint_design(dms)
        s = sets[sets.subject == sb].iloc[0]
        rows.append(dict(subject=sb, n_runs=len(dms), n_aware=int(s.n_aware_pooled),
                         n_unaware=int(s.n_unaware_pooled), n_reg_joint=J.shape[1],
                         var_joint=cvar(J, C2), var_perrun_ffx=v_ffx))
    d = pd.DataFrame(rows)
    d["ratio_joint_over_perrun"] = d.var_joint / d.var_perrun_ffx
    d.to_csv("Analysis/h3_joint_design_qc.csv", index=False)
    assert d.var_joint.notna().all(), "joint H3 contrast not estimable for some participant"
    gain = (1 / d.var_joint).sum() / (1 / d.var_perrun_ffx).sum() - 1
    print(f"H3 participants: {len(d)}; joint contrast estimable in all")
    print(f"variance ratio joint / per-run: median {d.ratio_joint_over_perrun.median():.3f} "
          f"[{d.ratio_joint_over_perrun.min():.3f}, {d.ratio_joint_over_perrun.max():.3f}]")
    print(f"total precision gain (sum of 1/var): {100 * gain:.1f}%")
    print("wrote Analysis/h3_joint_design_qc.csv")


if __name__ == "__main__":
    main()
