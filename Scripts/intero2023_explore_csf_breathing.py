"""EXPLORATORY E11 (deviations.md, 2026-11-28): are the fMRIPrep CSF / white-matter
regressors locked to paced vs free epochs in the breathing runs?

Nuisance columns only; no BOLD. For each registered run, each tissue signal (percent of
its run mean) is regressed on the run's Model A design with csf and white_matter
removed, and the H1 C1 contrast (mean of paced conditions present minus free) is taken.
Run: ~/miniforge3/envs/neuro/bin/python Scripts/intero2023_explore_csf_breathing.py
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g          # noqa: E402
import intero2023_twomodel_qc as q         # noqa: E402

TISSUE = ["csf", "white_matter"]


def main():
    man = g.run_list()
    rows = []
    for r in man.itertuples():
        ev, cf, _ = q.read_run(r.subject, r.bids_entity, int(r.n_vol))
        tis = cf[TISSUE].astype(float)
        dm = q.model_a(ev, cf.drop(columns=TISSUE), int(r.n_vol))
        X = dm.values.astype(float)
        c = g.spec_vector(dm, g.c1_spec(dm))
        for col in TISSUE:
            y = 100 * (tis[col].values / tis[col].values.mean() - 1)
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            rows.append(dict(subject=r.subject, run=r.bids_entity, tissue=col,
                             C1=float(c @ b),
                             paced=float(np.mean([b[list(dm.columns).index(k)]
                                                  for k in q.PACED_SIDE if k in dm.columns])),
                             free=float(b[list(dm.columns).index("free")])))
    d = pd.DataFrame(rows)
    sub = d.groupby(["tissue", "subject"])[["C1", "paced", "free"]].mean().reset_index()
    rng = np.random.default_rng(20261128)
    out = []
    for t, gg in sub.groupby("tissue"):
        for k in ["C1", "paced", "free"]:
            v = gg[k].values
            m = v.mean()
            flips = rng.choice([-1, 1], size=(10000, len(v)))
            null = (flips * v).mean(axis=1)
            p = (1 + (np.abs(null) >= abs(m)).sum()) / 10001
            out.append(dict(tissue=t, term=k, N=len(v), mean_psc=m,
                            dz=m / v.std(ddof=1), pct_neg=100 * (v < 0).mean(), p_signflip=p))
    res = pd.DataFrame(out)
    res.to_csv(os.path.join(g.REPO, "Analysis", "results",
                            "EXPLORATORY_E11_csf_paced_free.csv"), index=False)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
