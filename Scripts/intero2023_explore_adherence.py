#!/usr/bin/env python
"""
intero2023_explore_adherence.py -- EXPLORATORY E12 (deviations.md, "After the update", 2026-11-28):
epoch-by-epoch pacer adherence in the breathing runs, and its brain correlates.

  python Scripts/intero2023_explore_adherence.py gate             # Step 0: behaviour only (no BOLD, no brain value)
  python Scripts/intero2023_explore_adherence.py qc               # Step 1 design QC (no BOLD)
  python Scripts/intero2023_explore_adherence.py subject sub-057  # Step 1 joint two-run first level
  python Scripts/intero2023_explore_adherence.py group            # Step 1 group tests (H1 masks, whole brain)
  python Scripts/intero2023_explore_adherence.py ridge            # Step 2 LOSO ridge on the T4 LSS betas
  Rscript Scripts/intero2023_explore_adherence_mixed.R            # Step 3 epoch mixed model on the T4 index

Measures (per paced epoch of the 90 registered runs):
  plv   paced PLV: belt breath phase vs the epoch's actual pacer trace, over the epoch from +3 s
        (stage C1, intero2023_loc_behaviour.breathing_runs; cache stageC1_measures.pkl key BEP).
        Exists only in runs whose belt and timing pass the C1 rule (58 runs, 33 participants).
  mae   adherence_mae_s (Data/resp_epochs.csv): mean |measured breath period - instructed period| over
        the 3-4 breaths wholly inside the epoch; belt grade good / degraded only. Lower = closer.
Change type (E7 rule, logged): salience high = step, low = ramp, No Change epochs = none.
Direction: the instructed change (Faster / Slower; none for No Change).

Step 0 gate (fixed in the E12 row; the plan's rule, made stricter): in paced epochs with both measures,
  r_raw     = Pearson r of plv and mae after centring each within participant x run (repeated-measures r);
  r_partial = the same after residualising each measure, within participant, on run intercepts, a
              No Change indicator, step vs ramp and direction (the covariate span of the Step 1 model).
  PASS iff r_raw <= -.40 AND r_partial <= -.40 (higher PLV must go with smaller period error).
  If it fails, E12 Steps 1-3 are not run.
Reported with it (descriptive): variance shares (between participant, between run within participant,
within run), odd/even split-half reliability of participant means (Spearman-Brown), run-1 vs run-2
ICC(3,1) of participant means.
"""
import os, sys, json, pickle, warnings
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
TAG = "E12"
CACHE = os.path.join(g.ROOT, "derivatives", "localizer_qc_cache", "stageC1_measures.pkl")
BELT_OK = ("good", "degraded")
GATE_R = -0.40
SEED = 20261128


# --------------------------------------------------------------------------- epoch table
def epoch_table():
    """Every paced epoch of the 90 registered runs, with both adherence measures and the stimulus codes."""
    man = g.run_list()
    ev = pd.concat([pd.read_csv(f"{g.REPO}/{q.FL}/{s}_{e}_desc-model_events.tsv", sep="\t", na_values=["n/a"])
                    .assign(subject=s, bids_entity=e) for s, e in zip(man.subject, man.bids_entity)], ignore_index=True)
    P = ev[ev.condition.isin(q.PACED_MOD)].copy()
    P["epoch_index"] = P.epoch_index.astype(int)
    B = pickle.load(open(CACHE, "rb"))["BEP"]
    B = B[B.kind == "paced"][["subject", "bids_entity", "epoch_index", "plv", "window_s"]].copy()
    B["epoch_index"] = B.epoch_index.astype(int)
    R = pd.read_csv(os.path.join(g.REPO, "Data", "resp_epochs.csv"))
    R = R[R.trial_type == "paced_breathing"][["subject", "bids_entity", "epoch_index", "belt_grade", "adherence_mae_s",
                                              "actual_ratio", "instructed_ratio", "n_contained"]]
    P = P.merge(B, on=["subject", "bids_entity", "epoch_index"], how="left", validate="1:1") \
         .merge(R, on=["subject", "bids_entity", "epoch_index"], how="left", validate="1:1")
    P["mae"] = P.adherence_mae_s.where(P.belt_grade.isin(BELT_OK))
    assert P.plv.isna().sum() == len(P) - 406 and P.plv.notna().sum() == 406
    assert not (P.plv.notna() & ~P.belt_grade.isin(BELT_OK)).any()
    P["run"] = P.bids_entity.str[-1].astype(int)
    P["change_type"] = np.where(P.condition == "paced_nochange", "none",
                                P.salience.map({"high": "step", "low": "ramp"}))
    assert P.change_type.notna().all()
    P["direction"] = np.where(P.condition == "paced_nochange", "none", P.change_direction)
    assert P.direction.isin(["Faster", "Slower", "none"]).all()
    return P.sort_values(["subject", "run", "epoch_index"]).reset_index(drop=True)


def covariates(d):
    """The Step 1 covariate span (within participant): run intercepts, No Change, step vs ramp, direction."""
    X = [(d.run == r).astype(float).values for r in sorted(d.run.unique())]
    X.append((d.change_type == "none").astype(float).values)
    X.append(np.select([d.change_type == "step", d.change_type == "ramp"], [0.5, -0.5], 0.0))
    X.append(np.select([d.direction == "Faster", d.direction == "Slower"], [0.5, -0.5], 0.0))
    return np.column_stack(X)


def resid_within(d, col):
    out = pd.Series(np.nan, index=d.index)
    for s, u in d.groupby("subject"):
        X = covariates(u); y = u[col].values
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        out[u.index] = y - X @ b
    return out


def sb(r):
    return 2 * r / (1 + r)


def icc31(a, b):
    X = np.c_[a, b]; n, k = X.shape
    gm = X.mean(); rm = X.mean(1); cm = X.mean(0)
    msr = k * ((rm - gm) ** 2).sum() / (n - 1)
    mse = ((X - rm[:, None] - cm[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse)


def gate():
    P = epoch_table()
    P.drop(columns=["onset", "duration", "confidence", "arousal", "response", "adherence_mae_s"], errors="ignore") \
     .round(6).to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_epochs.csv"), index=False)
    rows = []
    D = P[P.plv.notna() & P.mae.notna()].copy()
    pr = D.groupby(["subject", "run"])
    xc = D.plv - pr.plv.transform("mean"); yc = D.mae - pr.mae.transform("mean")
    n, k = len(D), pr.ngroups
    r = float(np.corrcoef(xc, yc)[0, 1]); df = n - k - 1
    p = float(2 * stats.t.sf(abs(r) * np.sqrt(df / (1 - r ** 2)), df))
    rows.append(dict(measure="gate: r(plv, mae), centred within participant x run", n_epochs=n,
                     n_participants=D.subject.nunique(), n_runs=k, r=r, df=df, p=p, criterion=f"<= {GATE_R}",
                     passes=r <= GATE_R))
    rp = resid_within(D, "plv"); rm_ = resid_within(D, "mae")
    npar = sum(covariates(u).shape[1] for _, u in D.groupby("subject"))
    r2 = float(np.corrcoef(rp, rm_)[0, 1]); df2 = n - npar - 1
    p2 = float(2 * stats.t.sf(abs(r2) * np.sqrt(df2 / (1 - r2 ** 2)), df2))
    rows.append(dict(measure="gate: r(plv, mae), residualised within participant on run, No Change, step-ramp, direction",
                     n_epochs=n, n_participants=D.subject.nunique(), n_runs=k, r=r2, df=df2, p=p2,
                     criterion=f"<= {GATE_R}", passes=r2 <= GATE_R))
    rho = stats.spearmanr(xc, yc).statistic
    rows.append(dict(measure="descriptive: Spearman of the participant x run centred values", n_epochs=n, r=float(rho)))
    xs = D.plv - D.groupby("subject").plv.transform("mean"); ys = D.mae - D.groupby("subject").mae.transform("mean")
    rows.append(dict(measure="descriptive: r(plv, mae), centred within participant only", n_epochs=n,
                     r=float(np.corrcoef(xs, ys)[0, 1])))
    mm = D.groupby("subject")[["plv", "mae"]].mean()
    rows.append(dict(measure="descriptive: r(plv, mae) of participant means (between participants)",
                     n_participants=len(mm), r=float(np.corrcoef(mm.plv, mm.mae)[0, 1])))
    # variance shares and reliability, each measure on all its epochs
    for col in ("plv", "mae"):
        A = P[P[col].notna()].copy()
        tot = ((A[col] - A[col].mean()) ** 2).sum()
        wp = ((A[col] - A.groupby("subject")[col].transform("mean")) ** 2).sum()
        wpr = ((A[col] - A.groupby(["subject", "run"])[col].transform("mean")) ** 2).sum()
        res = resid_within(A, col)
        rows.append(dict(measure=f"{col}: variance shares", n_epochs=len(A), n_participants=A.subject.nunique(),
                         share_between_participant=float(1 - wp / tot), share_between_run_within_participant=float((wp - wpr) / tot),
                         share_within_run=float(wpr / tot), share_within_after_covariates=float((res ** 2).sum() / tot),
                         mean=float(A[col].mean()), sd=float(A[col].std()), median=float(A[col].median()),
                         sd_within_run=float(np.sqrt(wpr / (len(A) - A.groupby(['subject', 'run']).ngroups)))))
        A["k"] = A.groupby("subject").cumcount()
        h = A.assign(half=A.k % 2).pivot_table(index="subject", columns="half", values=col, aggfunc="mean").dropna()
        rr = float(np.corrcoef(h[0], h[1])[0, 1])
        rows.append(dict(measure=f"{col}: odd/even split-half of participant means", n_participants=len(h), r=rr,
                         spearman_brown=sb(rr)))
        w = A.pivot_table(index="subject", columns="run", values=col, aggfunc="mean").dropna()
        rows.append(dict(measure=f"{col}: run 1 vs run 2 participant means, ICC(3,1)", n_participants=len(w),
                         icc31=float(icc31(w[1].values, w[2].values)), r=float(np.corrcoef(w[1], w[2])[0, 1])))
        # stimulus-driven part (mechanical): mean by change type x direction
        for (ct, di), u in A.groupby(["change_type", "direction"]):
            rows.append(dict(measure=f"{col}: mean by change type x direction", change_type=ct, direction=di,
                             n_epochs=len(u), mean=float(u[col].mean()), sd=float(u[col].std())))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_gate.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(out.round(3).to_string())
    print("\nGATE", "PASSES" if bool(out.passes.iloc[0]) and bool(out.passes.iloc[1]) else "FAILS")


if __name__ == "__main__":
    a = sys.argv[1]
    {"gate": gate}[a]()
