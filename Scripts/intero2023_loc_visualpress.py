"""EXPLORATORY (localizer_deviations.md, 2026-11-28): circle-tracking presses recovered
from the censored Extero / Match logs, and their reliability. Behaviour only.

Each Extero / Match half-cycle row stores the FIRST key in a window opening at the phase
onset (StimOnset). The scanner trigger (key 5) arrives every 2 s on the same input, so a
real press is stored only if it beats the next trigger. The trigger grid is recovered per
run, which gives every trial a known censoring time; press latency is then estimated by
Kaplan-Meier (censoring set by the scanner, independent of behaviour).

Primary F1 = KM P(press within 1.0 s of phase onset); secondary RMST1.5 = restricted mean
time to press over 0-1.5 s. Gate: Extero F1 run-1 vs run-2 ICC(3,1) >= .40.

Run: ~/miniforge3/envs/neuro/bin/python Scripts/intero2023_loc_visualpress.py
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L          # noqa: E402
import intero2023_localizer_qc as qc       # noqa: E402  read_log / log_attempts / secs

RES = os.path.join(L.REPO, "Analysis", "results")
T_STAR, TAU = 1.0, 1.5
MIN_TRIALS = 10
BLOCKS = {"ActiveExtero": "extero", "ActiveMatch": "match"}
SEED = 20261128


# --------------------------------------------------------------------------- KM
def km(time, event):
    """Kaplan-Meier survivor as (unique event times, S just after each)."""
    time, event = np.asarray(time, float), np.asarray(event, bool)
    ts = np.unique(time[event])
    S, s = [], 1.0
    for t in ts:
        n = (time >= t).sum()
        d = ((time == t) & event).sum()
        s *= 1 - d / n
        S.append(s)
    return ts, np.array(S)


def surv_at(ts, S, t):
    k = np.searchsorted(ts, t, side="right")
    return 1.0 if k == 0 else S[k - 1]


def scores(df):
    if len(df) < MIN_TRIALS:
        return dict(n_trials=len(df), F1=np.nan, RMST15=np.nan)
    ts, S = km(df.time, df.event)
    grid = np.concatenate([[0], ts[ts < TAU], [TAU]])
    rmst = sum(surv_at(ts, S, a) * (b - a) for a, b in zip(grid[:-1], grid[1:]))
    return dict(n_trials=len(df), n_events=int(df.event.sum()),
                F1=1 - surv_at(ts, S, T_STAR), RMST15=rmst)


# --------------------------------------------------------------------------- trials
def run_trials(A, sub, run):
    T = A[A.Block.isin(BLOCKS) & (A.Description == "Trial")].copy()
    T = T[pd.to_numeric(T.TrialNum, errors="coerce") > 0]
    T["onset"] = (pd.to_datetime(T.StimOnset) - pd.to_datetime(T.RunStartTime)).dt.total_seconds()
    T["rt"] = pd.to_numeric(T.RT, errors="coerce")
    T["window"] = T.t - T.onset
    five = T[T.Response == "5"]
    ph = ((five.onset + five.rt) % 2.0).values
    phi = float(np.median(ph))
    iqr = float(np.subtract(*np.percentile(ph, [75, 25])))
    T["c_trig"] = (phi - T.onset) % 2.0
    dev = float(np.median(np.abs(five.rt - T.loc[five.index, "c_trig"])))
    T["press"] = T.Response.isin(["1", "2"])
    T["event"] = T.press
    T["time"] = np.where(T.press, T.rt,
                         np.where(T.Response == "5", T.rt, np.minimum(T.c_trig, T.window)))
    T["correct"] = np.where(T.press, pd.to_numeric(T.Accuracy, errors="coerce"), np.nan)
    T["condition"] = T.Block.map(BLOCKS)
    # Block number within the run (instruction rows delimit blocks)
    ins = A[(A.Description == "Instruction") & A.Block.isin(BLOCKS)]
    T["block_k"] = np.searchsorted(ins.index.values, T.index.values, side="right")
    qc_row = dict(subject=sub, run=run, n_trials=len(T), phi=phi, grid_iqr_s=iqr,
                  trig_latency_dev_s=dev,
                  n_press_after_trigger=int((T.press & (T.rt > T.c_trig + 0.01)).sum()),
                  n_none_with_trigger_in_window=int(((~T.press) & (T.Response != "5")
                                                     & (T.c_trig < T.window - 0.01)).sum()))
    qc_row["run_ok"] = (iqr <= 0.02) and (dev <= 0.01)
    T = T.assign(subject=sub, run=run)
    return T[["subject", "run", "condition", "block_k", "time", "event", "correct", "Response"]], qc_row


def icc(x, y):
    X = np.column_stack([x, y]); n, k = X.shape
    gm = X.mean()
    msr = k * ((X.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((X.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((X - X.mean(1, keepdims=True) - X.mean(0, keepdims=True) + gm) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse), (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


def main():
    u = L.run_list()
    rec = pd.read_csv(os.path.join(qc.CACHE, "events_qc_full.csv"))
    atts = qc.log_attempts()
    trials, qrows = [], []
    for r in u.itertuples():
        sel = rec[(rec.subject == r.subject) & (rec.run == r.run)].iloc[0]
        A = [a for a in atts if a["pid"] == int(r.subject[4:]) and a["run"] == r.run
             and a["file"] == sel.selected_file and a["attempt"] == sel.selected_attempt]
        assert len(A) == 1, (r.subject, r.run)
        t, q = run_trials(A[0]["rows"], r.subject, r.run)
        trials.append(t); qrows.append(q)
    T = pd.concat(trials, ignore_index=True)
    Q = pd.DataFrame(qrows)
    ok = Q[Q.run_ok][["subject", "run"]]
    print(f"runs {len(Q)}, grid OK {Q.run_ok.sum()}; grid IQR max {Q.grid_iqr_s.max():.4f} s; "
          f"latency dev max {Q.trig_latency_dev_s.max():.4f} s; presses after trigger "
          f"{Q.n_press_after_trigger.sum()}; no-key trials with trigger in window "
          f"{Q.n_none_with_trigger_in_window.sum()}")
    T = T.merge(ok, on=["subject", "run"])

    # per run x condition
    R = []
    for (s, rn, c), g in T.groupby(["subject", "run", "condition"]):
        R.append(dict(subject=s, run=rn, condition=c, **scores(g),
                      stored_press_rate=g.event.mean(), acc_stored=g.correct.mean()))
    R = pd.DataFrame(R).merge(Q, on=["subject", "run"], suffixes=("", "_run"))
    R.to_csv(os.path.join(RES, "EXPLORATORY_LOC_visualpress_runs.csv"), index=False)

    # per participant x condition (pooled runs), plus pooled Extero + Match
    P = []
    for (s, c), g in pd.concat([T, T.assign(condition="pooled")]).groupby(["subject", "condition"]):
        P.append(dict(subject=s, condition=c, n_runs=g.run.nunique(), **scores(g),
                      stored_press_rate=g.event.mean(), acc_stored=g.correct.mean()))
    P = pd.DataFrame(P)
    P.to_csv(os.path.join(RES, "EXPLORATORY_LOC_visualpress_participants.csv"), index=False)

    # reliability
    rng = np.random.default_rng(SEED)
    rel = []
    Rp = pd.concat([R, R.groupby(["subject", "run"]).apply(
        lambda g: pd.Series(scores(T[(T.subject == g.name[0]) & (T.run == g.name[1])])),
        include_groups=False).reset_index().assign(condition="pooled")], ignore_index=True)
    for c in ["extero", "match", "pooled"]:
        for m in ["F1", "RMST15"]:
            w = Rp[Rp.condition == c].pivot(index="subject", columns="run", values=m).dropna()
            x, y = w[1].values, w[2].values
            i3, i2 = icc(x, y)
            bs = []
            for _ in range(2000):
                k = rng.integers(0, len(x), len(x)); bs.append(icc(x[k], y[k])[0])
            # odd vs even block across runs
            g = T if c == "pooled" else T[T.condition == c]
            # block order within participant x condition (pooled: each condition split separately)
            g = g.copy()
            g["bkey"] = g.run * 100 + g.block_k
            g["ord"] = g.groupby(["subject", "condition"]).bkey.rank(method="dense").astype(int) - 1
            half = []
            for s, gg in g.groupby("subject"):
                a, b = scores(gg[gg.ord % 2 == 0])[m], scores(gg[gg.ord % 2 == 1])[m]
                half.append((a, b))
            h = np.array([v for v in half if np.all(np.isfinite(v))])
            rh = np.corrcoef(h[:, 0], h[:, 1])[0, 1]
            rel.append(dict(condition=c, measure=m, N_tworun=len(x), icc31=i3,
                            icc31_ci_lo=np.percentile(bs, 2.5), icc31_ci_hi=np.percentile(bs, 97.5),
                            icc21=i2, spearman_r1r2=pd.Series(x).corr(pd.Series(y), method="spearman"),
                            N_split=len(h), splithalf_r=rh, spearman_brown=2 * rh / (1 + rh),
                            mean_run1=x.mean(), mean_run2=y.mean()))
    rel = pd.DataFrame(rel)
    rel.to_csv(os.path.join(RES, "EXPLORATORY_LOC_visualpress_reliability.csv"), index=False)

    # descriptive relations with breath-tracking measures (behaviour only)
    c1 = pd.read_csv(os.path.join(RES, "EXPLORATORY_LOC_C1_participants.csv"))[
        ["subject", "trk_tracking_plv", "loc_match_plv"]]
    ta = pd.read_csv(os.path.join(L.REPO, "Analysis", "localizer_qc", "tracking_accuracy_runs.csv"))
    ta = ta.merge(u[["subject", "run"]], on=["subject", "run"])
    ta = ta[ta.usable_breath_timing.astype(bool)].groupby("subject").lag_median_s.mean().rename("trk_lag_median_s")
    W = P.pivot(index="subject", columns="condition", values=["F1", "RMST15"])
    W.columns = [f"{a}_{b}" for a, b in W.columns]
    W = W.join(c1.set_index("subject")).join(ta)
    cols = ["F1_extero", "F1_match", "F1_pooled", "RMST15_extero", "trk_tracking_plv",
            "trk_lag_median_s", "loc_match_plv"]
    cr = []
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            m = W[[a, b]].dropna()
            cr.append(dict(a=a, b=b, N=len(m), spearman=m[a].corr(m[b], method="spearman")))
    cr = pd.DataFrame(cr)
    for f in ["F1_extero", "F1_pooled"]:
        m = W[["trk_tracking_plv", f, "loc_match_plv"]].dropna()
        Z = (m - m.mean()) / m.std()
        Rm = np.corrcoef(Z.values.T)
        vif = np.diag(np.linalg.inv(Rm))
        for n_, v in zip(m.columns, vif):
            cr = pd.concat([cr, pd.DataFrame([dict(a=f"VIF[{f} model]", b=n_, N=len(m), spearman=np.nan, vif=v)])])
    cr.to_csv(os.path.join(RES, "EXPLORATORY_LOC_visualpress_correlations.csv"), index=False)

    # group KM curve (descriptive)
    for c in ["extero", "match"]:
        g = T[T.condition == c]; ts, S = km(g.time, g.event)
        print(c, "group F at 0.25/0.5/1.0/1.5 s:",
              [round(1 - surv_at(ts, S, t), 3) for t in (0.25, 0.5, 1.0, 1.5)],
              "acc stored", round(g.correct.mean(), 3))
    pd.set_option("display.width", 200)
    print(P.groupby("condition")[["n_trials", "F1", "RMST15", "stored_press_rate", "acc_stored"]].describe().T.round(3))
    print(rel.round(3).to_string(index=False))
    print(cr.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
