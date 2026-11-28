"""Descriptive response and adherence quality by localizer condition (localizer_deviations.md,
2026-11-28). Behaviour only; no inference.

Registered run set (102 runs, 54 participants), attempt of record. Per participant, then the
group median [IQR] and N. Extero / Match / Self keys come from the raw logs: each row stores
the LAST key in its window, so a stored button is a press made after that window's scanner
trigger (key 5), a selected subset.

Run: ~/miniforge3/envs/neuro/bin/python Scripts/intero2023_loc_response_summary.py
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L          # noqa: E402
import intero2023_localizer_qc as qc       # noqa: E402

RES = os.path.join(L.REPO, "Analysis", "results")
QCD = os.path.join(L.REPO, "Analysis", "localizer_qc")
BLK = {"ActiveExtero": "extero", "ActiveMatch": "match", "SelfRating": "self"}


def log_trials():
    u = L.run_list()
    rec = pd.read_csv(os.path.join(qc.CACHE, "events_qc_full.csv"))
    atts = qc.log_attempts()
    out = []
    for r in u.itertuples():
        sel = rec[(rec.subject == r.subject) & (rec.run == r.run)].iloc[0]
        A = [a for a in atts if a["pid"] == int(r.subject[4:]) and a["run"] == r.run
             and a["file"] == sel.selected_file and a["attempt"] == sel.selected_attempt]
        assert len(A) == 1
        A = A[0]["rows"]
        T = A[A.Block.isin(BLK) & (A.Description == "Trial")].copy()
        T = T[pd.to_numeric(T.TrialNum, errors="coerce") > 0]
        T["rt"] = pd.to_numeric(T.RT, errors="coerce")
        T["kind"] = np.where(T.Response == "5", "trigger",
                             np.where(T.Response.isin(["1", "2", "3", "4"]), "button", "none"))
        T["condition"] = T.Block.map(BLK)
        # circle tasks: phase onset, trigger grid, last trigger in the response window
        C = T.condition.isin(["extero", "match"])
        T.loc[C, "onset"] = (pd.to_datetime(T.loc[C, "StimOnset"])
                             - pd.to_datetime(T.loc[C, "RunStartTime"])).dt.total_seconds()
        T.loc[C, "win"] = T.loc[C, "t"] - T.loc[C, "onset"]
        f = T[C & (T.kind == "trigger")]
        ph = ((f.onset + f.rt) % 2.0).values
        phi = float(np.median(ph)) if len(ph) else np.nan
        grid_ok = len(ph) > 0 and float(np.subtract(*np.percentile(ph, [75, 25]))) <= 0.02
        c = (phi - T.onset) % 2.0
        T["c_last"] = np.where(c + 2 < T.win, c + 2, c)
        T["grid_ok"] = grid_ok
        T["correct"] = np.where((T.kind == "button") & C, pd.to_numeric(T.Accuracy, errors="coerce"), np.nan)
        T = T.assign(subject=r.subject, run=r.run, protocol_version=sel.protocol_version)
        out.append(T[["subject", "run", "protocol_version", "condition", "kind", "Response", "rt",
                      "c_last", "grid_ok", "correct"]])
    return pd.concat(out, ignore_index=True), u


def main():
    T, u = log_trials()
    subs = sorted(u.subject.unique())
    P = pd.DataFrame(index=pd.Index(subs, name="subject"))

    # --- Intero (existing tables)
    comp = pd.read_csv(os.path.join(QCD, "compliance.csv"))
    comp = comp.merge(u[["subject", "run"]], on=["subject", "run"])
    g = comp.groupby("subject")
    P["intero_presses_per_block"] = g.n_presses.median()
    P["intero_pct_compliant_blocks"] = 100 * g.compliant.mean()
    P["intero_ipi_cv"] = g.ipi_cv.median()
    P["intero_key_alternation"] = g.alternation.median()
    c1 = pd.read_csv(os.path.join(RES, "EXPLORATORY_LOC_C1_participants.csv")).set_index("subject")
    for a, b in [("trk_tracking_plv", "intero_tracking_plv"),
                 ("trk_tracking_plv_chance_mean", "intero_tracking_plv_chance"),
                 ("trk_pair_cv", "intero_pair_score_cv"), ("trk_pair_cv_chance", "intero_pair_score_chance"),
                 ("loc_extero_plv", "extero_breath_plv_to_circle"),
                 ("loc_match_plv", "match_breath_plv_to_circle"),
                 ("loc_null_plv", "virtual_pacer_null_plv"),
                 ("cd_paced_plv", "breathingrun_paced_plv"), ("cd_free_plv", "breathingrun_free_plv_virtual"),
                 ("adhere_r", "breathingrun_adhere_r")]:
        P[b] = c1[a]
    ta = pd.read_csv(os.path.join(QCD, "tracking_accuracy_runs.csv")).merge(u[["subject", "run"]], on=["subject", "run"])
    ta = ta[ta.usable_breath_timing.astype(bool)]
    P["intero_press_lead_s"] = ta.groupby("subject").press_lead_s.mean()

    # --- Extero / Match / Self (raw logs)
    for c in ["extero", "match", "self"]:
        g = T[T.condition == c].groupby("subject")
        for k in ["button", "trigger", "none"]:
            P[f"{c}_pct_stored_{k}"] = 100 * g.kind.apply(lambda s: (s == k).mean())
        P[f"{c}_n_trials"] = g.size()
        if c != "self":
            b = T[(T.condition == c) & (T.kind == "button")]
            P[f"{c}_stored_button_pct_correct"] = 100 * b.groupby("subject").correct.mean()
            bg = b[b.grid_ok]
            P[f"{c}_stored_button_rt_after_onset_s"] = bg.groupby("subject").rt.median()
            P[f"{c}_stored_button_rt_after_trigger_s"] = (bg.rt - bg.c_last).groupby(bg.subject).median()
    s = T[(T.condition == "self") & (T.kind == "button")]
    answers = s.Response.value_counts().to_dict()
    P.to_csv(os.path.join(RES, "EXPLORATORY_LOC_response_summary_participants.csv"))

    rows = []
    for col in P.columns:
        v = P[col].dropna()
        rows.append(dict(measure=col, N=len(v), median=v.median(), q25=v.quantile(.25),
                         q75=v.quantile(.75), min=v.min(), max=v.max()))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(RES, "EXPLORATORY_LOC_response_summary.csv"), index=False)
    pd.set_option("display.width", 200)
    print(S.round(3).to_string(index=False))
    print("Self stored answers:", answers, " trials pooled by condition x kind:")
    print(pd.crosstab(T.condition, T.kind))
    print("grid-ok runs:", T.groupby(["subject", "run"]).grid_ok.first().sum(), "of", T.groupby(["subject", "run"]).ngroups)
    v = T[(T.condition != "self") & (T.kind == "button")]
    print("pooled stored-button accuracy by condition:", v.groupby("condition").correct.mean().round(3).to_dict())


if __name__ == "__main__":
    main()
