#!/usr/bin/env python
"""
intero2023_loc_physio_relations.py -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §10.2):
physiology-informed between-participant relations that stage C1 did not already cover.

Choices: Analysis/localizer_stageC3_implementation_notes.md §2 (committed before any value).
Family (BH-FDR over 4), using stage C1's corr_row (Spearman, permutation p, bootstrap CI,
Pearson, JZS BF10):
  1 tracking PLV x C1 SAL-DMN background-FC difference (Intero - Extero; base, median, 4 s)
  2 tracking PLV x participant mean paced-epoch T4 index (dot, own template)
  3 tracking PLV x participant aware - unaware T4 index (>= 1 epoch of each)
  4 breathing-run adherence (adhere_r) x participant mean paced-epoch T4 index
  python Scripts/intero2023_loc_physio_relations.py
Table: Analysis/results/EXPLORATORY_LOC_C3_physio_relations.csv
"""
import os, sys, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
from intero2023_loc_behaviour import corr_row, bh
RD = os.path.join(g.REPO, "Analysis", "results")


def fc_c1():
    vals = {}
    for f in glob.glob(os.path.join(g.OUT, "loc_conn", "sub-*_locfc.csv")):
        d = pd.read_csv(f)
        d = d[(d.task == "background") & (d.physio == "base") & (d.summary == "median") & (d["shift"] == 4.0)
              & (d.pair == "SAL-DMN") & d.condition.isin(["intero", "extero"])]
        if set(d.condition) != {"intero", "extero"}:
            continue
        if ((d.n_epochs < 2) | (d.n_vol < 20)).any():          # the C2 volume rule
            continue
        z = d.set_index("condition").z
        vals[d.subject.iloc[0]] = z["intero"] - z["extero"]
    return pd.Series(vals)


def main():
    C1 = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C1_participants.csv")).set_index("subject")
    s, _, _ = __import__("intero2023_loc_common").participant_sets()
    E = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_epochs.csv"))
    P = E[E.epoch_type == "paced"]
    t4_mean = P.groupby("subject").dot_own.mean()
    au = P[P.outcome.isin(["aware", "unaware"])].pivot_table(index="subject", columns="outcome",
                                                            values="dot_own", aggfunc="mean").dropna()
    t4_au = au.aware - au.unaware
    fc = fc_c1()
    fc = fc[fc.index.isin(s["intero"])]
    trk = C1.trk_tracking_plv
    rows = [corr_row(trk, fc, x="tracking PLV", y="C1 SAL-DMN background FC (Intero - Extero)"),
            corr_row(trk, t4_mean, x="tracking PLV", y="T4 mean paced-epoch index (dot, own)"),
            corr_row(trk, t4_au, x="tracking PLV", y="T4 aware - unaware index (dot, own)"),
            corr_row(C1.adhere_r, t4_mean, x="breathing-run adherence (adhere_r)", y="T4 mean paced-epoch index (dot, own)")]
    T = pd.DataFrame(rows)
    T["q_fdr4"] = bh(T.p_perm.values)
    T.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_physio_relations.csv"), index=False)
    print(T.round(3).to_string())


if __name__ == "__main__":
    main()
