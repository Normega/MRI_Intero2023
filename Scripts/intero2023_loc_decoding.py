#!/usr/bin/env python
"""
intero2023_loc_decoding.py -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §10.3):
leave-one-subject-out decoding of Intero vs Extero from per-run condition maps, with
cross-decoding to Match vs Extero. Everything this script writes is exploratory.

Choices: Analysis/localizer_stageC3_implementation_notes.md §3 (committed before any value).

  python Scripts/intero2023_loc_decoding.py [n_jobs]

Samples: per-run 4 mm `intero` and `extero` effect maps (each vs Rest) of the localizer
Intero set (48 participants, 92 runs); Match and Self maps of the same runs for cross-decoding.
Masks: whole brain (group mask), CAREN SAL / SMN / DMN, insula (AInsL u AInsR u PIns).
Features: raw PSC (primary) and pattern (each map centred and scaled over the mask voxels).
Linear SVM (SVC on a precomputed linear kernel, C = 1, untuned), LOSO, balanced accuracy,
participant-bootstrap CI, 5,000 within-run label-swap permutations (seed 20261129).
Table: Analysis/results/EXPLORATORY_LOC_C3_decoding.csv (+ _participants.csv, _null.npz outside the repo).
"""
import os, sys, warnings
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
OD = os.path.join(g.OUT, "loc_c3", "decoding")
N_PERM, SEED, N_BOOT, BOOT_SEED, C = 5000, 20261129, 10000, 20261127, 1.0
CONDS = ["intero", "extero", "match", "C6"]       # C6 = Self vs Rest


def masks():
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    R = lambda n: nib.load(os.path.join(g.REPO, "ROIs", f"{n}.nii")).get_fdata() > 0
    M = {"whole_brain": gm, "SAL": R("net_SAL"), "SMN": R("net_SMN"), "DMN": R("net_DMN"),
         "insula": R("caren_AInsL") | R("caren_AInsR") | R("caren_PIns")}
    return {k: v[gm] for k, v in M.items()}, gm


def load_maps():
    s, _, _ = L.participant_sets()
    _, gm = masks()
    rows, X = [], []
    for x in s["intero"]:
        for k in ("run1", "run2"):
            ps = {c: L.fl_name(x, k, 4, "base", c, "effect") for c in CONDS}
            if not all(os.path.exists(p) for p in ps.values()):
                continue
            for c, p in ps.items():
                X.append(nib.load(p).get_fdata(dtype=np.float32)[gm])
                rows.append(dict(subject=x, run=k, condition=c))
    return pd.DataFrame(rows), np.array(X, np.float32)


def kernel(X, kind):
    X = X.astype(np.float64)
    if kind == "pattern":
        X = X - X.mean(1, keepdims=True); X = X / X.std(1, keepdims=True)
    return X @ X.T


def loso(K, meta, y_train, test_sets):
    """y_train: label per row (+1 intero / -1 extero / 0 not used). test_sets: name -> (rows mask, true label).
    Returns per-subject counts {name: (n_correct, n)}, plus predictions of every test row."""
    from sklearn.svm import SVC
    subj = meta.subject.values
    pred = np.zeros(len(meta))
    for x in np.unique(subj):
        tr = (subj != x) & (y_train != 0)
        te = subj == x
        clf = SVC(kernel="precomputed", C=C).fit(K[np.ix_(tr, tr)], y_train[tr])
        pred[te] = clf.predict(K[np.ix_(te, tr)])
    return pred


def bal_acc(pred, truth, rows):
    out = []
    for lab in (+1, -1):
        sel = rows & (truth == lab)
        if sel.any():
            out.append((pred[sel] == lab).mean())
    return float(np.mean(out))


def swap_labels(meta, rng):
    """Swap the Intero/Extero labels of each run with probability .5."""
    flip = {k: rng.random() < .5 for k in meta.groupby(["subject", "run"]).groups}
    f = np.array([flip[(a, b)] for a, b in zip(meta.subject, meta.run)])
    return f


def one_analysis(args):
    mask_name, kind, K, meta, n_vox, seed = args
    c = meta.condition.values
    base = np.where(c == "intero", 1, np.where(c == "extero", -1, 0)).astype(int)
    rows_ie = (c == "intero") | (c == "extero")
    # cross-decoding truth: match -> +1 (Intero class), extero -> -1
    truth_x = np.where(c == "match", 1, np.where(c == "extero", -1, 0))
    rows_x = (c == "match") | (c == "extero")
    pred = loso(K, meta, base, None)
    obs = bal_acc(pred, base, rows_ie)
    obs_x = bal_acc(pred, truth_x, rows_x)
    match_as_intero = float((pred[c == "match"] == 1).mean())
    self_as_intero = float((pred[c == "C6"] == 1).mean())
    # participant-level
    part = []
    for x, d in meta.groupby("subject"):
        ix = d.index.values
        part.append(dict(subject=x, mask=mask_name, features=kind,
                         acc_ie=float((pred[ix][rows_ie[ix]] == base[ix][rows_ie[ix]]).mean()),
                         n_ie=int(rows_ie[ix].sum()),
                         acc_cross=float((pred[ix][rows_x[ix]] == truth_x[ix][rows_x[ix]]).mean()),
                         match_as_intero=float((pred[ix][c[ix] == "match"] == 1).mean())))
    P = pd.DataFrame(part)
    rng = np.random.default_rng(BOOT_SEED)
    bs = np.array([P.acc_ie.values[rng.integers(0, len(P), len(P))].mean() for _ in range(N_BOOT)])
    bx = np.array([P.acc_cross.values[rng.integers(0, len(P), len(P))].mean() for _ in range(N_BOOT)])
    # permutations
    rng = np.random.default_rng(SEED + seed)
    null, null_x = np.empty(N_PERM), np.empty(N_PERM)
    for i in range(N_PERM):
        f = swap_labels(meta, rng)
        yp = np.where(f & rows_ie, -base, base)
        pp = loso(K, meta, yp, None)
        null[i] = bal_acc(pp, yp, rows_ie)          # within-decoding: train and test labels permuted
        null_x[i] = bal_acc(pp, truth_x, rows_x)     # cross-decoding: training labels permuted only
    os.makedirs(OD, exist_ok=True)
    np.savez_compressed(os.path.join(OD, f"null_{mask_name}_{kind}.npz"), null=null, null_x=null_x)
    rec = dict(mask=mask_name, features=kind, n_participants=meta.subject.nunique(),
               n_runs=meta.groupby(["subject", "run"]).ngroups, n_voxels=n_vox,
               bal_acc=obs, ci_lo=float(np.quantile(bs, .025)), ci_hi=float(np.quantile(bs, .975)),
               null_mean=float(null.mean()), null_q95=float(np.quantile(null, .95)),
               p_perm=float((1 + (null >= obs - 1e-12).sum()) / (N_PERM + 1)),
               share_participants_acc_gt_half=float((P.acc_ie > .5).mean()),
               cross_bal_acc=obs_x, cross_ci_lo=float(np.quantile(bx, .025)), cross_ci_hi=float(np.quantile(bx, .975)),
               cross_null_mean=float(null_x.mean()), cross_null_q95=float(np.quantile(null_x, .95)),
               cross_p_perm=float((1 + (null_x >= obs_x - 1e-12).sum()) / (N_PERM + 1)),
               match_called_intero=match_as_intero, extero_called_extero=float((pred[c == "extero"] == -1).mean()),
               self_called_intero=self_as_intero)
    return rec, P


def main(n_jobs=10):
    meta, X = load_maps()
    M, _ = masks()
    jobs = []
    for m in ("whole_brain", "SAL", "SMN", "DMN", "insula"):
        for k in ("raw", "pattern"):
            jobs.append((m, k, kernel(X[:, M[m]], k), meta, int(M[m].sum()), len(jobs)))
    del X
    with ProcessPoolExecutor(n_jobs) as ex:
        res = list(ex.map(one_analysis, jobs))
    T = pd.DataFrame([r for r, _ in res])
    T["p_bonf10"] = np.minimum(T.p_perm * len(T), 1)
    T["cross_p_bonf10"] = np.minimum(T.cross_p_perm * len(T), 1)
    T.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_decoding.csv"), index=False)
    pd.concat([p for _, p in res]).round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_decoding_participants.csv"), index=False)
    print(T[["mask", "features", "bal_acc", "ci_lo", "ci_hi", "p_perm", "cross_bal_acc", "cross_p_perm",
             "match_called_intero", "self_called_intero"]].round(3).to_string())


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
