#!/usr/bin/env python
"""
intero2023_explore_awareness_decoder.py -- EXPLORATORY E13 (deviations.md, "After the update", 2026-11-28):
one pre-specified multivariate test of aware vs unaware paced epochs in the breathing runs.

  python Scripts/intero2023_explore_awareness_decoder.py counts   # label-free counts and timing (no true label used)
  python Scripts/intero2023_explore_awareness_decoder.py run      # the logged tests

Data. The stage-C3 T4 LSS betas (Scripts/intero2023_loc_t4.py lss): one beta pattern per paced and free
epoch, swc7p Model A by-condition LSS, 4 mm, PSC, CAREN SAL u SMN (13,189 voxels), 44 participants
(localizer Intero set x breathing set). No new fitting; the mask is fixed a priori.

PRIMARY (one test).
  Features  each paced epoch's beta pattern minus the participant's mean pattern over ALL of their paced
            epochs (aware, unaware, wrong direction, No Change): label-free, removes the base-rate route.
  Samples   aware (paced_aware) vs unaware (paced_unaware) epochs; wrong direction and No Change left out.
  Model     linear SVM (sklearn SVC, precomputed linear kernel), C = 1, class_weight = 'balanced', no tuning,
            no feature scaling (as the stage-C3 decoder).
  CV        leave one participant out (every participant is a test fold once; participants with one
            class still train the other folds).
  Metric    ROC AUC of the decision values WITHIN each held-out participant with both classes, averaged
            over those participants (42). Pooled accuracy is not used.
  Null      5,000 permutations of the aware / unaware labels within participant x run, the whole LOSO
            re-run for each; one-sided p = (1 + #{null >= observed}) / 5,001. Seed 20261128.
  CI        participant bootstrap of the mean within-participant AUC (10,000 resamples, seed 20261129).
CONTROLS (pre-set, reported beside the primary, not corrected):
  C-conf    as the primary, with salience (step / ramp) and direction (Faster / Slower) regressed out of
            every voxel inside each training fold (OLS with intercept on the training epochs); the test
            participant's features are residualised with the training coefficients. Same null.
  C-phys    physiology-only benchmark: the same LOSO scheme and metric on four per-epoch breathing features:
            adherence_mae_s; the direction-aligned deviation of the measured from the instructed period
            ratio, s x (log actual_ratio - log instructed_ratio), s = +1 Slower, -1 Faster (positive = the
            breathing changed more than instructed); paced PLV (stage C1); epoch mean rvt_z from the
            model confound sheet (volumes whose mid-time falls in the epoch). Belt good / degraded only;
            complete cases. Each feature is centred on the participant's mean over their change epochs
            with complete features (label-free) and scaled by the training-fold SD. L2 logistic
            regression, C = 1, class_weight balanced. Same null and CI. The brain decoder's within-
            participant AUC on the same epochs and participants is given beside it, with a paired
            participant-bootstrap 95% CI of brain - physiology.
SECONDARY (Holm over the secondary family; m = 1 because the parcel decoder is not run):
  S-conf    within-participant z-scored confidence (1-7) of every paced epoch, regressed on the same
            centred features by linear kernel ridge, lambda fixed at the mean of the training kernel's
            diagonal (label-free, no tuning), no intercept (features and target are centred). LOSO.
            Metric: Spearman rho(predicted, observed) within each held-out participant, averaged.
            Null: 5,000 permutations of confidence within participant x run, LOSO re-run; one-sided.
            Participants with constant confidence are left out.
Not run: whole-brain parcel features (needs a whole-brain LSS refit and an atlas that is not on disk).
Outputs: Analysis/results/EXPLORATORY_E13_*.csv. Per-epoch decision values stay outside the repo
($OUT/group/EXPLORATORY_E13/).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd
from scipy import stats
from joblib import Parallel, delayed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
q = g.q
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
TAG = "E13"
T4D = os.path.join(g.OUT, "loc_c3", "t4")
OD = os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}")
N_PERM, SEED, BOOT_SEED, N_BOOT = 5000, 20261128, 20261129, 10000
N_JOBS = int(os.environ.get("E13_NJOBS", "8"))


# --------------------------------------------------------------------------- data
def load():
    """Centred paced-epoch patterns (all paced epochs) and the epoch table."""
    subs = sorted(f[:-8] for f in os.listdir(T4D) if f.endswith("_lss.npz"))
    assert len(subs) == 44, len(subs)
    X, rows = [], []
    for s in subs:
        z = np.load(os.path.join(T4D, f"{s}_lss.npz")); B = z["beta"].astype(np.float64)
        E = pd.read_csv(os.path.join(T4D, f"{s}_lss_epochs.csv"))
        assert len(E) == len(B)
        pk = E.condition.isin(q.PACED_MOD).values
        Bp = B[pk] - B[pk].mean(0, keepdims=True)
        X.append(Bp); rows.append(E[pk])
    E = pd.concat(rows, ignore_index=True)
    E["epoch_index"] = E.epoch_index.astype(int)
    X = np.vstack(X)
    assert len(E) == 595 and X.shape == (595, 13189)
    return X, E


def physio_features(E):
    """Per-epoch breathing features (see docstring); NaN where unavailable."""
    import intero2023_explore_adherence as e12
    A = e12.epoch_table()[["subject", "bids_entity", "epoch_index", "plv", "mae", "actual_ratio", "instructed_ratio",
                           "belt_grade"]].rename(columns={"bids_entity": "entity"})
    F = E[["subject", "entity", "epoch_index", "condition", "change_direction", "onset", "duration"]].merge(
        A, on=["subject", "entity", "epoch_index"], how="left", validate="1:1")
    ok = F.belt_grade.isin(e12.BELT_OK)
    s = F.change_direction.map({"Slower": 1.0, "Faster": -1.0})
    F["ratio_dev"] = (s * (np.log(F.actual_ratio) - np.log(F.instructed_ratio))).where(ok)
    man = g.run_list().set_index(["subject", "bids_entity"])
    rvt = []
    for r in F.itertuples():
        m = man.loc[(r.subject, r.entity)]
        if str(m.has_rvt).upper() != "TRUE":
            rvt.append(np.nan); continue
        cf = pd.read_csv(f"{g.REPO}/{q.FL}/{r.subject}_{r.entity}_desc-model_confounds.tsv", sep="\t")
        t = g.TR * (np.arange(len(cf)) + 0.5)
        w = (t >= r.onset) & (t < r.onset + r.duration)
        rvt.append(float(cf.rvt_z.values[w].mean()) if w.any() else np.nan)
    F["rvt"] = rvt
    return F[["plv", "mae", "ratio_dev", "rvt"]].reset_index(drop=True)


# --------------------------------------------------------------------------- metrics
def within_auc(y, d, groups):
    from sklearn.metrics import roc_auc_score
    out = {}
    for s in np.unique(groups):
        k = groups == s
        if 0 < y[k].sum() < k.sum():
            out[s] = roc_auc_score(y[k], d[k])
    return out


def within_spearman(y, p, groups, min_n=3):
    out = {}
    for s in np.unique(groups):
        k = groups == s
        if k.sum() >= min_n and np.std(y[k]) > 0 and np.std(p[k]) > 0:
            out[s] = stats.spearmanr(y[k], p[k]).statistic
    return out


def perm_labels(y, strata, rng):
    y2 = y.copy()
    for st in np.unique(strata):
        k = np.flatnonzero(strata == st)
        y2[k] = y[rng.permutation(k)]
    return y2


def boot_ci(v, seed=BOOT_SEED, n=N_BOOT):
    rng = np.random.default_rng(seed); v = np.asarray(v)
    m = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(n)])
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


# --------------------------------------------------------------------------- SVM LOSO on fold kernels
def fold_kernels(X, groups, conf=None):
    """Per held-out participant: (train idx, test idx, K_train, K_test). With conf, the training-fold
    OLS of every feature on [1, conf] is removed from train and test features."""
    folds = []
    for s in np.unique(groups):
        te = np.flatnonzero(groups == s); tr = np.flatnonzero(groups != s)
        Xtr, Xte = X[tr], X[te]
        if conf is not None:
            Ctr = np.c_[np.ones(len(tr)), conf[tr]]; Cte = np.c_[np.ones(len(te)), conf[te]]
            Bc = np.linalg.lstsq(Ctr, Xtr, rcond=None)[0]
            Xtr = Xtr - Ctr @ Bc; Xte = Xte - Cte @ Bc
        folds.append((tr, te, Xtr @ Xtr.T, Xte @ Xtr.T))
    return folds


def svm_decision(folds, y, n):
    from sklearn.svm import SVC
    d = np.full(n, np.nan)
    for tr, te, Ktr, Kte in folds:
        clf = SVC(kernel="precomputed", C=1.0, class_weight="balanced").fit(Ktr, y[tr])
        d[te] = clf.decision_function(Kte)
    return d


def logit_decision(F, y, groups):
    from sklearn.linear_model import LogisticRegression
    d = np.full(len(y), np.nan)
    for s in np.unique(groups):
        te = groups == s; tr = ~te
        sd = F[tr].std(0, ddof=1); sd[sd == 0] = 1.0
        clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000).fit(F[tr] / sd, y[tr])
        d[te] = clf.decision_function(F[te] / sd)
    return d


def run_perms(fn, y, strata, n_perm=N_PERM, seed=SEED):
    """fn(y) -> scalar; returns the null distribution (permutations within strata, drawn in order from
    one seeded generator, then evaluated in N_JOBS chunks)."""
    rng = np.random.default_rng(seed)
    perms = [perm_labels(y, strata, rng) for _ in range(n_perm)]
    chunks = np.array_split(np.arange(n_perm), N_JOBS)
    res = Parallel(n_jobs=N_JOBS)(delayed(lambda ix: [fn(perms[i]) for i in ix])(ix) for ix in chunks)
    return np.array([v for r in res for v in r])


# --------------------------------------------------------------------------- kernel ridge (secondary)
def kridge_all(X, Y, groups):
    """LOSO linear kernel ridge for many targets at once (columns of Y); lambda = mean diag(K_train)."""
    P = np.full(Y.shape, np.nan)
    for s in np.unique(groups):
        te = np.flatnonzero(groups == s); tr = np.flatnonzero(groups != s)
        K = X[tr] @ X[tr].T; lam = float(np.mean(np.diag(K)))
        w, V = np.linalg.eigh(K)
        A = V @ ((V.T @ Y[tr]) / (w + lam)[:, None])
        P[te] = (X[te] @ X[tr].T) @ A
    return P


# --------------------------------------------------------------------------- commands
def sets(E):
    A = E.condition.isin(["paced_aware", "paced_unaware"]).values
    return A


def counts():
    """Label-free: counts, feature availability, and timing of one LOSO pass with RANDOM labels."""
    import time
    X, E = load(); A = sets(E)
    Ea = E[A]
    c = Ea.groupby(["subject", "condition"]).size().unstack(fill_value=0)
    print(f"aware/unaware epochs {A.sum()} in {Ea.subject.nunique()} participants; "
          f"both classes: {((c > 0).all(1)).sum()}; aware {int(c.paced_aware.sum())}, unaware {int(c.paced_unaware.sum())}")
    F = physio_features(E)
    cc = F.notna().all(1).values & A
    c2 = E[cc].groupby(["subject", "condition"]).size().unstack(fill_value=0)
    print("feature availability on aware/unaware epochs:", F[A].notna().sum().to_dict())
    print(f"physiology complete cases: {cc.sum()} epochs, {E[cc].subject.nunique()} participants, both classes {((c2 > 0).all(1)).sum()}")
    cv = E.groupby("subject").confidence.std()
    print(f"confidence: paced epochs {len(E)}, participants {E.subject.nunique()}, constant: {sorted(cv[cv == 0].index)}")
    rng = np.random.default_rng(1)
    yr = rng.integers(0, 2, A.sum())
    t0 = time.time(); folds = fold_kernels(X[A], E.subject.values[A]); t1 = time.time()
    svm_decision(folds, yr, A.sum()); t2 = time.time()
    print(f"timing (random labels): kernels {t1 - t0:.1f} s, one LOSO pass {t2 - t1:.2f} s")


def run():
    X, E = load(); A = sets(E)
    Xa, Ea = X[A], E[A].reset_index(drop=True)
    y = (Ea.condition == "paced_aware").astype(int).values
    grp = Ea.subject.values
    strata = (Ea.subject + "_" + Ea.run.astype(str)).values
    n = len(y); rows = []; os.makedirs(OD, exist_ok=True)

    def summarise(name, role, aucs, null, extra=None):
        v = np.array(list(aucs.values())); obs = float(v.mean()); lo, hi = boot_ci(v)
        r = dict(analysis=f"EXPLORATORY {TAG}", test=name, role=role, n_participants_scored=len(v),
                 mean_within_auc=obs, ci95_lo=lo, ci95_hi=hi, median_within_auc=float(np.median(v)),
                 n_above_05=int((v > 0.5).sum()), null_mean=float(null.mean()), null_sd=float(null.std()),
                 null_95th=float(np.percentile(null, 95)), p_one_sided=float((1 + (null >= obs).sum()) / (1 + len(null))),
                 n_perm=len(null))
        if extra: r.update(extra)
        rows.append(r); print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in r.items()}, flush=True)
        return obs

    # PRIMARY
    folds = fold_kernels(Xa, grp)
    d = svm_decision(folds, y, n)
    auc = within_auc(y, d, grp)
    null = run_perms(lambda yp: np.mean(list(within_auc(yp, svm_decision(folds, yp, n), grp).values())), y, strata)
    summarise("primary: SVM, SAL u SMN LSS betas, participant-centred", "primary", auc, null,
              dict(n_epochs=n, n_participants_train=len(np.unique(grp)), n_aware=int(y.sum()), n_unaware=int(n - y.sum())))
    np.save(os.path.join(OD, "null_primary.npy"), null)
    # C-conf
    conf = np.c_[Ea.salience.map({"high": 0.5, "low": -0.5}).values, Ea.change_direction.map({"Faster": 0.5, "Slower": -0.5}).values]
    assert np.isfinite(conf).all()
    folds_c = fold_kernels(Xa, grp, conf)
    dc = svm_decision(folds_c, y, n)
    nullc = run_perms(lambda yp: np.mean(list(within_auc(yp, svm_decision(folds_c, yp, n), grp).values())), y, strata)
    summarise("control: salience and direction regressed out in each training fold", "control", within_auc(y, dc, grp), nullc)
    # C-phys
    Fall = physio_features(E)
    ok = Fall.notna().all(1).values & E.condition.isin(q.CHANGE_ONLY).values
    mu = Fall[ok].groupby(E.subject[ok]).mean()          # label-free: all change epochs with complete features
    F = Fall[A].reset_index(drop=True)
    cc = F.notna().all(1).values
    Fc = (F[cc].values - mu.loc[Ea.subject[cc]].values)
    yp_, gp_, sp_ = y[cc], grp[cc], strata[cc]
    dp = logit_decision(Fc, yp_, gp_)
    nullp = run_perms(lambda yy: np.mean(list(within_auc(yy, logit_decision(Fc, yy, gp_), gp_).values())), yp_, sp_)
    ap = within_auc(yp_, dp, gp_)
    ab = within_auc(yp_, d[cc], gp_)
    common = sorted(set(ap) & set(ab))
    diff = np.array([ab[s] - ap[s] for s in common]); dlo, dhi = boot_ci(diff)
    summarise("control: physiology-only benchmark (MAE, ratio deviation, PLV, RVT), logistic", "control (benchmark)", ap, nullp,
              dict(n_epochs=int(cc.sum()), n_participants_train=len(np.unique(gp_)),
                   brain_auc_same_epochs=float(np.mean([ab[s] for s in common])),
                   brain_minus_phys=float(diff.mean()), brain_minus_phys_ci95_lo=dlo, brain_minus_phys_ci95_hi=dhi))
    # benchmark coefficients on all complete cases (descriptive: which feature carries it)
    from sklearn.linear_model import LogisticRegression
    sd = Fc.std(0, ddof=1)
    coef = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000).fit(Fc / sd, yp_).coef_[0]
    feat = pd.DataFrame(dict(analysis=f"EXPLORATORY {TAG} (descriptive)", feature=["plv", "mae", "ratio_dev", "rvt"],
                             std_coef_full_fit=coef,
                             within_auc_single_feature=[np.mean(list(within_auc(yp_, Fc[:, j], gp_).values())) for j in range(4)],
                             n_epochs=int(cc.sum())))
    feat.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_physio_features.csv"), index=False)
    print(feat.round(3).to_string(index=False))
    # SECONDARY: confidence
    cv = E.groupby("subject").confidence.transform("std")
    kc = (cv > 0).values
    Ec = E[kc].reset_index(drop=True); Xc = X[kc]
    zc = (Ec.confidence - Ec.groupby("subject").confidence.transform("mean")) / Ec.groupby("subject").confidence.transform("std")
    zc = zc.values; gc = Ec.subject.values; sc = (Ec.subject + "_" + Ec.run.astype(str)).values
    rng = np.random.default_rng(SEED)
    Y = np.column_stack([zc] + [perm_labels(zc, sc, rng) for _ in range(N_PERM)])
    Pr = kridge_all(Xc, Y, gc)
    vals = np.array([np.mean(list(within_spearman(Y[:, j], Pr[:, j], gc).values())) for j in range(Y.shape[1])])
    rho = within_spearman(zc, Pr[:, 0], gc); v = np.array(list(rho.values())); lo, hi = boot_ci(v)
    obs = float(vals[0]); nl = vals[1:]
    p = float((1 + (nl >= obs).sum()) / (1 + len(nl)))
    rows.append(dict(analysis=f"EXPLORATORY {TAG}", test="secondary: kernel ridge of within-participant z confidence",
                     role="secondary", n_participants_scored=len(v), n_epochs=len(zc), mean_within_spearman=obs,
                     ci95_lo=lo, ci95_hi=hi, median_within_spearman=float(np.median(v)), n_above_0=int((v > 0).sum()),
                     null_mean=float(nl.mean()), null_sd=float(nl.std()), null_95th=float(np.percentile(nl, 95)),
                     p_one_sided=p, p_holm=p, holm_m=1, n_perm=len(nl)))
    print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_decoder.csv"), index=False)
    pp = pd.DataFrame([dict(subject=s, primary_auc=auc.get(s), confreg_auc=within_auc(y, dc, grp).get(s),
                            physio_auc=ap.get(s), brain_auc_physio_epochs=ab.get(s), conf_rho=rho.get(s),
                            n_aware=int(y[grp == s].sum()), n_unaware=int((1 - y[grp == s]).sum()))
                       for s in sorted(np.unique(E.subject))])
    pp.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_participants.csv"), index=False)
    Ea.assign(decision_primary=d, decision_confreg=dc).drop(columns=["onset", "duration"]) \
      .to_csv(os.path.join(OD, "decision_values.csv"), index=False)
    np.save(os.path.join(OD, "null_confreg.npy"), nullc); np.save(os.path.join(OD, "null_physio.npy"), nullp)
    np.save(os.path.join(OD, "null_confidence.npy"), nl)


if __name__ == "__main__":
    {"counts": counts, "run": run}[sys.argv[1]]()
