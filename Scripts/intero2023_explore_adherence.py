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


# =========================================================================== Steps 1-3
# Run after the PI's post-hoc override of the failed gate (deviations.md, 2026-11-28), exactly as the E12 row
# specifies; the open details are fixed in the override row.
H1 = ["caren_PIns", "caren_COper", "caren_SM"]
PM_ADH, PM_SR, PM_DIR = "paced_adh", "paced_steprampPM", "paced_dirPM"
MODELS = {  # name -> (measure, collapse change epochs, physio columns, variant label for file names)
    "E12adh": ("plv", False, [], "base"),            # primary
    "E12collapsed": ("plv", True, [], "base"),       # sensitivity (a): no outcome split
    "E12adhrvt": ("plv", False, g.PHYSIO_RVT, "rvt"),  # sensitivity (b): + RVT x RRF, S2b
    "E12mae": ("mae", False, [], "base"),            # sensitivity (c): MAE as the modulator
}


def run_events(sub, ent, n, physio, P):
    ev, cf, _ = g.read_run_physio(sub, ent, n, physio)
    ev = ev.reset_index(drop=True)
    ev["epoch_index"] = pd.to_numeric(ev.epoch_index, errors="coerce")
    m = P[(P.subject == sub) & (P.bids_entity == ent)].set_index("epoch_index")
    pk = ev.condition.isin(q.PACED_MOD)
    for c in ("plv", "mae", "change_type", "direction"):
        ev[c] = np.where(pk, ev.epoch_index.map(m[c]), np.nan) if c in ("plv", "mae") else \
            ev.epoch_index.map(m[c]).where(pk)
    assert ev.loc[pk, "change_type"].notna().all()
    return ev, cf


def e12_design(ev, cf, n, measure, collapse):
    base = ev[["onset", "duration", "condition"]].rename(columns={"condition": "trial_type"}).assign(modulation=1.0)
    if collapse:
        base.loc[base.trial_type.isin(q.CHANGE_ONLY), "trial_type"] = "paced_change"
    parts = [base]
    a = ev[ev.condition.isin(q.PACED_MOD) & ev[measure].notna()]
    if len(a) >= 2 and a[measure].std() > 0:
        parts.append(a[["onset", "duration"]].assign(trial_type=PM_ADH, modulation=(a[measure] - a[measure].mean()).values))
    ch = ev[ev.condition.isin(q.CHANGE_ONLY)]
    for name, val in ((PM_SR, ch.change_type.map({"step": 0.5, "ramp": -0.5})),
                      (PM_DIR, ch.direction.map({"Faster": 0.5, "Slower": -0.5}))):
        assert val.notna().all()
        if len(ch) >= 2 and val.std() > 0:
            parts.append(ch[["onset", "duration"]].assign(trial_type=name, modulation=(val - val.mean()).values))
    return q.design(pd.concat(parts, ignore_index=True), cf, n)


def joint_designs(sub, runs, P, model):
    measure, collapse, physio, _ = MODELS[model]
    dms = []
    for r in runs.itertuples():
        ev, cf = run_events(sub, r.bids_entity, int(r.n_vol), physio, P)
        dms.append(e12_design(ev, cf, int(r.n_vol), measure, collapse))
    return dms, hj.build_joint_design(dms)


def eligible(sub, runs, P, model, s2b):
    measure, _, physio, _ = MODELS[model]
    if physio and sub not in s2b:
        return False
    return P[(P.subject == sub) & P[measure].notna()].shape[0] >= 2


def qc():
    from intero2023_explore_fir import vif, rank_def
    import intero2023_group as G
    P = epoch_table(); man = g.run_list(); s2b = set(G.subject_sets()[0]["S2b"]); rows = []
    for sub, runs in man.groupby("subject"):
        for model in MODELS:
            if not eligible(sub, runs, P, model, s2b):
                continue
            dms, J = joint_designs(sub, runs, P, model)
            v = vif(J, [PM_ADH, PM_SR, PM_DIR])
            Xc = J[[c for c in (PM_ADH, PM_SR, PM_DIR) if c in J]]
            rows.append(dict(subject=sub, model=model, n_runs=len(dms), n_reg=J.shape[1], rank_deficient=int(rank_def(J)),
                             n_adh_epochs=int(P[(P.subject == sub) & P[MODELS[model][0]].notna()].shape[0]),
                             var_adh=hj.cvar(J, {PM_ADH: 1.0}),
                             vif_adh=v.get(PM_ADH), vif_steprampPM=v.get(PM_SR), vif_dirPM=v.get(PM_DIR),
                             r_adh_stepramp=float(np.corrcoef(Xc[PM_ADH], Xc[PM_SR])[0, 1]) if PM_SR in Xc else np.nan,
                             r_adh_dir=float(np.corrcoef(Xc[PM_ADH], Xc[PM_DIR])[0, 1]) if PM_DIR in Xc else np.nan))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_joint.csv"), index=False)
    for m, u in d.groupby("model", sort=False):
        ok = u.var_adh.notna()
        print(f"{m:13s} N {ok.sum():2d} (eligible {len(u)}); rank-deficient {u.rank_deficient.sum()}; "
              f"var median {u.var_adh.median():.3f} [{u.var_adh.min():.3f}, {u.var_adh.max():.3f}]; "
              f"VIF adh median {u.vif_adh.median():.2f} max {u.vif_adh.max():.2f}; "
              f"VIF stepramp {u.vif_steprampPM.median():.2f}, dir {u.vif_dirPM.median():.2f}; "
              f"|r| adh-stepramp median {u.r_adh_stepramp.abs().median():.2f}, adh-dir {u.r_adh_dir.abs().median():.2f}")


def subject(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    import intero2023_group as G
    P = epoch_table(); man = g.run_list(); runs = man[man.subject == sub]; assert len(runs), sub
    s2b = set(G.subject_sets()[0]["S2b"])
    todo = [m for m in MODELS if eligible(sub, runs, P, m, s2b)]
    od = os.path.join(g.OUT, "firstlevel", sub); os.makedirs(od, exist_ok=True)
    if not todo:
        print(sub, "not eligible"); return
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    Ys = []
    for r in runs.itertuples():
        Y = masker.transform(g.bold_path(sub, r.bids_entity)); assert Y.shape[0] == int(r.n_vol)
        Ys.append(mean_scaling(Y, axis=0)[0])
    Y = np.vstack(Ys); log = dict(subject=sub, models={})
    for m in todo:
        _, J = joint_designs(sub, runs, P, m)
        if not np.isfinite(hj.cvar(J, {PM_ADH: 1.0})):
            log["models"][m] = "not estimable"; continue
        out, dof = fit(Y, J, {m: {PM_ADH: 1.0}})
        e, v = out[m]
        if not (v > 0).all():
            raise ValueError(f"{sub} {m}: non-positive variances")
        var = MODELS[m][3]
        for arr, kind in ((e, "effect"), (v, "variance")):
            g.save_float32(masker.inverse_transform(arr), os.path.join(od, f"{sub}_var-{var}_con-{m}joint_{kind}.nii.gz"))
        log["models"][m] = dict(dof=dof, n_reg=J.shape[1])
    json.dump(log, open(os.path.join(od, f"{sub}_e12.json"), "w"))
    print(sub, "done", log["models"], flush=True)


def group():
    import nibabel as nib
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.second_level import non_parametric_inference
    import intero2023_group as G
    from intero2023_explore_awareness import jzs_bf10
    from intero2023_explore_rating import atlas_labels
    man = g.run_list(); allsubs = sorted(man.subject.unique()); rows = []
    roles = {"E12adh": ("primary", True), "E12collapsed": ("sensitivity (a): no outcome split", False),
             "E12adhrvt": ("sensitivity (b): + RVT x RRF, S2b", False), "E12mae": ("sensitivity (c): MAE modulator", False)}
    for name, (role, wb) in roles.items():
        var = MODELS[name][3]
        subs = [s for s in allsubs if os.path.exists(G.fl(s, var, f"{name}joint", "effect"))]
        files = [G.fl(s, var, f"{name}joint", "effect") for s in subs]
        imgs = [nib.load(f) for f in files]
        dm = pd.DataFrame({"intercept": np.ones(len(files))})
        json.dump(dict(N=len(subs), subjects=subs), open(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_subjects.json"), "w"))
        for roi in H1:
            rmask = nib.load(G.ROI(roi)).get_fdata() > 0
            lp = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="intercept", mask=G.ROI(roi),
                                          n_perm=G.N_PERM, two_sided_test=True, random_state=SEED, n_jobs=G.N_JOBS)
            tmap = G.ols_t(imgs, dm, "intercept", G.ROI(roi))
            lpv = lp.get_fdata()[rmask]; tv = tmap.get_fdata()[rmask]
            k = int(np.nanargmax(np.abs(tv))); ijk = np.argwhere(rmask)[k]
            pmin = float(10 ** -np.nanmax(lpv))
            mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
            y = np.array([mk.transform(f).mean() for f in files])
            row = dict(analysis=f"EXPLORATORY {TAG} (post-hoc override of the gate)", contrast=name, role=role, roi=roi,
                       N=len(files), peak_abs_t_signed=round(float(tv[k]), 3),
                       peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(tmap.affine, ijk)],
                       p_fwe=pmin, p_fwe_bonf3=min(1.0, 3 * pmin), n_vox_p_lt_05=int((lpv > -np.log10(0.05)).sum()),
                       mask_mean=float(y.mean()), mask_mean_dz=float(y.mean() / y.std(ddof=1)),
                       mask_mean_t=float(stats.ttest_1samp(y, 0).statistic))
            for r in (0.5, 0.707, 1.0):
                row[f"BF10_r{r}"] = jzs_bf10(row["mask_mean_t"], len(y), r)
            row["BF01_r0.707"] = 1 / row["BF10_r0.707"]
            rows.append(row); print({k_: v_ for k_, v_ in row.items()}, flush=True)
        if wb:
            od = os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}_{name}"); os.makedirs(od, exist_ok=True)
            res = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept", mask=g.GROUP_MASK,
                                           n_perm=G.N_PERM, two_sided_test=True, threshold=0.002,
                                           random_state=SEED, n_jobs=G.N_JOBS)
            for kk in ("t", "logp_max_size"):
                g.save_float32(res[kk], os.path.join(od, f"wholebrain_{kk}.nii.gz"))
            cl = G.cluster_table(res["t"], res["logp_max_size"], len(files) - 1)
            if len(cl):
                cl = cl.assign(labels=[atlas_labels(m) for m in cl.peak_mni])
            cl.insert(0, "contrast", name)
            cl.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_wholebrain.csv"), index=False)
            print(name, "whole brain:", cl.head(15).to_string() if len(cl) else "no cluster", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_h1masks.csv"), index=False)


# --------------------------------------------------------------------------- Step 2
def ridge():
    """LOSO linear kernel ridge on the T4 LSS betas: z(PLV) within participant x run, residualised on
    [1, No Change, step-ramp, direction] in each training fold; mean within-participant Spearman; one-sided."""
    sys.path.insert(0, HERE)
    from intero2023_explore_awareness_decoder import perm_labels, within_spearman, boot_ci, T4D
    P = epoch_table(); P = P[P.plv.notna()]
    X, rows = [], []
    for s in sorted(P.subject.unique()):
        f = os.path.join(T4D, f"{s}_lss.npz")
        if not os.path.exists(f):
            continue
        B = np.load(f)["beta"].astype(np.float64)
        E = pd.read_csv(os.path.join(T4D, f"{s}_lss_epochs.csv")); E["epoch_index"] = E.epoch_index.astype(int)
        E = E.reset_index().merge(P[P.subject == s][["bids_entity", "epoch_index", "plv", "change_type", "direction"]]
                                  .rename(columns={"bids_entity": "entity"}), on=["entity", "epoch_index"])
        E = E[E.condition.isin(q.PACED_MOD)]
        for rn, u in E.groupby("run"):
            if len(u) < 3 or u.plv.std() == 0:
                continue
            Bx = B[u["index"].values]; X.append(Bx - Bx.mean(0, keepdims=True))
            rows.append(u.assign(z=(u.plv - u.plv.mean()) / u.plv.std()))
    E = pd.concat(rows, ignore_index=True); X = np.vstack(X)
    grp = E.subject.values; strata = (E.subject + "_" + E.run.astype(str)).values
    C = np.column_stack([np.ones(len(E)), (E.change_type == "none").astype(float),
                         np.select([E.change_type == "step", E.change_type == "ramp"], [0.5, -0.5], 0.0),
                         np.select([E.direction == "Faster", E.direction == "Slower"], [0.5, -0.5], 0.0)])
    rng = np.random.default_rng(SEED)
    Y = np.column_stack([E.z.values] + [perm_labels(E.z.values, strata, rng) for _ in range(5000)])
    Pr = np.full(Y.shape, np.nan); Yo = np.full(Y.shape, np.nan)
    for s in np.unique(grp):
        te = np.flatnonzero(grp == s); tr = np.flatnonzero(grp != s)
        Bc = np.linalg.lstsq(C[tr], Y[tr], rcond=None)[0]
        Ytr = Y[tr] - C[tr] @ Bc; Yo[te] = Y[te] - C[te] @ Bc
        K = X[tr] @ X[tr].T; lam = float(np.mean(np.diag(K))); w, V = np.linalg.eigh(K)
        Pr[te] = (X[te] @ X[tr].T) @ (V @ ((V.T @ Ytr) / (w + lam)[:, None]))
    vals = np.array([np.mean(list(within_spearman(Yo[:, j], Pr[:, j], grp).values())) for j in range(Y.shape[1])])
    rho = np.array(list(within_spearman(Yo[:, 0], Pr[:, 0], grp).values())); lo, hi = boot_ci(rho)
    out = pd.DataFrame([dict(analysis=f"EXPLORATORY {TAG} Step 2 (post-hoc override of the gate)",
                             test="LOSO kernel ridge, T4 LSS betas -> z(PLV) net of change type and direction",
                             n_participants_scored=len(rho), n_epochs=len(E), mean_within_spearman=float(vals[0]),
                             ci95_lo=lo, ci95_hi=hi, median_within_spearman=float(np.median(rho)), n_above_0=int((rho > 0).sum()),
                             null_mean=float(vals[1:].mean()), null_95th=float(np.percentile(vals[1:], 95)),
                             p_one_sided=float((1 + (vals[1:] >= vals[0]).sum()) / len(vals)), n_perm=len(vals) - 1)])
    out.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_ridge.csv"), index=False)
    print(out.T.to_string())


if __name__ == "__main__":
    a = sys.argv[1]
    {"gate": gate, "qc": qc, "subject": lambda: subject(sys.argv[2]), "group": group, "ridge": ridge}[a]()
