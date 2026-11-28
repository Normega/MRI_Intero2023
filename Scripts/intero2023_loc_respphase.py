#!/usr/bin/env python
"""
intero2023_loc_respphase.py -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §10.6):
inhale- vs exhale-locked responses by condition, in the physio subset, with per-participant
belt polarity. Everything this script writes is exploratory and descriptive.

Choices: Analysis/localizer_stageC3_implementation_notes.md §6 (committed before any value).

  python Scripts/intero2023_loc_respphase.py segments       # raw belt -> per-run breath tables (cache)
  python Scripts/intero2023_loc_respphase.py check          # design only: collinearity of the two models
  python Scripts/intero2023_loc_respphase.py subject <sub>  # per-run fits in the a-priori masks
  python Scripts/intero2023_loc_respphase.py group          # tables

Belt, polarity and breath detection are the step-1 / stage-C1 functions of
intero2023_localizer_qc (physio_segments, bandpass, breath_table): polarity +1 means the
belt rises on inhalation, inhale onset = trough, exhale onset = peak. Phase = angle of the
Hilbert transform of the polarity-corrected band-passed belt (0 at the end of inhalation),
as stage C1. Volume k is sampled at k * TR + StartTime (fMRIPrep slice-timing reference).

Two models, both added to the registered per-run localizer design (intero2023_loc_common.design,
6 mm, 200 s DCT, swc7p nuisance sheet, AR(1)):
  onset model  (HRF)  per condition c: impulse events at the inhale onsets and at the exhale
                      onsets inside c's task window, convolved with the SPM canonical HRF
                      (inh_c, exh_c); contrast inh_c - exh_c.
  phase model  (none) per condition c: cos(phi) and sin(phi) at each volume inside c's task
                      window (0 elsewhere), not convolved (RETROICOR-like, first order);
                      inhale - exhale = mean over the inhale half-cycle minus the exhale half =
                      -(4 / pi) * b_sin; amplitude = sqrt(b_cos^2 + b_sin^2).
Imaging outputs: $OUT/loc_c3/respphase/<sub>.json (mask values only). Tables:
Analysis/results/EXPLORATORY_LOC_C3_respphase*.csv.
"""
import os, sys, json, pickle, warnings
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
OD = os.path.join(g.OUT, "loc_c3", "respphase")
CACHE = os.path.join(OD, "belt_segments.pkl")
CONDS = ["intero", "extero", "match", "self", "rest"]
MASKS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "caren_COper", "caren_SM", "ieat_ACCsphere"]
START_TIME = 0.958          # fMRIPrep slice-timing reference (s after the volume trigger)
SKIP = 3.0                  # first 3 s of every block dropped (stage C1 rule)
VIF_MAX = 10.0              # pre-set feasibility rule for the onset model (notes §6)


def bold_json(sub, run):
    return os.path.join(g.FMRIPREP, sub, "func",
                        f"{sub}_{L.entity(run)}_{g.SPACE}_desc-preproc_bold.json")


# --------------------------------------------------------------------------- belt
def segments():
    import intero2023_localizer_qc as qc
    from scipy.signal import hilbert
    x, blocks, presses, U, seg, POL, polarity = qc.physio_segments()
    runs = L.run_list()
    out = {}
    Ui = U.set_index(["subject", "run"])
    for r in runs.itertuples():
        k = (r.subject, int(r.run))
        if k not in seg:
            continue
        u = Ui.loc[k]
        pol = polarity.get(r.subject, 0)
        ok = (u.belt_grade in ("good", "degraded")) and (u.grade_time in ("exact", "good")) \
            and int(u.pulse_delta) == 0 and pol != 0
        if not ok:
            continue
        sg = seg[k]; fs = sg["fs"]
        xf = qc.bandpass(sg["belt"], fs)
        bt = qc.breath_table(xf, fs, pol)
        ph = np.angle(hilbert(xf * pol))
        out[k] = dict(fs=fs, breaths=bt, phase=ph.astype(np.float32), polarity=pol,
                      grade_time=u.grade_time, belt_grade=u.belt_grade)
    os.makedirs(OD, exist_ok=True)
    pickle.dump(dict(runs=out), open(CACHE, "wb"))
    print(f"segments: {len(out)} usable runs, {len({s for s, _ in out})} participants")


def task_windows(sub, run):
    """Task windows per condition from the events file (Rest included; Intero = compliant blocks)."""
    ev = pd.read_csv(L.events_path(sub, run), sep="\t")
    it = ev.trial_type == "intero"
    nc = it & (pd.to_numeric(ev.n_presses, errors="coerce") < L.MIN_PRESSES)
    ev.loc[nc, "trial_type"] = "intero_nc"
    return {c: list(zip(ev.onset[ev.trial_type == c], (ev.onset + ev.duration)[ev.trial_type == c]))
            for c in CONDS}


def extra_regressors(sub, run, n_vol, info):
    """Onset-model and phase-model columns for one run (both as DataFrames, n_vol rows)."""
    from nilearn.glm.first_level import compute_regressor
    st = START_TIME
    if os.path.exists(bold_json(sub, run)):
        st = float(json.load(open(bold_json(sub, run))).get("StartTime", START_TIME))
    ft = TRk = L.TR * np.arange(n_vol) + st
    win = task_windows(sub, run)
    bt = info["breaths"]; ph = info["phase"]; fs = info["fs"]
    on, phm, counts = {}, {}, {}
    ii = np.clip(np.round(ft * fs).astype(int), 0, len(ph) - 1)
    phv = ph[ii]
    for c in CONDS:
        tin, tex = [], []
        inwin = np.zeros(n_vol, bool)
        for a, b in win[c]:
            a2 = a + SKIP
            tin += [t for t in bt.t_in if a2 <= t < b]
            tex += [t for t in bt.t_ex if a2 <= t < b]
            inwin |= (ft >= a2) & (ft < b)
        counts[c] = dict(n_inhale=len(tin), n_exhale=len(tex), n_vol=int(inwin.sum()))
        for lab, ts in (("inh", tin), ("exh", tex)):
            if len(ts):
                cond = np.vstack([ts, np.zeros(len(ts)), np.ones(len(ts))])
                on[f"{lab}_{c}"] = compute_regressor(cond, "spm", TR_frame(n_vol))[0].ravel()
        if inwin.sum() >= 10:
            phm[f"cos_{c}"] = np.where(inwin, np.cos(phv), 0.0)
            phm[f"sin_{c}"] = np.where(inwin, np.sin(phv), 0.0)
    return pd.DataFrame(on), pd.DataFrame(phm), counts


def TR_frame(n_vol):
    return L.TR * (np.arange(n_vol) + 0.5)          # the registered design's frame times


def base_design(sub, run, n_vol):
    ev, _ = L.read_events(sub, run, n_vol)
    cf = L.read_confounds(sub, run, n_vol)
    return L.design(ev, cf, n_vol)


def vif(X, j):
    y = X[:, j]; Z = np.delete(X, j, axis=1)
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    r2 = 1 - ((y - Z @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return 1 / max(1 - r2, 1e-12)


def check():
    """Design only (no BOLD): VIF of each added column and r(inh_c, exh_c) in every usable run."""
    S = pickle.load(open(CACHE, "rb"))["runs"]
    runs = L.run_list().set_index(["subject", "run"])
    rows = []
    for (sub, run), info in sorted(S.items()):
        n = int(runs.loc[(sub, run), "n_vol"])
        dm = base_design(sub, run, n)
        on, phm, counts = extra_regressors(sub, run, n, info)
        for model, extra in (("onset", on), ("phase", phm)):
            X = pd.concat([dm.reset_index(drop=True), extra], axis=1)
            Xs = X.values.astype(float)
            sd = Xs.std(0); sd[sd < 1e-12] = 1; Xs = Xs / sd
            for c in CONDS:
                cols = [f"inh_{c}", f"exh_{c}"] if model == "onset" else [f"cos_{c}", f"sin_{c}"]
                if not all(k in X.columns for k in cols):
                    continue
                j = [list(X.columns).index(k) for k in cols]
                rows.append(dict(subject=sub, run=run, model=model, condition=c, **counts[c],
                                 vif_a=vif(Xs, j[0]), vif_b=vif(Xs, j[1]),
                                 r_ab=float(np.corrcoef(X[cols[0]], X[cols[1]])[0, 1])))
    D = pd.DataFrame(rows)
    D.round(4).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_respphase_design_check.csv"), index=False)
    s = D.groupby(["model", "condition"]).agg(n_runs=("run", "size"), vif_a_median=("vif_a", "median"),
                                             vif_b_median=("vif_b", "median"), r_ab_median=("r_ab", "median"),
                                             n_inhale_median=("n_inhale", "median"))
    print(s.round(3).to_string())
    print("participants:", D.subject.nunique(), "runs:", D.groupby(["subject", "run"]).ngroups)


# --------------------------------------------------------------------------- fits
def subject(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    import intero2023_group as G
    import nibabel as nib
    S = pickle.load(open(CACHE, "rb"))["runs"]
    runs = L.run_list(); runs = runs[runs.subject == sub]
    model_feasible = json.load(open(os.path.join(OD, "feasibility.json")))
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    mvec = {m: (nib.load(G.ROI(m)).get_fdata() > 0)[gm] for m in MASKS}
    out = dict(subject=sub, runs={})
    for r in runs.itertuples():
        k = (sub, int(r.run))
        if k not in S:
            continue
        n = int(r.n_vol)
        Y = masker.transform(g.bold_path(sub, L.entity(r.run)))
        Y, _ = mean_scaling(Y, axis=0)
        Ym = np.column_stack([Y[:, mvec[m]].mean(1) for m in MASKS])       # mask-mean PSC series
        dm = base_design(sub, int(r.run), n)
        on, phm, counts = extra_regressors(sub, int(r.run), n, S[k])
        rec = dict(counts=counts, models={})
        for model, extra in (("onset", on), ("phase", phm)):
            if not model_feasible.get(model, False):
                continue
            X = pd.concat([dm.reset_index(drop=True), extra], axis=1)
            specs = {}
            for c in CONDS:
                if model == "onset" and {f"inh_{c}", f"exh_{c}"} <= set(X.columns):
                    specs[f"ie_{c}"] = {f"inh_{c}": 1.0, f"exh_{c}": -1.0}
                if model == "phase" and {f"cos_{c}", f"sin_{c}"} <= set(X.columns):
                    specs[f"cos_{c}"] = {f"cos_{c}": 1.0}
                    specs[f"sin_{c}"] = {f"sin_{c}": 1.0}
            res, dof = fit(Ym, X, specs)
            rec["models"][model] = {kk: dict(effect=e.tolist(), variance=v.tolist()) for kk, (e, v) in res.items()}
            rec["models"][model]["dof"] = dof
        out["runs"][int(r.run)] = rec
    out["masks"] = MASKS
    json.dump(out, open(os.path.join(OD, f"{sub}.json"), "w"))
    print(f"{sub}: {len(out['runs'])} runs", flush=True)


# --------------------------------------------------------------------------- group
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def group():
    s, _, _ = L.participant_sets()
    feas = json.load(open(os.path.join(OD, "feasibility.json")))
    rows = []
    for f in sorted(os.listdir(OD)):
        if not (f.startswith("sub-") and f.endswith(".json")):
            continue
        d = json.load(open(os.path.join(OD, f)))
        for run, rec in d["runs"].items():
            for model, res in rec["models"].items():
                for k, v in res.items():
                    if k == "dof":
                        continue
                    kind, c = k.split("_", 1)
                    for m, e, var in zip(d["masks"], v["effect"], v["variance"]):
                        rows.append(dict(subject=d["subject"], run=int(run), model=model, term=kind,
                                         condition=c, mask=m, effect=e, variance=var))
    R = pd.DataFrame(rows)
    # Intero requires the participant to be in the Intero set (compliance rule, as every Intero contrast)
    R = R[(R.condition != "intero") | R.subject.isin(s["intero"])]
    # participant value = inverse-variance FFX over runs
    def ffx(d):
        w = 1 / d.variance
        return pd.Series(dict(effect=(w * d.effect).sum() / w.sum(), n_runs=len(d)))
    P = R.groupby(["subject", "model", "term", "condition", "mask"]).apply(ffx).reset_index()
    P = P.pivot_table(index=["subject", "model", "condition", "mask"], columns="term", values="effect").reset_index()
    if "sin" in P:
        P["ie"] = np.where(P.model == "phase", -(4 / np.pi) * P["sin"], P.get("ie"))
        P["amp"] = np.where(P.model == "phase", np.hypot(P["cos"], P["sin"]), np.nan)
    P.round(6).to_csv(os.path.join(OD, "participant_values.csv"), index=False)
    out = []
    for (model, mask), d in P.groupby(["model", "mask"]):
        W = d.pivot_table(index="subject", columns="condition", values="ie")
        for c in CONDS:
            if c not in W:
                continue
            y = W[c].dropna()
            t = stats.ttest_1samp(y, 0)
            out.append(dict(model=model, mask=mask, test=f"inhale-exhale {c}", N=len(y), mean=y.mean(),
                            sd=y.std(ddof=1), dz=y.mean() / y.std(ddof=1), t=t.statistic, p=t.pvalue))
        for a, b in (("intero", "extero"), ("match", "extero"), ("intero", "rest"), ("match", "intero")):
            if a in W and b in W:
                y = (W[a] - W[b]).dropna(); t = stats.ttest_1samp(y, 0)
                out.append(dict(model=model, mask=mask, test=f"inhale-exhale {a} - {b}", N=len(y), mean=y.mean(),
                                sd=y.std(ddof=1), dz=y.mean() / y.std(ddof=1), t=t.statistic, p=t.pvalue))
        if model == "phase":
            A = d.pivot_table(index="subject", columns="condition", values="amp")
            for a, b in (("intero", "extero"), ("match", "extero"), ("intero", "rest")):
                y = (A[a] - A[b]).dropna(); t = stats.ttest_1samp(y, 0)
                out.append(dict(model=model, mask=mask, test=f"amplitude {a} - {b}", N=len(y), mean=y.mean(),
                                sd=y.std(ddof=1), dz=y.mean() / y.std(ddof=1), t=t.statistic, p=t.pvalue))
    O = pd.DataFrame(out)
    O["family"] = np.where(O.test.str.startswith("amplitude"), "amplitude", np.where(O.test.str.contains(" - "), "ie_condition_difference", "ie_condition"))
    O["q_fdr"] = np.nan
    for (model, fam), ix in O.groupby(["model", "family"]).groups.items():
        O.loc[ix, "q_fdr"] = bh(O.loc[ix, "p"].values)
    O.round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_respphase.csv"), index=False)
    print(O[O.q_fdr < .05].round(3).to_string())
    print(O.groupby(["model", "family"]).p.apply(lambda p: (p < .05).sum()))


if __name__ == "__main__":
    if sys.argv[1] == "subject":
        subject(sys.argv[2])
    else:
        {"segments": segments, "check": check, "group": group}[sys.argv[1]]()
