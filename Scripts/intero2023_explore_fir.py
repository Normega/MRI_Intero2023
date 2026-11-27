#!/usr/bin/env python
"""
intero2023_explore_fir.py -- EXPLORATORY E7 (deviations.md, "After the update", 2026-11-26):
is there a transient BOLD response to the pacer change, in the three H3 salience masks?

  python Scripts/intero2023_explore_fir.py qc               # design only, no BOLD
  python Scripts/intero2023_explore_fir.py subject sub-057  # joint two-run FIR first level
  python Scripts/intero2023_explore_fir.py group            # logged group tests + time courses

Change type. In every pacer trace, salience 'high' is a STEP (one abrupt ~20% change of the
half-cycle period at pulse 5, ~9.5 s after epoch onset; later (11.5 s, or pulse 7 at 15.9 s) in the long v1
epochs) and salience 'low' is a RAMP (the period starts drifting at pulse 3, ~4.5 s, by ~7%
per pulse). `qc` re-derives both from Data/events/*_pacertrace.tsv and checks the mapping.

Design (per run). Model A's non-paced regressors (free, rating, subrating; 'spm + derivative')
and its nuisance columns, unchanged (built by twomodel_qc.design on the same events). The
paced epochs are modelled by FIR instead of the convolved boxcars: 15 bins of 2 s (0-30 s
from epoch onset), separately for step, ramp (aware + unaware + wrong-direction, no awareness
split) and No Change. Bin b is 1 at the volume whose acquisition time (TR*(k+0.5)) falls in
[onset + 2b, onset + 2b + 2). Epochs longer than 25 s (16 v1 epochs in 4 runs: sub-052 both
runs, sub-053 both runs) have different change times, so they do not enter the FIR; they
are modelled by one shared convolved boxcar 'paced_long' (+ derivative), of no interest.

First level: joint two-run model as for H3/E3 (build_joint_design: shared task columns incl.
the FIR bins, run-specific nuisance and intercept), PSC per run, AR(1), 6 mm on the group mask.

Contrasts (windows in s from epoch onset; each window = mean of its 3 bins):
  early W1 = [6, 12): ramp change (4.5 s) + ~5 s HRF peak; before the step response.
  late  W2 = [12, 18): step change (9.5 s) + ~5 s HRF peak; ends before the earliest
                       rating-onset response (shortest epoch 17.5 s + ~2 s).
  PRIMARY  E7timing  = (step - ramp)[W2] - (step - ramp)[W1]   (change-locked timing shift)
  secondary E7stepNC = (step - No Change)[W2]
  secondary E7rampNC = (ramp - No Change)[W1]
Group: sign-flip permutation voxel max-t FWE in each H3 mask, two-sided, 10,000 permutations,
seed 20261126, Bonferroni x3; mask-mean dz and JZS BF10 (r = 0.707; 0.5, 1.0 robustness);
whole-brain two-sided cluster FWE (p < .001 per tail) on the primary only.
"""
import os, re, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
import intero2023_glm_common as g
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
TYPES = ("step", "ramp", "nochange")
NBIN, BINW, LONG_S = 15, 2.0, 25.0
W1, W2 = [3, 4, 5], [6, 7, 8]            # bin indices: [6,12) s and [12,18) s
CHANGE = q.CHANGE_ONLY
TAG = "E7"


def col(t, b):
    return f"paced_{t}_fir{b:02d}"


def epoch_type(row):
    if row.condition == "paced_nochange":
        return "nochange"
    return {"high": "step", "low": "ramp"}[row.salience]


def joint_design(sub, runs):
    """Per-run FIR designs, joined as for H3 (build_joint_design), censored FIR bins dropped."""
    dms, cnts = [], []
    for r in runs.itertuples():
        dm, cnt = fir_design(sub, r.bids_entity, int(r.n_vol)); dms.append(dm); cnts.append(cnt)
    J, gone = drop_censored(hj.build_joint_design(dms))
    return J, dms, cnts, gone


def fir_design(sub, ent, n_vol):
    """Model A with the paced boxcars replaced by FIR (step / ramp / No Change) + paced_long."""
    ev, cf, _ = g.read_run_physio(sub, ent, n_vol, [])
    paced = ev[ev.condition.isin(q.PACED_MOD)].copy()
    paced["etype"] = [epoch_type(r) for r in paced.itertuples()]
    long_ = paced[paced.duration > LONG_S]
    fir = paced[paced.duration <= LONG_S]
    other = ev[~ev.condition.isin(q.PACED_MOD)][["onset", "duration", "condition"]] \
        .rename(columns={"condition": "trial_type"})
    if len(long_):
        other = pd.concat([other, long_[["onset", "duration"]].assign(trial_type="paced_long")])
    dm = q.design(other.assign(modulation=1.0), cf, n_vol)
    ft = g.TR * (np.arange(n_vol) + 0.5)
    cols = {}
    for t in TYPES:
        on = fir[fir.etype == t].onset.values
        if not len(on):
            continue
        for b in range(NBIN):
            x = np.zeros(n_vol)
            for o in on:
                rel = ft - o
                x[(rel >= b * BINW) & (rel < (b + 1) * BINW)] = 1.0
            cols[col(t, b)] = x
    firdf = pd.DataFrame(cols, index=dm.index)
    dm = pd.concat([firdf, dm], axis=1)
    counts = {t: int((fir.etype == t).sum()) for t in TYPES}
    counts["long"] = len(long_)
    return dm, counts


def drop_censored(dm):
    """Drop FIR bins whose every volume is censored by a spike / non-steady-state regressor
    (such a bin is not estimable: all its data points are absorbed by the censoring columns)."""
    cens = [c for c in dm.columns if re.search(r"(spike|non_steady_state)", c)]
    hit = (dm[cens].values != 0).any(1) if cens else np.zeros(len(dm), bool)
    fir = [c for c in dm.columns if "_fir" in c]
    gone = [c for c in fir if not ((dm[c].values != 0) & ~hit).any()]
    return dm.drop(columns=gone), gone


def window_specs(columns):
    """The logged contrasts, with each window averaged over the bins present for EVERY type the
    contrast involves (a bin censored for one type is dropped for all, so the window stays
    matched). A contrast with an empty window is not estimable for that participant."""
    cols = set(columns)
    def win(types, bins):
        return [b for b in bins if all(col(t, b) in cols for t in types)]
    out = {}
    w1, w2 = win(("step", "ramp"), W1), win(("step", "ramp"), W2)
    if w1 and w2:
        p = {}
        for b in w2:
            p[col("step", b)] = 1 / len(w2); p[col("ramp", b)] = -1 / len(w2)
        for b in w1:
            p[col("step", b)] = -1 / len(w1); p[col("ramp", b)] = 1 / len(w1)
        out["E7timing"] = p
    w = win(("step", "nochange"), W2)
    if w:
        out["E7stepNC"] = {**{col("step", b): 1 / len(w) for b in w}, **{col("nochange", b): -1 / len(w) for b in w}}
    w = win(("ramp", "nochange"), W1)
    if w:
        out["E7rampNC"] = {**{col("ramp", b): 1 / len(w) for b in w}, **{col("nochange", b): -1 / len(w) for b in w}}
    return out


def vif(dm, names):
    X = dm.drop(columns=[c for c in dm.columns if c == "constant"]).values.astype(float)
    keep = X.std(0) > 1e-12
    X = X[:, keep]; cn = list(np.array([c for c in dm.columns if c != "constant"])[keep])
    Z = (X - X.mean(0)) / X.std(0)
    R = np.corrcoef(Z, rowvar=False)
    v = np.diag(np.linalg.pinv(R, rcond=1e-12))
    return {n: float(v[cn.index(n)]) for n in names if n in cn}


def rank_def(dm):
    X = dm.values.astype(float); sd = X.std(0); sd[sd < 1e-12] = 1.0
    Xs = X / sd; sv = np.linalg.svd(Xs, compute_uv=False)
    return int((sv > max(Xs.shape) * np.finfo(float).eps * sv.max()).sum()) < Xs.shape[1]


def change_times():
    """Per paced epoch: the first pulse (>= 3) whose period differs from pulse 2's by > 2%."""
    man = g.run_list(); rows = []
    for r in man.itertuples():
        pt = pd.read_csv(f"{g.REPO}/Data/events/{r.subject}_{r.bids_entity}_pacertrace.tsv", sep="\t",
                         na_values="n/a")
        ev = pd.read_csv(f"{q.FL}/{r.subject}_{r.bids_entity}_desc-model_events.tsv", sep="\t", na_values="n/a")
        for e in ev[ev.condition.isin(q.PACED_MOD)].itertuples():
            p = pt[pt.epoch_index == e.epoch_index].sort_values("pulse_index")
            per = p.period_s.values; base = per[1]; rel = np.abs(per / base - 1)
            idx = [i for i in range(2, len(per)) if np.isfinite(per[i]) and rel[i] > 0.02]
            k = idx[0] if idx else None
            fin = per[np.isfinite(per)]
            rows.append(dict(
                subject=r.subject, entity=r.bids_entity, epoch=int(e.epoch_index), condition=e.condition,
                salience=e.salience if isinstance(e.salience, str) else "n/a",
                pacer_salience=p.salience.iloc[0] if isinstance(p.salience.iloc[0], str) else "n/a",
                direction=e.change_direction, duration=e.duration, long_epoch=e.duration > LONG_S,
                onset_diff=abs(e.onset - p.onset.iloc[0]),
                change_pulse=k + 1 if k is not None else np.nan,
                change_s=p.onset.iloc[k] - p.onset.iloc[0] if k is not None else np.nan,
                first_step_rel=rel[k] if k is not None else np.nan,
                frac_total_at_change=(abs(per[k] - base) / abs(fin[-1] - base)) if k is not None else np.nan,
                max_rel_nochange=float(np.nanmax(rel[2:])) if k is None else np.nan))
    return pd.DataFrame(rows)


def qc():
    os.chdir(g.REPO)
    ct = change_times()
    assert (ct.onset_diff < 0.01).all()
    assert (ct.salience == ct.pacer_salience).all()
    assert ct[ct.condition == "paced_nochange"].change_pulse.isna().all()
    assert ct[ct.condition != "paced_nochange"].change_pulse.notna().all()
    ct.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_change_times.csv"), index=False)
    print("=== change time by type (s from epoch onset) ===")
    print(ct.groupby(["salience", "long_epoch", "change_pulse"], dropna=False).agg(
        n=("change_s", "size"), median_s=("change_s", "median"), min_s=("change_s", "min"),
        max_s=("change_s", "max"), first_step_rel=("first_step_rel", "median"),
        frac_total=("frac_total_at_change", "median"), dur=("duration", "median")).round(3).to_string())
    print(f"No Change: max |period/base - 1| = {ct.max_rel_nochange.max():.4f}")
    man = g.run_list(); rows = []; jrows = []
    firnames = [col(t, b) for t in TYPES for b in range(NBIN)]
    for sub, runs in man.groupby("subject"):
        dms = []
        for r in runs.itertuples():
            dm, cnt = fir_design(sub, r.bids_entity, int(r.n_vol)); dms.append(dm)
            v = vif(dm, firnames + ["free", "rating", "subrating", "paced_long"])
            fv = [v[n] for n in firnames if n in v]
            rec = dict(subject=sub, entity=r.bids_entity, n_reg=dm.shape[1], rank_deficient=int(rank_def(dm)),
                       **{f"n_{k}": c for k, c in cnt.items()},
                       vif_fir_median=np.median(fv) if fv else np.nan, vif_fir_max=max(fv) if fv else np.nan,
                       vif_free=v.get("free"), vif_rating=v.get("rating"), vif_subrating=v.get("subrating"))
            dmc, gone_r = drop_censored(dm)
            rec["n_fir_censored"] = len(gone_r); rec["rank_deficient_after_drop"] = int(rank_def(dmc))
            rows.append(rec)
        J, _, _, gone = joint_design(sub, runs)
        C = window_specs(J.columns)
        v = vif(J, firnames + ["free", "rating", "subrating"])
        fv = [v[n] for n in firnames if n in v]
        jr = dict(subject=sub, n_runs=len(dms), n_reg=J.shape[1], rank_deficient=int(rank_def(J)),
                  fir_censored=";".join(gone),
                  n_step=int(sum(r_["n_step"] for r_ in rows if r_["subject"] == sub)),
                  n_ramp=int(sum(r_["n_ramp"] for r_ in rows if r_["subject"] == sub)),
                  n_nochange=int(sum(r_["n_nochange"] for r_ in rows if r_["subject"] == sub)),
                  vif_fir_median=np.median(fv) if fv else np.nan, vif_fir_max=max(fv) if fv else np.nan,
                  vif_rating=v.get("rating"), vif_free=v.get("free"))
        for name in ("E7timing", "E7stepNC", "E7rampNC"):
            jr[f"var_{name}"] = hj.cvar(J, C[name]) if name in C else np.nan
            jr[f"nbins_{name}"] = sum(1 for k in C.get(name, {}) if k.startswith("paced_step") or
                                      (name == "E7rampNC" and k.startswith("paced_ramp"))) // (2 if name == "E7timing" else 1)
        jrows.append(jr)
    d = pd.DataFrame(rows); j = pd.DataFrame(jrows)
    d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_runs.csv"), index=False)
    j.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_joint.csv"), index=False)
    print(f"\n=== per-run designs: {len(d)} runs, {d.subject.nunique()} participants ===")
    print(f"rank-deficient: {d.rank_deficient.sum()} (after dropping fully censored FIR bins: "
          f"{d.rank_deficient_after_drop.sum()}); regressors median {d.n_reg.median():.0f}")
    print(f"epochs per run (median): step {d.n_step.median()}, ramp {d.n_ramp.median()}, "
          f"No Change {d.n_nochange.median()}; long {d.n_long.sum()} in {(d.n_long > 0).sum()} runs")
    print("runs lacking a type:", {t: int((d[f'n_{t}'] == 0).sum()) for t in TYPES})
    print("VIF per run: FIR median-of-medians %.2f, FIR max %.2f; free median %.2f; rating median %.2f [%.2f, %.2f]"
          % (d.vif_fir_median.median(), d.vif_fir_max.max(), d.vif_free.median(), d.vif_rating.median(),
             d.vif_rating.min(), d.vif_rating.max()))
    print(f"\n=== joint designs: {len(j)} participants ===")
    print(f"rank-deficient (after dropping censored FIR bins): {j.rank_deficient.sum()}; "
          f"participants with a censored bin: {(j.fir_censored != '').sum()} -> {dict(zip(j.subject[j.fir_censored != ''], j.fir_censored[j.fir_censored != '']))}")
    print("VIF joint: FIR median %.2f, FIR max %.2f; rating median %.2f"
          % (j.vif_fir_median.median(), j.vif_fir_max.max(), j.vif_rating.median()))
    for name in ("E7timing", "E7stepNC", "E7rampNC"):
        ok = j[f"var_{name}"].notna()
        print(f"{name}: estimable in {ok.sum()} / {len(j)}; variance multiplier median "
              f"{j.loc[ok, f'var_{name}'].median():.3f} [{j.loc[ok, f'var_{name}'].min():.3f}, "
              f"{j.loc[ok, f'var_{name}'].max():.3f}]; not estimable: {sorted(j.subject[~ok])}")


def subject(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    man = g.run_list(); runs = man[man.subject == sub]; assert len(runs), sub
    od = os.path.join(g.OUT, "firstlevel", sub); os.makedirs(od, exist_ok=True)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    J, dms, cnts, gone = joint_design(sub, runs)
    specs = {n: s for n, s in window_specs(J.columns).items() if np.isfinite(hj.cvar(J, s))}
    for t in TYPES:
        for b in range(NBIN):
            if col(t, b) in J.columns:
                specs[col(t, b)] = {col(t, b): 1.0}
    log = dict(subject=sub, counts=cnts, n_reg=J.shape[1], fir_censored=gone,
               specs={k: v for k, v in specs.items() if k.startswith("E7")})
    if not any(k.startswith("E7") for k in specs):
        log["contrasts"] = []
        json.dump(log, open(os.path.join(od, f"{sub}_e7.json"), "w"))
        print(sub, "no FIR contrast estimable", log); return
    Ys = []
    for r in runs.itertuples():
        Y = masker.transform(g.bold_path(sub, r.bids_entity)); assert Y.shape[0] == int(r.n_vol)
        Ys.append(mean_scaling(Y, axis=0)[0])
    out, dof = fit(np.vstack(Ys), J, specs)
    def save(arr, name):
        if "variance" in name and not (arr > 0).all():
            raise ValueError(f"{sub} {name}: {(arr <= 0).sum()} non-positive variances")
        g.save_float32(masker.inverse_transform(arr), os.path.join(od, name))
    for name in [k for k in out if k.startswith("E7")]:
        e, v = out[name]
        save(e, f"{sub}_var-base_con-{name}joint_effect.nii.gz")
        save(v, f"{sub}_var-base_con-{name}joint_variance.nii.gz")
    for t in TYPES:
        if any(col(t, b) in out for b in range(NBIN)):
            nv = next(len(out[col(t, b)][0]) for b in range(NBIN) if col(t, b) in out)
            E = np.vstack([out[col(t, b)][0] if col(t, b) in out else np.full(nv, np.nan, np.float32)
                           for b in range(NBIN)])          # a censored bin is saved as NaN
            save(E, f"{sub}_var-base_con-E7fir{t}joint_effect.nii.gz")      # 4-D, 15 bins
    log.update(dof=dof, contrasts=sorted(k for k in out if k.startswith("E7")))
    json.dump(log, open(os.path.join(od, f"{sub}_e7.json"), "w"))
    print(sub, "done", log["contrasts"], "dof", dof)


def group():
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.second_level import non_parametric_inference
    import intero2023_group as G
    from intero2023_explore_awareness import jzs_bf10
    man = g.run_list(); allsubs = sorted(man.subject.unique())
    rows = []; ref = nib.load(g.GROUP_MASK)
    for name, role in (("E7timing", "primary"), ("E7stepNC", "secondary"), ("E7rampNC", "secondary")):
        subs = [s for s in allsubs if os.path.exists(G.fl(s, "base", f"{name}joint", "effect"))]
        files = [G.fl(s, "base", f"{name}joint", "effect") for s in subs]
        imgs = [nib.load(f) for f in files]
        dm = pd.DataFrame({"intercept": np.ones(len(files))})
        json.dump(dict(N=len(subs), subjects=subs),
                  open(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_subjects.json"), "w"))
        for roi in G.SETS["H3"]:
            rmask = nib.load(G.ROI(roi)).get_fdata() > 0
            lp = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="intercept",
                                          mask=G.ROI(roi), n_perm=G.N_PERM, two_sided_test=True,
                                          random_state=G.SEED, n_jobs=G.N_JOBS)
            tmap = G.ols_t(imgs, dm, "intercept", G.ROI(roi))
            lpv = lp.get_fdata()[rmask]; tv = tmap.get_fdata()[rmask]
            k = int(np.nanargmax(np.abs(tv))); ijk = np.argwhere(rmask)[k]
            pmin = float(10 ** -np.nanmax(lpv))
            mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
            y = np.array([mk.transform(f).mean() for f in files])
            row = dict(analysis=f"EXPLORATORY {TAG}", contrast=name, role=role, roi=roi, N=len(files),
                       peak_abs_t_signed=round(float(tv[k]), 3),
                       peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(tmap.affine, ijk)],
                       p_fwe=pmin, p_fwe_bonf3=min(1.0, 3 * pmin), n_vox_p_lt_05=int((lpv > -np.log10(0.05)).sum()),
                       mask_mean=float(y.mean()), mask_mean_dz=float(y.mean() / y.std(ddof=1)),
                       mask_mean_t=float(stats.ttest_1samp(y, 0).statistic))
            for r in (0.5, 0.707, 1.0):
                row[f"BF10_r{r}"] = jzs_bf10(row["mask_mean_t"], len(y), r)
            rows.append(row); print(pd.Series(row).to_dict())
        if role == "primary":
            od = os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}_{name}"); os.makedirs(od, exist_ok=True)
            wb = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept",
                                          mask=g.GROUP_MASK, n_perm=G.N_PERM, two_sided_test=True,
                                          threshold=0.002, random_state=G.SEED, n_jobs=G.N_JOBS)
            for kk in ("t", "logp_max_size"):
                g.save_float32(wb[kk], os.path.join(od, f"wholebrain_{kk}.nii.gz"))
            cl = G.cluster_table(wb["t"], wb["logp_max_size"], len(files) - 1)
            cl.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_wholebrain.csv"), index=False)
            print("whole brain:", cl.head(10).to_string() if len(cl) else "no cluster")
    pd.DataFrame(rows).to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_fir_h3masks.csv"), index=False)
    # descriptive time courses: per-bin group mean +- SE of the mask-mean FIR betas
    subs = [s for s in allsubs if all(os.path.exists(G.fl(s, "base", f"E7fir{t}joint", "effect")) for t in TYPES)]
    tc = []
    for roi in G.SETS["H3"]:
        mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
        B = {t: np.vstack([np.nanmean(mk.transform(G.fl(s, "base", f"E7fir{t}joint", "effect")), 1) for s in subs])
             for t in TYPES}                                    # participants x bins
        series = dict(B, step_minus_nochange=B["step"] - B["nochange"],
                      ramp_minus_nochange=B["ramp"] - B["nochange"], step_minus_ramp=B["step"] - B["ramp"])
        for name, M in series.items():
            for b in range(NBIN):
                x = M[:, b]; x = x[np.isfinite(x)]
                tc.append(dict(roi=roi, series=name, bin=b, t_start_s=b * BINW, t_mid_s=b * BINW + 1, N=len(x),
                               mean=float(x.mean()), se=float(x.std(ddof=1) / np.sqrt(len(x)))))
    tc = pd.DataFrame(tc)
    tc.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_timecourse_h3masks.csv"), index=False)
    plot(tc)


def plot(tc):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rois = list(dict.fromkeys(tc.roi))
    fig, ax = plt.subplots(2, len(rois), figsize=(4.2 * len(rois), 7), sharex=True)
    for j, roi in enumerate(rois):
        for row, sers in enumerate((TYPES, ("step_minus_nochange", "ramp_minus_nochange", "step_minus_ramp"))):
            a = ax[row, j]
            for s, c in zip(sers, ("#c0392b", "#2471a3", "#7f8c8d")):
                d = tc[(tc.roi == roi) & (tc.series == s)]
                a.plot(d.t_mid_s, d["mean"], color=c, label=s); a.fill_between(d.t_mid_s, d["mean"] - d.se, d["mean"] + d.se, color=c, alpha=.2)
            for x, ls in ((4.5, ":"), (9.5, "--")):
                a.axvline(x, color="k", ls=ls, lw=.8)
            for lo, hi in ((6, 12), (12, 18)):
                a.axvspan(lo, hi, color="0.9", zorder=0) if row == 1 and lo == 12 else None
            a.axhline(0, color="k", lw=.5); a.set_title(roi if row == 0 else ""); a.legend(fontsize=7)
            if row == 1:
                a.set_xlabel("s from epoch onset")
        ax[0, 0].set_ylabel("FIR beta (PSC), no implicit baseline"); ax[1, 0].set_ylabel("difference (PSC)")
    fig.suptitle("EXPLORATORY E7: FIR time courses (mean ± SE); dotted ramp start 4.5 s, dashed step 9.5 s")
    fig.tight_layout()
    od = os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}_E7timing"); os.makedirs(od, exist_ok=True)
    fig.savefig(os.path.join(od, "E7_timecourse_h3masks.png"), dpi=130)


if __name__ == "__main__":
    a = sys.argv[1]
    {"qc": qc, "subject": lambda: subject(sys.argv[2]), "group": group}[a]()
