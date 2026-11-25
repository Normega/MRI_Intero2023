#!/usr/bin/env python
"""
intero2023_designqc.py -- ANALYSIS_HANDOFF.md step 3.3, design-matrix QC.

Builds every first-level design matrix exactly as registered (prereg S7) and
reports collinearity. Touches NO BOLD data: nothing here is parameter
estimation, so it is legal before the OSF freeze.

Registered design, per run:
  SPM canonical HRF + temporal derivative, drift_model=None (the cosine columns
  in the confound sheet ARE the high-pass), TR = 2 s, 7 condition boxcars,
  confidence parametric modulator on paced epochs (mean-centred over the trials
  actually modelled), the base nuisance columns of the confound sheet (physio
  columns dropped: they enter sensitivity models only). No global signal, ever.

Fixed 2026-11-25: the SE multiplier scaled the contrast by the column SDs
instead of dividing by them, so SEs were wrong whenever the contrasted columns
differed in scale; and physio columns were in the design.

Reports, per run:
  * VIF of each task regressor AMONG THE TASK REGRESSORS (handoff step 3.3).
    The intercept is deliberately excluded from the reference set: this design
    tiles the run (paced + free + rating ~ constant), so a VIF taken against
    the constant diverges for every run by construction and diagnoses nothing.
  * SE multiplier sqrt(c' (X'X)^-1 c) for each REGISTERED contrast, computed in
    the full design. This is the number that decides anything: it says how
    precisely the contrast we actually test can be estimated, and it catches
    rank deficiency that per-regressor VIF hides.
  * r(confidence PM, aware-minus-unaware contrast direction), the specific
    confusion H2 vs H3 would suffer from.

Usage:  python Scripts/intero2023_designqc.py [registered|candidates|both]
"""
import sys, warnings
import numpy as np, pandas as pd
from nilearn.glm.first_level import make_first_level_design_matrix

warnings.simplefilter("ignore")
TR = 2.0
FL = "FirstLevel"
CONDITIONS = ["paced_aware", "paced_unaware", "paced_wrongdir", "paced_nochange",
              "free", "rating", "subrating"]
PACED_MOD = ["paced_aware", "paced_unaware", "paced_wrongdir", "paced_nochange"]
# Physiological columns belong to the sensitivity models only (prereg S6-S7).
# Design QC describes the BASE model, which is what every hypothesis is fitted
# in, so they are dropped here (fixed 2026-11-25; earlier versions kept them).
PHYSIO = ["rvt_rrf_z", "rvt_z", "resp_rate_c",
          "petco2_hrf_lag00", "petco2_hrf_lag10", "petco2_hrf_lag20", "petco2_hrf_lag30"]


def build(subject, entity, n_vol):
    """Return (design_matrix, task_column_names). Registered spec, no deviations."""
    ev = pd.read_csv(f"{FL}/{subject}_{entity}_desc-model_events.tsv", sep="\t")
    cf = pd.read_csv(f"{FL}/{subject}_{entity}_desc-model_confounds.tsv", sep="\t")
    assert len(cf) == n_vol, f"{subject} {entity}: sheet {len(cf)} rows vs n_vol {n_vol}"
    cf = cf.drop(columns=[c for c in PHYSIO if c in cf.columns])

    ev = ev[ev.condition.isin(CONDITIONS)].copy()
    # Trailing events can over-run the acquisition (scanner stopped mid-rating).
    # Truncate to the scan and log it; never drop silently. See VERIFICATION_REPORT.
    scan_end = n_vol * TR
    over = (ev.onset + ev.duration) > scan_end
    n_trunc = int(over.sum())
    ev.loc[over, "duration"] = (scan_end - ev.loc[over, "onset"]).clip(lower=0)
    ev = ev[ev.duration > 0]

    events = ev[["onset", "duration", "condition"]].rename(columns={"condition": "trial_type"})
    events["modulation"] = 1.0

    # Confidence PM on paced epochs, mean-centred over the trials ACTUALLY modelled.
    pm = ev[ev.condition.isin(PACED_MOD) & ev.confidence.notna()].copy()
    if len(pm) >= 2 and pm.confidence.std() > 0:
        pm = pm[["onset", "duration", "confidence"]].rename(columns={"confidence": "modulation"})
        pm["modulation"] = pm.modulation - pm.modulation.mean()
        pm["trial_type"] = "conf_pm"
        events = pd.concat([events, pm], ignore_index=True)

    frame_times = TR * (np.arange(n_vol) + 0.5)
    dm = make_first_level_design_matrix(
        frame_times, events=events, hrf_model="spm + derivative",
        drift_model=None, add_regs=cf.values, add_reg_names=list(cf.columns))
    task = [c for c in dm.columns
            if c in CONDITIONS + ["conf_pm"] or c.endswith("_derivative")]
    return dm, task, n_trunc


def vifs(dm, cols):
    """VIF of each named column against the OTHER TASK columns (no intercept)."""
    ref = [c for c in cols]
    X = dm[ref].values.astype(float)
    names = ref
    out = {}
    for c in cols:
        j = names.index(c)
        y = X[:, j]
        Z = np.delete(X, j, axis=1)
        if Z.shape[1] == 0:
            out[c] = 1.0
            continue
        beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
        resid = y - Z @ beta
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else 0.0
        out[c] = np.inf if r2 >= 1 - 1e-12 else 1.0 / (1.0 - r2)
    return out


CONTRASTS = {
    # prereg S7 / handoff S4. paced side of H1 is aware+unaware+nochange, i.e.
    # wrongdir and noresponse are regressors of no interest and stay out of it.
    "paced_gt_free":     {"paced_aware": 1/3, "paced_unaware": 1/3,
                          "paced_nochange": 1/3, "free": -1.0},
    "confidence_mod":    {"conf_pm": 1.0},
    "aware_gt_unaware":  {"paced_aware": 1.0, "paced_unaware": -1.0},
}

# Amendment candidate (see VERIFICATION_REPORT / prereg S4 correction): the
# registered paced_gt_free is a FIXED 3-way average, so it is undefined in any
# run lacking one of the three paced conditions -- which is 18 of 91 runs,
# almost all of them missing paced_unaware because the participant never
# answered "No Change" on a change trial. The adaptive form averages over the
# paced conditions actually present in the run. Reported alongside, never
# substituted silently.
PACED_SIDE = ["paced_aware", "paced_unaware", "paced_nochange"]


def contrast_se(dm):
    """sqrt(c' (X'X)^-1 c) per registered contrast; nan if not estimable here."""
    names = list(dm.columns)
    X = dm.values.astype(float)
    # Column-standardise before any rank/estimability algebra. The raw design
    # mixes boxcars (~1) with csf/white_matter (~10^3); without this, pinv's
    # singular-value cutoff reports spurious rank deficiency and spurious
    # "not estimable" verdicts. Scaling changes no contrast's meaning here
    # because the SE multiplier is rescaled back below.
    sd = X.std(axis=0)
    sd[sd < 1e-12] = 1.0
    Xs = X / sd
    u, sv, vt = np.linalg.svd(Xs, full_matrices=False)
    tol = max(Xs.shape) * np.finfo(float).eps * sv.max()
    rank = int((sv > tol).sum())
    full_rank = rank == Xs.shape[1]
    XtX_inv = np.linalg.pinv(Xs.T @ Xs, rcond=1e-12)
    Vr = vt[:rank]                       # orthonormal basis of the row space
    out = {}
    for cname, spec in CONTRASTS.items():
        if not set(spec) <= set(names):
            out[cname] = np.nan          # a needed condition is absent this run
            continue
        c = np.zeros(len(names))
        for k, w in spec.items():
            c[names.index(k)] = w
        cs = c / sd                      # contrast in standardised coordinates:
        # beta_s = sd * beta, so c'beta = (c/sd)'beta_s (fixed 2026-11-25; was c*sd)
        if full_rank:
            estimable = True
        else:
            proj = Vr.T @ (Vr @ cs)
            estimable = np.linalg.norm(cs - proj) <= 1e-8 * max(1.0, np.linalg.norm(cs))
        out[cname] = float(np.sqrt(cs @ XtX_inv @ cs)) if estimable else np.nan
    # adaptive paced_gt_free
    present = [k for k in PACED_SIDE if k in names]
    if present:
        c = np.zeros(len(names))
        for k in present:
            c[names.index(k)] = 1.0 / len(present)
        c[names.index("free")] = -1.0
        cs = c / sd
        ok = full_rank or (np.linalg.norm(cs - Vr.T @ (Vr @ cs))
                           <= 1e-8 * max(1.0, np.linalg.norm(cs)))
        out["paced_gt_free_adaptive"] = float(np.sqrt(cs @ XtX_inv @ cs)) if ok else np.nan
        out["n_paced_conditions"] = len(present)
    else:
        out["paced_gt_free_adaptive"] = np.nan
        out["n_paced_conditions"] = 0
    out["rank_deficient"] = int(not full_rank)
    return out


def main(which):
    mans = []
    if which in ("registered", "both"):
        m = pd.read_csv(f"{FL}/firstlevel_manifest.csv"); m["set"] = "registered"; mans.append(m)
    if which in ("candidates", "both"):
        m = pd.read_csv(f"{FL}/firstlevel_manifest_candidates.csv"); m["set"] = "candidate"; mans.append(m)
    man = pd.concat(mans, ignore_index=True)
    keep = man[man.motion_gate_fail == False]

    rows = []
    for _, r in keep.iterrows():
        dm, task, n_trunc = build(r.subject, r.bids_entity, int(r.n_vol))
        main_cols = [c for c in task if not c.endswith("_derivative")]
        v = vifs(dm, main_cols)
        se = contrast_se(dm)
        rec = dict(subject=r.subject, entity=r.bids_entity, set=r["set"],
                   n_reg=dm.shape[1], n_vol=int(r.n_vol),
                   df_resid=int(r.n_vol) - dm.shape[1], n_truncated_events=n_trunc)
        rec.update({f"vif_{k}": val for k, val in v.items()})
        rec.update({f"se_{k}": val for k, val in se.items()})
        if {"paced_aware", "paced_unaware", "conf_pm"} <= set(dm.columns):
            d = dm["paced_aware"].values - dm["paced_unaware"].values
            rec["r_confpm_awareVunaware"] = float(np.corrcoef(d, dm["conf_pm"].values)[0, 1])
        rows.append(rec)

    d = pd.DataFrame(rows)
    d.to_csv("../analysis/design_qc.csv", index=False)
    vcols = [c for c in d.columns if c.startswith("vif_")]
    scols = [c for c in d.columns if c.startswith("se_")
             and c not in ("se_rank_deficient", "se_n_paced_conditions")]

    print(f"runs: {len(d)}   sets: {d['set'].value_counts().to_dict()}")
    print(f"\nregressors per design: min {d.n_reg.min()} median {int(d.n_reg.median())} max {d.n_reg.max()}")
    print(f"residual df:           min {d.df_resid.min()} median {int(d.df_resid.median())} max {d.df_resid.max()}")
    print(f"runs with truncated trailing events: {int((d.n_truncated_events>0).sum())}")
    if (d.n_truncated_events > 0).any():
        for _, r in d[d.n_truncated_events > 0].iterrows():
            print(f"   {r.subject} {r.entity}: {r.n_truncated_events} event(s) clipped to scan end")

    print(f"rank-deficient designs: {int(d.se_rank_deficient.sum())} / {len(d)}")
    print("\n=== VIF of task regressors among task regressors (no intercept) ===")
    print(f"{'regressor':22s} {'n runs':>7s} {'median':>8s} {'p90':>8s} {'max':>9s} {'>5':>4s} {'>10':>4s}")
    for c in vcols:
        s = d[c].dropna()
        if not len(s): continue
        print(f"{c[4:]:22s} {len(s):7d} {s.median():8.2f} {s.quantile(.9):8.2f} "
              f"{s.max():9.2f} {int((s>5).sum()):4d} {int((s>10).sum()):4d}")

    print("\n=== SE multiplier sqrt(c'(X'X)^-1 c) for the registered contrasts ===")
    print(f"{'contrast':22s} {'n est':>6s} {'n NOT':>6s} {'median':>9s} {'p90':>9s} {'max':>9s}")
    for c in scols:
        s_ = d[c]
        ok = s_.dropna()
        print(f"{c[3:]:22s} {len(ok):6d} {int(s_.isna().sum()):6d} "
              f"{ok.median():9.4f} {ok.quantile(.9):9.4f} {ok.max():9.4f}")
        bad = d[s_.isna()]
        if len(bad) and len(bad) <= 25:
            print("     not estimable: " + ", ".join(f"{r.subject}/{r.entity[-5:]}" for _, r in bad.iterrows()))

    print("\n=== r(confidence PM, aware-minus-unaware direction) ===")
    s = d["r_confpm_awareVunaware"].dropna()
    print(f"  n={len(s)}  median |r|={s.abs().median():.3f}  p90 |r|={s.abs().quantile(.9):.3f}  max |r|={s.abs().max():.3f}")
    print(f"  runs with |r| > 0.5: {int((s.abs()>0.5).sum())}")
    worst = d.reindex(d["r_confpm_awareVunaware"].abs().sort_values(ascending=False).index).head(5)
    print(worst[["subject", "entity", "r_confpm_awareVunaware"]].to_string(index=False))
    print("\nwrote analysis/design_qc.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "both")
