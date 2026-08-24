#!/usr/bin/env python
"""
intero2023_twomodel_qc.py -- design QC for the two-model amendment.

The registered spec puts the awareness conditions and the confidence parametric
modulator in ONE first-level model. intero2023_designqc.py showed what that
costs: conf_pm is aliased against the awareness boxcars (median |r| with the
aware-minus-unaware direction = .437, |r| > .5 in 26/70 runs, max .942) and it
makes 9/91 designs rank-deficient on its own -- drop conf_pm and all 9 recover.

The amendment splits it in two:

  Model A (authoritative, H1/H3/H6): the registered condition boxcars, NO
    parametric modulator.
  Model B (supportive, H2): the paced epochs collapsed to a single `trial`
    regressor carrying the confidence parametric modulator, plus free / rating
    / subrating unchanged.

This script builds both, per run, and reports rank, contrast estimability and
participant-level Ns. It fits nothing and reads no BOLD: legal pre-freeze.

Model B has one genuinely ambiguous choice -- whether `trial` spans the change
trials only or all paced epochs including paced_nochange -- so BOTH variants are
computed and reported side by side.

Usage:  python Scripts/intero2023_twomodel_qc.py [registered|candidates|both]
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
PACED_SIDE = ["paced_aware", "paced_unaware", "paced_nochange"]   # H1 paced side
CHANGE_ONLY = ["paced_aware", "paced_unaware", "paced_wrongdir"]


def read_run(subject, entity, n_vol):
    """Events (clipped to the scan) and confounds for one run."""
    ev = pd.read_csv(f"{FL}/{subject}_{entity}_desc-model_events.tsv", sep="\t")
    cf = pd.read_csv(f"{FL}/{subject}_{entity}_desc-model_confounds.tsv", sep="\t")
    assert len(cf) == n_vol, f"{subject} {entity}: sheet {len(cf)} rows vs n_vol {n_vol}"
    ev = ev[ev.condition.isin(CONDITIONS)].copy()
    scan_end = n_vol * TR
    over = (ev.onset + ev.duration) > scan_end
    n_trunc = int(over.sum())
    ev.loc[over, "duration"] = (scan_end - ev.loc[over, "onset"]).clip(lower=0)
    ev = ev[ev.duration > 0]
    return ev, cf, n_trunc


def design(ev_frame, cf, n_vol):
    frame_times = TR * (np.arange(n_vol) + 0.5)
    return make_first_level_design_matrix(
        frame_times, events=ev_frame, hrf_model="spm + derivative",
        drift_model=None, add_regs=cf.values, add_reg_names=list(cf.columns))


def model_a(ev, cf, n_vol):
    e = ev[["onset", "duration", "condition"]].rename(columns={"condition": "trial_type"})
    e = e.assign(modulation=1.0)
    return design(e, cf, n_vol)


def model_b(ev, cf, n_vol, trial_set):
    """Paced epochs in `trial_set` collapsed to one regressor + confidence PM."""
    other = ev[~ev.condition.isin(PACED_MOD)][["onset", "duration", "condition"]] \
        .rename(columns={"condition": "trial_type"}).assign(modulation=1.0)
    tr = ev[ev.condition.isin(trial_set)].copy()
    if len(tr) == 0:
        return None
    trial = tr[["onset", "duration"]].assign(trial_type="trial", modulation=1.0)
    e = pd.concat([other, trial], ignore_index=True)
    pm = tr[tr.confidence.notna()].copy()
    if len(pm) >= 2 and pm.confidence.std() > 0:
        pm = pm[["onset", "duration", "confidence"]].rename(columns={"confidence": "modulation"})
        pm["modulation"] = pm.modulation - pm.modulation.mean()
        pm["trial_type"] = "conf_pm"
        e = pd.concat([e, pm], ignore_index=True)
    return design(e, cf, n_vol)


def se_of(dm, spec):
    """sqrt(c'(X'X)^-1 c); nan if a needed column is absent or c is not estimable.

    Column-standardised before the rank algebra for the reason documented in
    intero2023_designqc.py: the raw design mixes boxcars (~1) with csf (~10^3).
    """
    names = list(dm.columns)
    if not set(spec) <= set(names):
        return np.nan, None
    X = dm.values.astype(float)
    sd = X.std(axis=0); sd[sd < 1e-12] = 1.0
    Xs = X / sd
    u, sv, vt = np.linalg.svd(Xs, full_matrices=False)
    tol = max(Xs.shape) * np.finfo(float).eps * sv.max()
    rank = int((sv > tol).sum())
    full_rank = rank == Xs.shape[1]
    c = np.zeros(len(names))
    for k, w in spec.items():
        c[names.index(k)] = w
    cs = c * sd
    if not full_rank:
        Vr = vt[:rank]
        if np.linalg.norm(cs - Vr.T @ (Vr @ cs)) > 1e-8 * max(1.0, np.linalg.norm(cs)):
            return np.nan, full_rank
    XtX_inv = np.linalg.pinv(Xs.T @ Xs, rcond=1e-12)
    return float(np.sqrt(cs @ XtX_inv @ cs)), full_rank


def adaptive_h1(dm):
    """paced > free averaged over the paced conditions PRESENT in this run."""
    present = [k for k in PACED_SIDE if k in dm.columns]
    if not present or "free" not in dm.columns:
        return np.nan, 0
    spec = {k: 1.0 / len(present) for k in present}
    spec["free"] = -1.0
    se, _ = se_of(dm, spec)
    return se, len(present)


# Registered estimability rule (prereg S3): a contrast counts for a participant
# only if EVERY cell it involves has >= 3 trials, POOLED across that
# participant's gate-passing runs. This is stricter than matrix estimability --
# a single unaware trial makes the design matrix full rank but is not a cell --
# and it is what produced the registered Ns (C2 = 21 post-gate). Both are
# reported below; the registered rule is the one that governs.
MIN_TRIALS = 3


def main(which):
    mans = []
    if which in ("registered", "both"):
        m = pd.read_csv(f"{FL}/firstlevel_manifest.csv"); m["set"] = "registered"; mans.append(m)
    if which in ("candidates", "both"):
        m = pd.read_csv(f"{FL}/firstlevel_manifest_candidates.csv")
        # sub-056 is NOT admitted: run-2 fails the registered motion gate and
        # run-1 over-runs the acquisition by 28 s (VERIFICATION_REPORT S6).
        m = m[m.subject != "sub-056"]
        m["set"] = "candidate"; mans.append(m)
    man = pd.concat(mans, ignore_index=True)
    keep = man[man.motion_gate_fail == False]

    rows = []
    for _, r in keep.iterrows():
        ev, cf, n_trunc = read_run(r.subject, r.bids_entity, int(r.n_vol))
        rec = dict(subject=r.subject, entity=r.bids_entity, set=r["set"],
                   n_vol=int(r.n_vol), n_truncated_events=n_trunc)

        dmA = model_a(ev, cf, int(r.n_vol))
        _, fullA = se_of(dmA, {"free": 1.0})
        rec["A_rank_deficient"] = int(fullA is False)
        rec["A_n_reg"] = dmA.shape[1]
        rec["se_A_h1_adaptive"], rec["A_n_paced"] = adaptive_h1(dmA)
        rec["se_A_h1_registered"], _ = se_of(
            dmA, {"paced_aware": 1/3, "paced_unaware": 1/3, "paced_nochange": 1/3, "free": -1.0})
        rec["se_A_h3_aware_gt_unaware"], _ = se_of(
            dmA, {"paced_aware": 1.0, "paced_unaware": -1.0})

        for tag, tset in (("chg", CHANGE_ONLY), ("all", PACED_MOD)):
            dmB = model_b(ev, cf, int(r.n_vol), tset)
            if dmB is None:
                rec[f"B{tag}_rank_deficient"] = np.nan
                rec[f"se_B{tag}_h2_conf"] = np.nan
                continue
            _, fullB = se_of(dmB, {"free": 1.0})
            rec[f"B{tag}_rank_deficient"] = int(fullB is False)
            rec[f"B{tag}_n_reg"] = dmB.shape[1]
            rec[f"se_B{tag}_h2_conf"], _ = se_of(dmB, {"conf_pm": 1.0})
        rows.append(rec)

    d = pd.DataFrame(rows)
    d.to_csv("../analysis/design_qc_twomodel.csv", index=False)

    subs = sorted(d.subject.unique())
    maia = set(pd.read_csv("Data/maia_scores.csv").subject)

    # Pooled per-participant trial counts, for the registered >=3-per-cell rule.
    cnt = keep.groupby("subject")[["n_aware", "n_unaware", "n_wrongdir",
                                   "n_nochange", "n_free"]].sum()
    cnt["n_paced_h1"] = cnt.n_aware + cnt.n_unaware + cnt.n_nochange
    cnt["n_trial_chg"] = cnt.n_aware + cnt.n_unaware + cnt.n_wrongdir
    cnt["n_trial_all"] = cnt.n_trial_chg + cnt.n_nochange
    # confidence must also VARY, or the PM is a constant column
    confvar = {}
    for s_, g in keep.groupby("subject"):
        vals = []
        for _, r in g.iterrows():
            ev, _, _ = read_run(r.subject, r.bids_entity, int(r.n_vol))
            vals += list(ev[ev.condition.isin(PACED_MOD)].confidence.dropna())
        confvar[s_] = (len(vals) >= MIN_TRIALS) and (float(np.std(vals)) > 0)

    CELLS = {  # contrast -> pooled cells that must each reach MIN_TRIALS
        "se_A_h1_adaptive":            ["n_paced_h1", "n_free"],
        # C1's registered cells are "paced" and "free", not the three paced
        # sub-conditions; the fixed 1/3 form's extra requirement that all three
        # be PRESENT is enforced by matrix estimability, not by a trial count.
        "se_A_h1_registered":          ["n_paced_h1", "n_free"],
        "se_A_h3_aware_gt_unaware":    ["n_aware", "n_unaware"],
        "se_Bchg_h2_conf":             ["n_trial_chg"],
        "se_Ball_h2_conf":             ["n_trial_all"],
    }

    def n_with(col, registered=False):
        ok = sorted(d[d[col].notna()].subject.unique())
        if not registered:
            return ok
        cells = CELLS[col]
        ok = [s_ for s_ in ok if s_ in cnt.index
              and all(cnt.loc[s_, c] >= MIN_TRIALS for c in cells)]
        if col.endswith("h2_conf"):
            ok = [s_ for s_ in ok if confvar.get(s_, False)]
        return ok

    print(f"runs: {len(d)}   participants: {len(subs)}   "
          f"sets: {d['set'].value_counts().to_dict()}")
    print(f"\nregressors per design: A median {int(d.A_n_reg.median())}, "
          f"B(chg) median {int(d.Bchg_n_reg.median())}, B(all) median {int(d.Ball_n_reg.median())}")
    print("\n=== rank deficiency (the reason for the split) ===")
    print(f"  Model A                 : {int(d.A_rank_deficient.sum())} / {len(d)} runs")
    print(f"  Model B (change trials) : {int(d.Bchg_rank_deficient.sum(skipna=True))} / {int(d.Bchg_rank_deficient.notna().sum())} runs")
    print(f"  Model B (all paced)     : {int(d.Ball_rank_deficient.sum(skipna=True))} / {int(d.Ball_rank_deficient.notna().sum())} runs")
    print("  (registered single model, for comparison: 9 / 91 -- see analysis/design_qc.csv)")

    print("\n=== participant-level N per hypothesis ===")
    rowsN = [
        ("H1 paced>free, adaptive", "A", "se_A_h1_adaptive"),
        ("H1 paced>free, registered fixed 1/3", "A", "se_A_h1_registered"),
        ("H3 aware>unaware", "A", "se_A_h3_aware_gt_unaware"),
        ("H2 confidence PM, trial = change trials", "B(chg)", "se_Bchg_h2_conf"),
        ("H2 confidence PM, trial = all paced", "B(all)", "se_Ball_h2_conf"),
    ]
    print(f"{'contrast':42s} {'model':7s} {'runs':>5s} {'N matrix':>9s} {'N registered':>13s}")
    got = {}
    for label, mod, col in rowsN:
        loose = n_with(col)
        strict = n_with(col, registered=True)
        got[col] = strict
        print(f"{label:42s} {mod:7s} {int(d[col].notna().sum()):5d} "
              f"{len(loose):9d} {len(strict):13d}")

    h1 = got["se_A_h1_adaptive"]
    print(f"\nH6 (MAIA covariate on H1, Model A): N = {len([s for s in h1 if s in maia])}")

    h3 = got["se_A_h3_aware_gt_unaware"]
    print(f"\n=== is the H3 subsample selected? ===")
    man2 = cnt
    inn = man2.loc[[s for s in h3 if s in man2.index]]
    out = man2.loc[[s for s in subs if s not in h3 and s in man2.index]]
    print(f"  H3 in  (N={len(inn)}): median unaware {inn.n_unaware.median():.1f}, "
          f"median aware {inn.n_aware.median():.1f}, "
          f"median unaware share {(inn.n_unaware/(inn.n_aware+inn.n_unaware)).median():.2f}")
    print(f"  H3 out (N={len(out)}): median unaware {out.n_unaware.median():.1f}, "
          f"median aware {out.n_aware.median():.1f}, "
          f"median unaware share {(out.n_unaware/(out.n_aware+out.n_unaware)).median():.2f}")
    print(f"  H3 in : {sorted(inn.index)}")

    print("\n=== precision cost of the split, runs where both are estimable ===")
    both = d[d.se_A_h1_adaptive.notna() & d.se_Ball_h2_conf.notna()]
    print(f"  n={len(both)}  median SE(H1 | Model A) = {both.se_A_h1_adaptive.median():.4f}")
    print(f"           median SE(H2 | Model B all) = {both.se_Ball_h2_conf.median():.4f}")
    print("\nwrote analysis/design_qc_twomodel.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "both")
