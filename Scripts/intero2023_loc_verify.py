#!/usr/bin/env python
"""
intero2023_loc_verify.py -- verification of the localizer first-level inputs and designs.

  python Scripts/intero2023_loc_verify.py designs      # design only, no BOLD values read
  python Scripts/intero2023_loc_verify.py e2e sub-060  # after the first level of sub-060:
                                                       # fit() vs nilearn FirstLevelModel on run-1

`designs` checks every one of the 102 registered runs:
  files exist (preproc BOLD, brain mask, fMRIPrep confounds, events, nuisance sheet);
  n_vol (NIfTI header) == fMRIPrep confound rows == sheet rows;
  sheet motion / FD columns == fMRIPrep's (6 dp) and spike count == FD > 0.5 count;
  events start inside the scan; how many seconds / events are clipped at the scan end;
  every per-run design (base, hp128, hp320, +rvt, +co2 where available) is full rank,
  and so is every participant's joint design, with the contrasts that are estimable.
Writes Analysis/localizer_glm_verification.csv (run rows) and
Analysis/localizer_glm_verification_joint.csv (participant rows).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
g = L.g
warnings.simplefilter("ignore")
MOTION = [f"{a}_{b}" for a in ("trans", "rot") for b in ("x", "y", "z")]


def rank_info(dm):
    X = dm.values.astype(float)
    sd = X.std(axis=0); sd[sd < 1e-12] = 1.0
    Xs = X / sd
    sv = np.linalg.svd(Xs, compute_uv=False)
    rank = int((sv > max(Xs.shape) * np.finfo(float).eps * sv.max()).sum())
    return rank, X.shape[1], float(sv.max() / sv.min()) if sv.min() > 0 else np.inf


def max_task_vif(dm):
    """Largest VIF among the task regressors, against every other column of the design."""
    X = dm.values.astype(float); out = []
    for k in [c for c in L.TASK if c in dm.columns]:
        j = list(dm.columns).index(k); y = X[:, j]; Z = np.delete(X, j, axis=1)
        b, *_ = np.linalg.lstsq(Z, y, rcond=None)
        out.append(1.0 / max(1e-12, ((y - Z @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    return max(out)


def variants(r):
    v = {"base": ((), 200.0), "hp128": ((), 128.0), "hp320": ((), 320.0)}
    if r.has_rvt:
        v["rvt"] = (tuple(L.PHYSIO_RVT), 200.0)
    if r.has_petco2:
        v["co2"] = (tuple(L.PHYSIO_CO2), 200.0)
    return v


def designs():
    u = L.run_list()
    rows, jrows = [], []
    for r in u.itertuples():
        sub, run = r.subject, int(r.run); ent = L.entity(run)
        files = dict(bold=g.bold_path(sub, ent), mask=g.mask_path(sub, ent),
                     confounds=g.fmriprep_confounds_path(sub, ent), events=L.events_path(sub, run),
                     sheet=L.sheet_path(sub, run))
        rec = dict(subject=sub, run=run, protocol_version=r.protocol_version,
                   files_exist=all(os.path.exists(f) for f in files.values()))
        n_hdr = int(nib.load(files["bold"]).shape[3])                      # header only
        fp = pd.read_csv(files["confounds"], sep="\t", na_values=["n/a"])
        sh = pd.read_csv(files["sheet"], sep="\t")
        rec.update(n_vol_header=n_hdr, n_rows_fmriprep=len(fp), n_rows_sheet=len(sh),
                   volumes_match=n_hdr == len(fp) == len(sh) == int(r.n_vol))
        cols = MOTION + [m + "_derivative1" for m in MOTION] + ["framewise_displacement", "csf", "white_matter"]
        d = max(float(np.abs(sh[c].values - fp[c].fillna(0).round(6).values).max()) for c in cols)
        fd = fp.framewise_displacement.fillna(0).values
        n_spk = sum(c.startswith("spike_") for c in sh.columns)
        rec.update(sheet_max_abs_diff=d, n_spikes=n_spk, sheet_matches_fmriprep=(d <= 1e-6) and n_spk == int((fd > 0.5).sum()))
        raw = pd.read_csv(files["events"], sep="\t")
        end = (raw.onset + raw.duration).max(); scan = int(r.n_vol) * L.TR
        ev, n_clip = L.read_events(sub, run, int(r.n_vol))
        rec.update(scan_dur_s=scan, events_end_s=round(float(end), 2), s_past_scan_end=round(max(0.0, end - scan), 2),
                   n_events_clipped=n_clip, events_start_ok=bool((raw.onset >= 0).all() and (raw.onset < scan).all()),
                   events_fit_after_clip=bool(((ev.onset + ev.duration) <= scan + 1e-9).all()),
                   n_intero=int((ev.trial_type == "intero").sum()), n_intero_nc=int((ev.trial_type == "intero_nc").sum()),
                   n_extero=int((ev.trial_type == "extero").sum()), n_match=int((ev.trial_type == "match").sum()),
                   n_self=int((ev.trial_type == "self").sum()), n_instruction=int((ev.trial_type == "instruction").sum()))
        full = []
        for v, (phys, hp) in variants(r).items():
            dm = L.design(ev, L.read_confounds(sub, run, int(r.n_vol), phys), int(r.n_vol), hp)
            rk, p, cond = rank_info(dm)
            rec[f"rank_{v}"] = rk; rec[f"n_reg_{v}"] = p; rec[f"full_rank_{v}"] = rk == p
            full.append(rk == p)
            if v == "base":
                rec["n_drift_200"] = sum(c.startswith("drift") for c in dm.columns)
                rec["vif_max_task_base"] = round(max_task_vif(dm), 2)
                rec["estimable_base"] = ";".join(c for c, s in L.CONTRASTS.items() if L.estimable(dm, s))
        rec["design_full_rank_all_variants"] = all(full)
        rows.append(rec)
    R = pd.DataFrame(rows)
    R["ok"] = R.files_exist & R.volumes_match & R.sheet_matches_fmriprep & R.events_start_ok \
        & R.events_fit_after_clip & R.design_full_rank_all_variants
    R.to_csv(os.path.join(g.REPO, "Analysis", "localizer_glm_verification.csv"), index=False)
    # joint designs
    for sub, runs in u.groupby("subject"):
        rec = dict(subject=sub, n_runs=len(runs))
        V = {"base": ((), 200.0), "hp128": ((), 128.0), "hp320": ((), 320.0)}
        if runs.has_rvt.all():
            V["rvt"] = (tuple(L.PHYSIO_RVT), 200.0)
        if runs.has_petco2.all():
            V["co2"] = (tuple(L.PHYSIO_CO2), 200.0)
        for v, (phys, hp) in V.items():
            J = L.build_joint_design([dm for _, dm in L.run_designs(sub, runs, phys, hp)])
            rk, p, _ = rank_info(J)
            rec[f"full_rank_{v}"] = rk == p; rec[f"n_reg_{v}"] = p
            if v == "base":
                rec["n_vol"] = J.shape[0]; rec["dof_base"] = J.shape[0] - rk
                rec["estimable_base"] = ";".join(c for c, s in L.CONTRASTS.items() if L.estimable(J, s))
        rec["variants"] = ";".join(V)
        jrows.append(rec)
    Jt = pd.DataFrame(jrows)
    Jt.to_csv(os.path.join(g.REPO, "Analysis", "localizer_glm_verification_joint.csv"), index=False)
    s, _, _ = L.participant_sets()
    print(f"runs {len(R)}, participants {R.subject.nunique()}; all checks pass: {int(R.ok.sum())}/{len(R)}")
    for c in ["files_exist", "volumes_match", "sheet_matches_fmriprep", "events_start_ok", "events_fit_after_clip",
              "design_full_rank_all_variants"]:
        print(f"  {c}: {int(R[c].sum())}/{len(R)}")
    print(f"  runs with events clipped at the scan end: {int((R.n_events_clipped > 0).sum())} "
          f"({R[R.n_events_clipped > 0][['subject', 'run', 's_past_scan_end']].values.tolist()})")
    print(f"  runs with a non-compliant Intero block: {int((R.n_intero_nc > 0).sum())}; "
          f"runs with no compliant Intero block: {int((R.n_intero == 0).sum())}")
    print(f"  per-run designs, task-regressor VIF (base): median of run max {R.vif_max_task_base.median():.2f}, "
          f"max {R.vif_max_task_base.max():.2f}")
    fr = Jt[[c for c in Jt.columns if c.startswith("full_rank_")]]
    print(f"joint designs: {len(Jt)}; full rank in every fitted variant: {int(fr.fillna(True).all(axis=1).sum())}")
    ok_int = Jt[Jt.subject.isin(s["intero"])].estimable_base.str.contains("C1").all()
    print(f"  C1/C2/C5 estimable in every Intero-set participant (N {len(s['intero'])}): {ok_int}")
    print(f"  C4/C6 estimable in all 54: {Jt.estimable_base.str.contains('C4').all() and Jt.estimable_base.str.contains('C6').all()}")
    print(f"  variants: rvt {int(Jt.variants.str.contains('rvt').sum())}, co2 {int(Jt.variants.str.contains('co2').sum())}")


def e2e(sub):
    """fit() + the loc design vs nilearn FirstLevelModel on one run (base, 6 mm)."""
    from nilearn.glm.first_level import FirstLevelModel
    u = L.run_list(); r = u[u.subject == sub].iloc[0]
    run, n = int(r.run), int(r.n_vol)
    ev, _ = L.read_events(sub, run, n)
    cf = L.read_confounds(sub, run, n)
    flm = FirstLevelModel(t_r=L.TR, slice_time_ref=0.5, hrf_model="spm", drift_model="cosine",
                          high_pass=1 / 200.0, noise_model="ar1", signal_scaling=0, smoothing_fwhm=6,
                          mask_img=g.GROUP_MASK, standardize=False, minimize_memory=False)
    flm.fit(g.bold_path(sub, L.entity(run)), events=ev, confounds=cf)
    dm_flm = flm.design_matrices_[0]
    dm_loc = L.design(ev, cf, n)
    same_dm = list(dm_flm.columns) == list(dm_loc.columns) and np.allclose(dm_flm.values, dm_loc.values)
    out = {}
    mk = flm.masker_
    for c in ("C1", "C5", "C6"):
        spec = L.CONTRASTS[c]
        res = flm.compute_contrast(L.spec_vector(dm_flm, spec), output_type="all")
        rec = {}
        for kind, key in (("effect", "effect_size"), ("variance", "effect_variance"), ("t", "stat")):
            a = mk.transform(res[key]).ravel()
            b = mk.transform(nib.load(L.fl_name(sub, f"run{run}", 6, "base", c, kind))).ravel()
            rec[f"{kind}_max_abs_diff"] = float(np.abs(a - b).max())
            rec[f"{kind}_max_rel_diff"] = float((np.abs(a - b) / np.maximum(np.abs(a), 1e-6)).max())
        out[c] = rec
    js = json.load(open(os.path.join(L.fl_dir(sub), f"{sub}_task-localizer_firstlevel.json")))
    res = dict(subject=sub, run=run, design_identical=bool(same_dm), contrasts=out,
               saved_fits=sorted(js["fits"]), seconds=js["seconds"])
    json.dump(res, open(os.path.join(g.REPO, "Analysis", "localizer_glm_e2e_check.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "designs":
        designs()
    else:
        e2e(sys.argv[2])
