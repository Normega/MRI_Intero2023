#!/usr/bin/env python
"""
intero2023_loc_sheets.py -- first-level nuisance sheets for the localizer runs
(OSF p5ztc, plan §3: "the same nuisance set and FD spike rule" as swc7p).

The breathing-run sheets (FirstLevel/<sub>_task-breathing_run-<k>_desc-model_confounds.tsv)
were written by three R scripts, in this order:
  1. Scripts/intero2023_covariates_linux.R  (port of intero2023_covariates.R): from the
     fMRIPrep confounds TSV keep trans/rot x/y/z, their _derivative1, framewise_displacement,
     csf, white_matter, cosine*, non_steady_state_outlier*; non-finite -> 0; one spike_NN
     column per volume with FD > 0.5 mm; physio columns from Data/regressors;
  2. Scripts/intero2023_physio_sheets.R: physio columns rvt_rrf_z, rvt_z, resp_rate_c
     (from <sub>_<ent>_resp.tsv, padded/trimmed to n_vol, demeaned within run, NA -> 0,
     omitted if < 50% coverage; PMU fallback file if the Biopac rvt_rrf_z covers < 50%)
     and the 4-lag PetCO2 basis (whole or not at all);
  3. Scripts/intero2023_tissue_confounds.R: csf / white_matter re-taken from the fMRIPrep
     derivatives of record.
All three end at the same rule, which this script implements once, in Python, for the
localizer runs. The cosine_* columns are written, as in the breathing sheets, and dropped
at model time (intero2023_loc_common.read_confounds), so the sheet construction itself is
identical and can be checked against the breathing sheets:

  python Scripts/intero2023_loc_sheets.py validate   # rebuild all 95 breathing sheets in
                                                      # memory, compare with the committed ones
  python Scripts/intero2023_loc_sheets.py build      # write FirstLevel/localizer/*.tsv and
                                                      # FirstLevel/localizer/localizer_manifest.csv

No BOLD voxel is read (the NIfTI header gives n_vol only).
"""
import os, sys, glob
import numpy as np, pandas as pd, nibabel as nib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
g = L.g

MOTION = [f"{a}_{b}" for a in ("trans", "rot") for b in ("x", "y", "z")]
BELT = {"rvt_rrf_z": "rvt_rrf_z", "rvt_z": "rvt_z", "resp_rate_c": "resp_rate_brpm"}


def demean0(x):
    x = np.asarray(x, dtype=float).copy(); ok = np.isfinite(x)
    if not ok.any():
        return np.zeros(len(x))
    x[ok] -= x[ok].mean(); x[~ok] = 0.0
    return x


def fit_len(x, n):
    x = np.asarray(x, dtype=float)
    return x[:n] if len(x) >= n else np.concatenate([x, np.full(n - len(x), np.nan)])


def read_reg(f):
    return pd.read_csv(f, sep="\t", na_values=["n/a"])


def build_sheet(sub, ent):
    """Returns (sheet DataFrame, info dict) for one run, by the breathing-sheet rule."""
    fp = os.path.join(g.FMRIPREP, sub, "func", f"{sub}_{ent}_desc-confounds_timeseries.tsv")
    d = pd.read_csv(fp, sep="\t", na_values=["n/a"])
    n_vol = len(d)
    keep = MOTION + [m + "_derivative1" for m in MOTION] + ["framewise_displacement", "csf", "white_matter"] \
        + [c for c in d.columns if c.startswith("cosine")] \
        + [c for c in d.columns if c.startswith("non_steady_state_outlier")]
    keep = [c for c in keep if c in d.columns]
    cf = d[keep].apply(pd.to_numeric, errors="coerce")
    cf = cf.where(np.isfinite(cf), 0.0)
    fd = cf.framewise_displacement.values
    spikes = np.where(fd > L.FD_SPIKE)[0]
    for k, v in enumerate(spikes, 1):
        col = np.zeros(n_vol, dtype=int); col[v] = 1
        cf[f"spike_{k:02d}"] = col
    # physio (intero2023_physio_sheets.R physio_block)
    bf = os.path.join(L.REGDIR, f"{sub}_{ent}_resp.tsv")
    pf = os.path.join(L.REGDIR, f"{sub}_{ent}_desc-pmu_resp.tsv")
    bio_ok = os.path.exists(bf) and np.isfinite(fit_len(read_reg(bf).rvt_rrf_z, n_vol)).mean() >= L.PHYSIO_MIN_COV
    rf = bf if bio_ok else (pf if os.path.exists(pf) else (bf if os.path.exists(bf) else None))
    has_rvt, rvt_cov = False, 0.0
    if rf is not None:
        r = read_reg(rf)
        for nm, src in BELT.items():
            v = fit_len(r[src], n_vol)
            if nm == "rvt_rrf_z":
                rvt_cov = float(np.isfinite(v).mean())
            if np.isfinite(v).mean() >= L.PHYSIO_MIN_COV:
                cf[nm] = demean0(v); has_rvt = True
    co2f = os.path.join(L.REGDIR, f"{sub}_{ent}_petco2.tsv")
    has_co2 = False
    if os.path.exists(co2f):
        p = read_reg(co2f)
        cc = [c for c in p.columns if c.startswith("petco2_hrf_lag")]
        if cc and all(np.isfinite(fit_len(p[c], n_vol)).mean() >= L.PHYSIO_MIN_COV for c in cc):
            for c in cc:
                cf[c] = demean0(fit_len(p[c], n_vol))
            has_co2 = True
    cf = cf.round(6)
    src = "none" if not has_rvt else ("pmu" if rf == pf else "biopac")
    info = dict(n_vol=n_vol, fd_mean=float(fd.mean()), fd_frac_gt_thresh=float((fd > L.FD_SPIKE).mean()),
                n_spikes=len(spikes), has_rvt=has_rvt, rvt_rrf_coverage=round(rvt_cov, 3),
                has_rvt_rrf="rvt_rrf_z" in cf.columns, rvt_source=src, has_petco2=has_co2,
                n_nss=sum(c.startswith("non_steady_state_outlier") for c in cf.columns))
    return cf, info


def validate():
    """Rebuild every breathing sheet and compare with the committed one, column by column."""
    rows = []
    for f in sorted(glob.glob(os.path.join(g.REPO, "FirstLevel", "*_task-breathing_run-*_desc-model_confounds.tsv"))):
        b = os.path.basename(f); sub = b.split("_")[0]; ent = b[len(sub) + 1:].split("_desc-")[0]
        old = pd.read_csv(f, sep="\t", na_values=["n/a"])
        new, _ = build_sheet(sub, ent)
        same_cols = list(old.columns) == list(new.columns)
        common = [c for c in old.columns if c in new.columns]
        diff = max(float(np.nanmax(np.abs(old[c].values - new[c].values))) for c in common)
        worst = max(common, key=lambda c: float(np.nanmax(np.abs(old[c].values - new[c].values))))
        rows.append(dict(subject=sub, entity=ent, n_rows_old=len(old), n_rows_new=len(new),
                         same_columns=same_cols, only_old=";".join(sorted(set(old.columns) - set(new.columns))),
                         only_new=";".join(sorted(set(new.columns) - set(old.columns))),
                         max_abs_diff=diff, worst_column=worst))
    d = pd.DataFrame(rows)
    out = os.path.join(g.REPO, "Analysis", "localizer_sheet_validation.csv")
    d.to_csv(out, index=False)
    ok = d.same_columns & (d.n_rows_old == d.n_rows_new) & (d.max_abs_diff <= 1.5e-6)
    print(f"breathing sheets rebuilt: {len(d)}; identical columns and values (<= 1.5e-6): {int(ok.sum())}")
    if (~ok).any():
        print(d[~ok].to_string())
    print(f"wrote {out}")
    return ok.all()


def build():
    os.makedirs(L.SHEETS, exist_ok=True)
    eq = pd.read_csv(os.path.join(L.QCDIR, "events_qc.csv"))
    mg = pd.read_csv(os.path.join(L.QCDIR, "motion_gate.csv"))
    comp = pd.read_csv(os.path.join(L.QCDIR, "compliance.csv"))
    rows = []
    for r in eq.itertuples():
        sub, run = r.subject, int(r.run); ent = L.entity(run)
        cf, info = build_sheet(sub, ent)
        hdr = int(nib.load(g.bold_path(sub, ent)).shape[3])           # header only
        cf.to_csv(L.sheet_path(sub, run), sep="\t", index=False, na_rep="n/a")
        m = mg[(mg.subject == sub) & (mg.run == run)].iloc[0]
        gate = bool(info["fd_mean"] > L.FD_SPIKE or info["fd_frac_gt_thresh"] > L.FD_GATE_FRAC)
        assert gate == bool(m.motion_gate_fail), (sub, run)
        assert info["n_spikes"] == int(m.n_spikes), (sub, run)
        partial = r.recommended_use == "exclude_partial_log"
        cb = comp[(comp.subject == sub) & (comp.run == run)]
        reason = ";".join([x for x, f in (("partial_task", partial), ("motion_gate", gate)) if f])
        rows.append(dict(subject=sub, run=run, bids_entity=ent, protocol_version=r.protocol_version,
                         n_vol=info["n_vol"], n_vol_header=hdr, scan_dur_s=info["n_vol"] * L.TR,
                         recommended_use=r.recommended_use, s_lost_past_scan_end=r.s_lost_past_scan_end,
                         fd_mean=round(info["fd_mean"], 4), fd_frac_gt_thresh=round(info["fd_frac_gt_thresh"], 4),
                         n_spikes=info["n_spikes"], n_nss=info["n_nss"], motion_gate_fail=gate,
                         partial_task=partial, usable=not (gate or partial), exclusion=reason or "",
                         n_intero_compliant=int(cb.compliant.sum()), n_intero_nc=int((~cb.compliant.astype(bool)).sum()),
                         has_rvt=bool(info["has_rvt_rrf"]), rvt_rrf_coverage=info["rvt_rrf_coverage"],
                         rvt_source=info["rvt_source"], has_petco2=info["has_petco2"]))
    m = pd.DataFrame(rows)
    m.to_csv(L.MANIFEST, index=False)
    u = m[m.usable]
    print(f"sheets written: {len(m)}; usable runs {len(u)}, participants {u.subject.nunique()}")
    print("excluded:\n" + m[~m.usable][["subject", "run", "exclusion"]].to_string(index=False))
    print(f"usable runs with rvt_rrf_z: {int(u.has_rvt.sum())}; with PetCO2 basis: {int(u.has_petco2.sum())}; "
          f"rvt_source {u.rvt_source.value_counts().to_dict()}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "build"
    if cmd == "validate":
        sys.exit(0 if validate() else 1)
    build()
