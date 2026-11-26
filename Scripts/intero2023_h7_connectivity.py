#!/usr/bin/env python
"""
intero2023_h7_connectivity.py -- H7 between-network coupling (prereg v0.8 §3 H7, §7 final subsection).

Stage 1, per participant:  python Scripts/intero2023_h7_connectivity.py subject sub-057
  For each run, on the UNSMOOTHED fMRIPrep MNI res-2 BOLD inside the five CAREN networks:
   1. zero-phase Butterworth band-pass 0.01-0.1 Hz (nilearn.signal.butterworth,
      order 5, filtfilt) applied identically to the voxel series and to every regressor;
   2. OLS regression on motion 6+6, FD, CSF, WM (no cosines, no spikes) plus the task
      regressors: the seven Model A conditions (spm + derivative) or, for the FIR
      specification, 2 s bins over 0-32 s for paced / free / rating / subrating;
      +RVT / +PetCO2 specifications add those columns;
   3. network signal = median (or mean) over voxels at each volume; the fMRIPrep
      global_signal is filtered and regressed identically (interpretation gate).
  Windows [onset + shift, onset + duration + shift), shift 4 or 6 s, volume times
  TR (k + 0.5) as in the GLM designs; FD > 0.5 mm volumes dropped; each epoch
  segment demeaned; segments concatenated within condition across runs.
  Writes per-participant coupling (Fisher z) for every pair, condition and
  specification, with the global-signal and five-network-mean partials and the
  temporal SDs, to $OUT/h7/<sub>_h7.csv.

Stage 2:  python Scripts/intero2023_h7_connectivity.py group
  H7a paced vs free (N = 47), H7b aware vs unaware (N = 23, in_S7), SAL-SMN and
  SAL-ECN, two-sided paired t, Holm x2 within each part; the gate and every
  registered check; all 10 pairs FDR (exploratory). Analysis/results/H7_*.csv.
"""
import os, sys, itertools, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
from nilearn.signal import butterworth
from nilearn.glm.first_level import make_first_level_design_matrix
import intero2023_glm_common as g
q = g.q
warnings.simplefilter("ignore")
NETS = ["SAL", "SMN", "ECN", "DMN", "VIS"]
PAIRS_REG = [("SAL", "SMN"), ("SAL", "ECN")]
NUIS = [f"{a}_{b}" for a in ("trans", "rot") for b in "xyz"]
NUIS = NUIS + [c + "_derivative1" for c in NUIS] + ["framewise_displacement", "csf", "white_matter"]
PACED = q.PACED_MOD
COND = dict(paced=PACED, free=["free"], aware=["paced_aware"], unaware=["paced_unaware"])


def fir_design(ev, n_vol):
    e = ev[["onset", "duration", "condition"]].copy()
    e["trial_type"] = np.where(e.condition.isin(PACED), "paced", e.condition)
    e = e[e.trial_type.isin(["paced", "free", "rating", "subrating"])][["onset", "duration", "trial_type"]]
    ft = g.TR * (np.arange(n_vol) + 0.5)
    return make_first_level_design_matrix(ft, events=e, hrf_model="fir", fir_delays=list(range(16)),
                                          drift_model=None)


def run_signals(sub, ent, n_vol, masker, idx, task, physio):
    ev, cf, _ = g.read_run_physio(sub, ent, n_vol, physio)
    if task == "modelA":
        dm = q.model_a(ev, cf[NUIS + physio], n_vol)
    else:
        dm = fir_design(ev, n_vol).drop(columns=["constant"], errors="ignore").reset_index(drop=True)
        dm = pd.concat([dm, cf[NUIS + physio].reset_index(drop=True)], axis=1)
        assert len(dm) == n_vol
    X = dm.drop(columns=["constant"], errors="ignore").values.astype(float)
    Y = masker.transform(g.bold_path(sub, ent)).astype(np.float64)
    gs = pd.read_csv(g.fmriprep_confounds_path(sub, ent), sep="\t")["global_signal"].values[:, None]
    filt = lambda a: butterworth(a, sampling_rate=1 / g.TR, low_pass=0.1, high_pass=0.01, copy=True)
    Xf = np.column_stack([filt(X), np.ones(n_vol)])
    beta, *_ = np.linalg.lstsq(Xf, np.column_stack([filt(Y), filt(gs)]), rcond=None)
    R = np.column_stack([filt(Y), filt(gs)]) - Xf @ beta
    Rv, Rg = R[:, :-1], R[:, -1]
    out = {}
    for summ in ("median", "mean"):
        f = np.median if summ == "median" else np.mean
        out[summ] = np.column_stack([f(Rv[:, idx[n]], axis=1) for n in NETS] + [Rg])
    fd = cf["framewise_displacement"].values
    return out, ev, fd


def windows(ev, n_vol, fd, shift):
    t = g.TR * (np.arange(n_vol) + 0.5)
    segs = {}
    for c, labs in COND.items():
        for e in ev[ev.condition.isin(labs)].itertuples():
            k = np.where((t >= e.onset + shift) & (t < e.onset + e.duration + shift) & (fd <= 0.5))[0]
            if len(k) >= 2:
                segs.setdefault(c, []).append(k)
    return segs


def coupling(S):
    """S: concatenated, epoch-demeaned (T x 6: five networks + global)."""
    rows = {}
    net = S[:, :5]; gsig = S[:, 5]; five = net.mean(axis=1)
    def partial(a, b, c):
        Z = np.column_stack([c, np.ones_like(c)])
        ra = a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]
        rb = b - Z @ np.linalg.lstsq(Z, b, rcond=None)[0]
        return np.corrcoef(ra, rb)[0, 1]
    for (i, a), (j, b) in itertools.combinations(enumerate(NETS), 2):
        rows[f"{a}-{b}"] = dict(r=np.corrcoef(net[:, i], net[:, j])[0, 1],
                                r_gs=partial(net[:, i], net[:, j], gsig),
                                r_five=partial(net[:, i], net[:, j], five))
    sd = {n: net[:, i].std(ddof=1) for i, n in enumerate(NETS)}
    return rows, sd


def subject(sub):
    man = g.run_list()
    runs = man[man.subject == sub]
    union = nib.load(os.path.join(g.REPO, "ROIs", "net_SAL.nii"))
    masks = {n: nib.load(os.path.join(g.REPO, "ROIs", f"net_{n}.nii")).get_fdata() > 0 for n in NETS}
    U = np.any(np.stack(list(masks.values())), axis=0)
    masker = NiftiMasker(mask_img=nib.Nifti1Image(U.astype(np.uint8), union.affine)).fit()
    order = np.flatnonzero(U.ravel(order="C"))
    idx = {n: np.flatnonzero(masks[n].ravel(order="C")[order]) for n in NETS}
    specs = [("modelA", []), ("fir", []), ("modelA", g.PHYSIO_RVT), ("modelA", g.PHYSIO_CO2)]
    rows = []
    for task, physio in specs:
        ok = all((not physio) or (str(r.has_rvt if physio == g.PHYSIO_RVT else r.has_petco2).upper() == "TRUE")
                 for r in runs.itertuples())
        if not ok:
            continue
        sig = {}
        for r in runs.itertuples():
            s, ev, fd = run_signals(sub, r.bids_entity, int(r.n_vol), masker, idx, task, physio)
            sig[r.bids_entity] = (s, ev, fd, int(r.n_vol))
        ptag = "base" if not physio else ("rvt" if physio == g.PHYSIO_RVT else "co2")
        rvt_src = ";".join(sorted(set(runs.rvt_source)))
        for summ in ("median", "mean"):
            for shift in (4.0, 6.0):
                cat = {}
                for ent, (s, ev, fd, n) in sig.items():
                    for c, segs in windows(ev, n, fd, shift).items():
                        for k in segs:
                            seg = s[summ][k]; cat.setdefault(c, []).append(seg - seg.mean(axis=0))
                for c, parts in cat.items():
                    S = np.vstack(parts)
                    cp, sd = coupling(S)
                    for pair, v in cp.items():
                        rows.append(dict(subject=sub, task=task, physio=ptag, rvt_sources=rvt_src, summary=summ,
                                         shift=shift, condition=c, pair=pair, n_vol=len(S), n_epochs=len(parts),
                                         z=np.arctanh(v["r"]), z_gs=np.arctanh(v["r_gs"]),
                                         z_five=np.arctanh(v["r_five"]),
                                         **{f"sd_{n}": sd[n] for n in NETS}))
    od = os.path.join(g.OUT, "h7"); os.makedirs(od, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(od, f"{sub}_h7.csv"), index=False)
    print(f"{sub}: {len(rows)} rows")


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for rank, i in enumerate(o):
        run = max(run, min(1, (m - rank) * p[i])); adj[i] = run
    return adj


def fdr(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); prev = 1
    for rank in range(m - 1, -1, -1):
        i = o[rank]; prev = min(prev, p[i] * m / (rank + 1)); adj[i] = prev
    return adj


def paired(d, a, b, col, subs):
    x = d[(d.condition == a) & d.subject.isin(subs)].set_index("subject")[col]
    y = d[(d.condition == b) & d.subject.isin(subs)].set_index("subject")[col]
    k = x.index.intersection(y.index); diff = (x[k] - y[k]).values
    t, p = stats.ttest_1samp(diff, 0)
    return dict(N=len(k), mean_diff=float(diff.mean()), dz=float(diff.mean() / diff.std(ddof=1)), t=float(t), p=float(p))


def group():
    h3 = pd.read_csv(os.path.join(g.REPO, "Data", "h3_sets.csv"))
    s7 = list(h3[h3.in_S7].subject)
    man = g.run_list(); allsub = sorted(man.subject.unique())
    for c in ("has_rvt", "has_petco2"):
        man[c] = man[c].astype(str).str.upper() == "TRUE"
    gb = man.groupby("subject")
    s2 = list(gb.rvt_source.apply(lambda x: (x == "biopac").all()).pipe(lambda x: x[x]).index)
    s1 = list(gb.has_petco2.all().pipe(lambda x: x[x]).index)
    d = pd.concat([pd.read_csv(os.path.join(g.OUT, "h7", f"{s}_h7.csv")) for s in allsub])
    rd = os.path.join(g.REPO, "Analysis", "results"); os.makedirs(rd, exist_ok=True)
    out = []
    parts = dict(H7a=("paced", "free", allsub), H7b=("aware", "unaware", s7))
    for part, (a, b, subs) in parts.items():
        for task, ptag, summ, shift in itertools.product(("modelA", "fir"), ("base", "rvt", "co2"),
                                                         ("median", "mean"), (4.0, 6.0)):
            if part == "H7b" and task == "fir":
                continue
            if ptag != "base" and (task != "modelA" or summ != "median" or shift != 4.0):
                continue
            use = subs if ptag == "base" else [s for s in subs if s in (s2 if ptag == "rvt" else s1)]
            dd = d[(d.task == task) & (d.physio == ptag) & (d.summary == summ) & (d["shift"] == shift)]
            primary = task == "modelA" and ptag == "base" and summ == "median" and shift == 4.0
            for col in ("z", "z_gs", "z_five"):
                res = []
                for pair in sorted(dd.pair.unique()):
                    r = paired(dd[dd.pair == pair], a, b, col, use); r.update(pair=pair); res.append(r)
                res = pd.DataFrame(res)
                reg = res.pair.isin([f"{x}-{y}" for x, y in PAIRS_REG])
                res["p_holm2"] = np.nan; res.loc[reg, "p_holm2"] = holm(res.loc[reg, "p"])
                res["p_fdr10"] = fdr(res.p)
                res = res.assign(part=part, contrast=f"{a}-{b}", task=task, physio=ptag, summary=summ,
                                 shift=shift, measure=col, primary_spec=primary)
                out.append(res)
        # temporal SD by condition (amplitude check), primary spec
        dd = d[(d.task == "modelA") & (d.physio == "base") & (d.summary == "median") & (d["shift"] == 4.0)
               & (d.pair == "SAL-SMN")]
        for n in NETS:
            r = paired(dd, a, b, f"sd_{n}", subs)
            out.append(pd.DataFrame([dict(part=part, contrast=f"{a}-{b}", measure=f"sd_{n}", pair=n, **r,
                                          task="modelA", physio="base", summary="median", shift=4.0,
                                          primary_spec=True)]))
    # (iv) co-activation: network-mean Model A estimate for the same contrast
    for part, con, subs in (("H7a", "C1", allsub), ("H7b", "C2joint", s7)):
        for n in ("SAL", "SMN", "ECN"):
            mk = NiftiMasker(mask_img=os.path.join(g.REPO, "ROIs", f"net_{n}.nii")).fit()
            v = np.array([mk.transform(os.path.join(g.OUT, "firstlevel", s, f"{s}_var-base_con-{con}_effect.nii.gz")).mean()
                          for s in subs])
            t, p = stats.ttest_1samp(v, 0)
            out.append(pd.DataFrame([dict(part=part, contrast=con, measure="coactivation_mean_effect", pair=n,
                                          N=len(v), mean_diff=float(v.mean()), dz=float(v.mean() / v.std(ddof=1)),
                                          t=float(t), p=float(p), task="modelA", physio="base",
                                          summary="mean", shift=np.nan, primary_spec=True)]))
    res = pd.concat(out, ignore_index=True)
    res.to_csv(os.path.join(rd, "H7_coupling.csv"), index=False)
    print("wrote Analysis/results/H7_coupling.csv")


if __name__ == "__main__":
    if sys.argv[1] == "subject":
        subject(sys.argv[2])
    else:
        group()
