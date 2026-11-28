#!/usr/bin/env python
"""
intero2023_loc_connectivity.py -- localizer stage C2 (EXPLORATORY; OSF p5ztc plan §6, §8 T5):
network functional connectivity (FC) in the localizer runs, background FC, and T5
(fROI-node coupling in the breathing runs). gPPI is in intero2023_loc_gppi.py.
Choices: Analysis/localizer_stageC2_implementation_notes.md (committed before any
connectivity statistic was computed). Readings: Preregistration/localizer_deviations.md.

  python Scripts/intero2023_loc_connectivity.py subject sub-060   # localizer network FC, one participant
  python Scripts/intero2023_loc_connectivity.py t5subject sub-060 # breathing-run fROI-node coupling (T5)
  python Scripts/intero2023_loc_connectivity.py verify            # T5 network pairs == H7 per-subject CSVs;
                                                                  # window / volume counts (no statistic)
  python Scripts/intero2023_loc_connectivity.py group             # tests -> Analysis/results/EXPLORATORY_LOC_C2_*
  python Scripts/intero2023_loc_connectivity.py figure            # Figures/localizer/loc_C2_network_fc.png

Everything is the swc7p H7 method (Scripts/intero2023_h7_connectivity.py, imported, never
modified): unsmoothed fMRIPrep MNI res-2 BOLD; zero-phase Butterworth 0.01-0.1 Hz on voxels
and regressors; OLS on H7's nuisance set (motion 6+6, FD, CSF, WM; no cosines, no spikes);
median network signal; windows [onset + shift, onset + duration + shift) at volume times
TR (k + 0.5); FD > 0.5 mm volumes dropped; each segment demeaned; concatenated within
condition; Pearson r -> Fisher z; fMRIPrep global_signal and the five-network mean,
processed identically, partialled out (H7.coupling).
Localizer specifications:
  task = total       nuisance only (network FC, plan §6 first bullet; reading logged)
  task = background  + localizer task regressors (intero, extero, match, self, instruction,
                     intero_nc; Rest implicit), SPM canonical + temporal derivative as H7
  physio = base | rvt (+rvt_rrf_z) | co2 (+4-lag PetCO2) | nocsf (H7 nuisance without csf; added)
Per-participant rows: $OUT/loc_conn/<sub>_locfc.csv and <sub>_t5.csv (outside the repo).
"""
import os, sys, json, itertools, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
from nilearn.signal import butterworth
from nilearn.glm.first_level import make_first_level_design_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
import intero2023_h7_connectivity as H7
g = L.g
q = g.q
warnings.simplefilter("ignore")
NETS = H7.NETS                                   # SAL SMN ECN DMN VIS (H7 order)
NET_PAIRS = [f"{a}-{b}" for a, b in itertools.combinations(NETS, 2)]
LCOND = ["intero", "extero", "match", "self", "rest"]
CONTRASTS = {"C1": ("intero", "extero", "intero"), "C2": ("intero", "self", "intero"),
             "C4": ("match", "extero", "all"), "C5": ("match", "intero", "intero")}
PRIMARY = {"C1": ["SAL-SMN", "SAL-DMN"], "C2": ["SAL-DMN"]}   # plan §6 (DMN-SAL = H7's SAL-DMN label)
MIN_EPOCHS, MIN_VOLS = 2, 20
CD = os.path.join(g.OUT, "loc_conn")
RD = os.path.join(g.REPO, "Analysis", "results")
FIG = os.path.join(g.REPO, "Figures", "localizer")
PARCELS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere"]
SHORT = {"caren_AInsL": "AInsL", "caren_AInsR": "AInsR", "caren_dACC": "dACC", "caren_PIns": "PIns",
         "ieat_ACCsphere": "sphere"}
FILT = lambda a: butterworth(a, sampling_rate=1 / g.TR, low_pass=0.1, high_pass=0.01, copy=True)


# --------------------------------------------------------------------------- shared helpers
def denoise(Y, X, gs):
    """H7.run_signals step 2, verbatim: filter data, regressors and global signal; OLS on
    [filtered X, intercept]; return residual voxel series and residual global signal."""
    n = Y.shape[0]
    Xf = np.column_stack([FILT(X), np.ones(n)]) if X.shape[1] else np.ones((n, 1))
    D = np.column_stack([FILT(Y), FILT(gs)])
    beta, *_ = np.linalg.lstsq(Xf, D, rcond=None)
    R = D - Xf @ beta
    return R[:, :-1], R[:, -1]


def summarise(Rv, idx, nodes, summ):
    f = np.median if summ == "median" else np.mean
    return np.column_stack([f(Rv[:, idx[n]], axis=1) for n in nodes])


def partial_r(a, b, c):
    Z = np.column_stack([c, np.ones_like(c)])
    ra = a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]
    rb = b - Z @ np.linalg.lstsq(Z, b, rcond=None)[0]
    return np.corrcoef(ra, rb)[0, 1]


def union_masker(masks):
    ref = nib.load(g.GROUP_MASK)
    U = np.any(np.stack(list(masks.values())), axis=0)
    masker = NiftiMasker(mask_img=nib.Nifti1Image(U.astype(np.uint8), ref.affine)).fit()
    order = np.flatnonzero(U.ravel(order="C"))
    idx = {n: np.flatnonzero(m.ravel(order="C")[order]) for n, m in masks.items()}
    return masker, idx


def net_masks():
    return {n: nib.load(os.path.join(g.REPO, "ROIs", f"net_{n}.nii")).get_fdata() > 0 for n in NETS}


# --------------------------------------------------------------------------- localizer FC
def loc_events(sub, run, n_vol):
    """Every block, Rest included, with the registered coding (non-compliant Intero ->
    intero_nc) and clipping. Without Rest it equals L.read_events (asserted)."""
    ev = pd.read_csv(L.events_path(sub, run), sep="\t")
    ev = ev[["onset", "duration", "trial_type", "n_presses"]].copy()
    it = ev.trial_type == "intero"
    ev.loc[it & (pd.to_numeric(ev.n_presses, errors="coerce") < L.MIN_PRESSES), "trial_type"] = "intero_nc"
    end = n_vol * g.TR
    over = (ev.onset + ev.duration) > end
    ev.loc[over, "duration"] = (end - ev.loc[over, "onset"]).clip(lower=0)
    ev = ev[ev.duration > 0][["onset", "duration", "trial_type"]].reset_index(drop=True)
    model, _ = L.read_events(sub, run, n_vol)
    chk = ev[ev.trial_type != "rest"].reset_index(drop=True)
    assert len(chk) == len(model) and (chk.trial_type.values == model.trial_type.values).all() \
        and np.allclose(chk.onset, model.onset) and np.allclose(chk.duration, model.duration), (sub, run)
    return ev, model


def loc_windows(ev, n_vol, fd, shift):
    t = g.TR * (np.arange(n_vol) + 0.5)
    segs = {}
    for e in ev[ev.trial_type.isin(LCOND)].itertuples():
        k = np.where((t >= e.onset + shift) & (t < e.onset + e.duration + shift) & (fd <= 0.5))[0]
        if len(k) >= 2:
            segs.setdefault(e.trial_type, []).append(k)
    return segs


SPECS = [("total", "base"), ("background", "base"), ("total", "rvt"), ("background", "rvt"),
         ("total", "co2"), ("background", "co2"), ("total", "nocsf"), ("background", "nocsf")]
PHYSIO_COLS = {"base": [], "rvt": L.PHYSIO_RVT, "co2": L.PHYSIO_CO2, "nocsf": []}


def subject(sub):
    runs = L.run_list(); runs = runs[runs.subject == sub]
    masker, idx = union_masker(net_masks())
    cache = {}
    for r in runs.itertuples():
        n = int(r.n_vol)
        Y = masker.transform(g.bold_path(sub, L.entity(r.run))).astype(np.float64)
        assert Y.shape[0] == n
        gs = pd.read_csv(g.fmriprep_confounds_path(sub, L.entity(r.run)), sep="\t")["global_signal"].values[:, None]
        ev, model = loc_events(sub, r.run, n)
        cache[int(r.run)] = (Y, gs, ev, model, n)
    rows = []
    for task, ptag in SPECS:
        phys = PHYSIO_COLS[ptag]
        if ptag == "rvt" and not runs.has_rvt.all():
            continue
        if ptag == "co2" and not runs.has_petco2.all():
            continue
        sig = {}
        for run, (Y, gs, ev, model, n) in cache.items():
            cf = L.read_confounds(sub, run, n, phys)
            nuis = [c for c in H7.NUIS if not (ptag == "nocsf" and c == "csf")] + list(phys)
            C = cf[nuis]
            if task == "background":
                ft = g.TR * (np.arange(n) + 0.5)
                dm = make_first_level_design_matrix(ft, events=model, hrf_model="spm + derivative", drift_model=None,
                                                    add_regs=C.values, add_reg_names=nuis)
                X = dm.drop(columns=["constant"]).values.astype(float)
            else:
                X = C.values.astype(float)
            Rv, Rg = denoise(Y, X, gs)
            fd = cf["framewise_displacement"].values
            sig[run] = ({s: np.column_stack([summarise(Rv, idx, NETS, s), Rg]) for s in ("median", "mean")}, ev, fd, n)
        for summ in ("median", "mean"):
            for shift in (4.0, 6.0):
                cat = {}
                for run, (s, ev, fd, n) in sig.items():
                    for c, segs in loc_windows(ev, n, fd, shift).items():
                        for k in segs:
                            seg = s[summ][k]; cat.setdefault(c, []).append(seg - seg.mean(axis=0))
                for c, parts in cat.items():
                    S = np.vstack(parts)
                    cp, sd = H7.coupling(S)
                    for pair, v in cp.items():
                        rows.append(dict(subject=sub, task=task, physio=ptag, summary=summ, shift=shift, condition=c,
                                         pair=pair, n_vol=len(S), n_epochs=len(parts), z=np.arctanh(v["r"]),
                                         z_gs=np.arctanh(v["r_gs"]), z_five=np.arctanh(v["r_five"]),
                                         **{f"sd_{k}": sd[k] for k in NETS}))
    os.makedirs(CD, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(CD, f"{sub}_locfc.csv"), index=False)
    print(f"{sub}: {len(rows)} rows", flush=True)


# --------------------------------------------------------------------------- T5 (breathing runs)
def gate_primary():
    d = pd.read_csv(os.path.join(RD, "LOC_L4_gate.csv"))
    d = d[d.rule == "C1"].set_index("parcel").primary_froi
    out = {p: d[p] for p in PARCELS}
    assert out == {"caren_AInsL": "C1", "caren_AInsR": "C1loso", "caren_dACC": "C1", "caren_PIns": "C1",
                   "ieat_ACCsphere": "C1loso"}, out
    return out


def froi_masks(sub):
    """Gate-primary 10% fROI per parcel as 3-D boolean (group-mask grid)."""
    z = np.load(os.path.join(g.OUT, "froi", "froi_masks.npz"))
    rule = gate_primary(); out = {}
    for p in PARCELS:
        pm = nib.load(os.path.join(g.REPO, "ROIs", f"{p}.nii")).get_fdata() > 0
        v = np.zeros(pm.shape, bool); v[pm] = z[f"{sub}|{rule[p]}|0.10|{p}"]
        out[SHORT[p]] = v
    return out


FR = [SHORT[p] for p in PARCELS]
T5_PAIRS = [f"{a}-{b}" for a in FR for b in NETS] + [f"{a}-{b}" for a, b in itertools.combinations(FR, 2)]


def t5_subject(sub):
    man = g.run_list(); runs = man[man.subject == sub]
    assert len(runs), sub
    nm = net_masks(); fm = froi_masks(sub)
    fu = np.any(np.stack(list(fm.values())), axis=0)
    masks = dict(nm)                                              # full networks (H7 check; five-net mean)
    masks.update({f"{n}x": nm[n] & ~fu for n in NETS})             # networks without the fROI voxels
    masks.update(fm)
    ov = fm["dACC"] & fm["sphere"]                                 # dACC parcel and the sphere overlap
    masks["dACCx"] = fm["dACC"] & ~ov; masks["spherex"] = fm["sphere"] & ~ov
    masker, idx = union_masker(masks)
    nodes = list(masks)
    info = dict(subject=sub, overlap_dACC_sphere=int(ov.sum()), **{f"nvox_{k}": int(v.sum()) for k, v in masks.items()})
    specs = [("modelA", []), ("modelA", g.PHYSIO_RVT), ("modelA", g.PHYSIO_CO2)]
    rows = []
    for task, physio in specs:
        ok = all((not physio) or (str(r.has_rvt if physio == g.PHYSIO_RVT else r.has_petco2).upper() == "TRUE")
                 for r in runs.itertuples())
        if not ok:
            continue
        ptag = "base" if not physio else ("rvt" if physio == g.PHYSIO_RVT else "co2")
        sig = {}
        for r in runs.itertuples():
            n = int(r.n_vol)
            ev, cf, _ = g.read_run_physio(sub, r.bids_entity, n, physio)
            dm = q.model_a(ev, cf[H7.NUIS + physio], n)                       # as H7.run_signals
            X = dm.drop(columns=["constant"], errors="ignore").values.astype(float)
            Y = masker.transform(g.bold_path(sub, r.bids_entity)).astype(np.float64)
            gs = pd.read_csv(g.fmriprep_confounds_path(sub, r.bids_entity), sep="\t")["global_signal"].values[:, None]
            Rv, Rg = denoise(Y, X, gs)
            sig[r.bids_entity] = ({s: (summarise(Rv, idx, nodes, s), Rg) for s in ("median", "mean")}, ev,
                                  cf["framewise_displacement"].values, n)
        for summ in ("median", "mean"):
            for shift in (4.0, 6.0):
                cat = {}
                for ent, (s, ev, fd, n) in sig.items():
                    for c, segs in H7.windows(ev, n, fd, shift).items():
                        for k in segs:
                            A, gg = s[summ]; seg = np.column_stack([A[k], gg[k]])
                            cat.setdefault(c, []).append(seg - seg.mean(axis=0))
                for c, parts in cat.items():
                    S = np.vstack(parts); col = {nd: S[:, i] for i, nd in enumerate(nodes)}; gsig = S[:, -1]
                    five = np.column_stack([col[n] for n in NETS]).mean(axis=1)
                    pairs = [(f"{a}-{b}", col[a], col[b]) for a, b in itertools.combinations(NETS, 2)]
                    pairs += [(f"{a}-{b}", col[a], col[f"{b}x"]) for a in FR for b in NETS]
                    for a, b in itertools.combinations(FR, 2):
                        ka, kb = (f"{a}x", f"{b}x") if {a, b} == {"dACC", "sphere"} else (a, b)
                        pairs.append((f"{a}-{b}", col[ka], col[kb]))
                    for pair, x, y in pairs:
                        rows.append(dict(subject=sub, task=task, physio=ptag, summary=summ, shift=shift, condition=c,
                                         pair=pair, kind="network" if pair in NET_PAIRS else "froi",
                                         n_vol=len(S), n_epochs=len(parts), z=np.arctanh(np.corrcoef(x, y)[0, 1]),
                                         z_gs=np.arctanh(partial_r(x, y, gsig)), z_five=np.arctanh(partial_r(x, y, five))))
                    for nd in FR:
                        rows.append(dict(subject=sub, task=task, physio=ptag, summary=summ, shift=shift, condition=c,
                                         pair=nd, kind="sd", n_vol=len(S), n_epochs=len(parts), z=np.nan,
                                         z_gs=np.nan, z_five=np.nan, sd=col[nd].std(ddof=1)))
    os.makedirs(CD, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(CD, f"{sub}_t5.csv"), index=False)
    json.dump(info, open(os.path.join(CD, f"{sub}_t5_nodes.json"), "w"), indent=1)
    print(f"{sub}: T5 {len(rows)} rows", flush=True)


def t5_sets():
    fr = json.load(open(os.path.join(g.OUT, "froi", "froi_index.json")))["subjects"]
    man = g.run_list()
    h3 = pd.read_csv(os.path.join(g.REPO, "Data", "h3_sets.csv"))
    for c in ("has_rvt", "has_petco2"):
        man[c] = man[c].astype(str).str.upper() == "TRUE"
    gb = man.groupby("subject")
    s2 = set(gb.rvt_source.apply(lambda x: (x == "biopac").all()).pipe(lambda x: x[x]).index)   # as H7 group()
    s1 = set(gb.has_petco2.all().pipe(lambda x: x[x]).index)
    a = sorted(set(fr) & set(man.subject)); b = sorted(set(a) & set(h3[h3.in_S7].subject))
    return dict(T5a=a, T5b=b, S2=s2, S1=s1)


# --------------------------------------------------------------------------- verification (no statistic)
def verify():
    """(1) T5 network pairs (full networks) reproduce the H7 per-subject Fisher z exactly;
    (2) volume / window counts per localizer condition."""
    S = t5_sets(); out = {}
    diffs = []
    for s in S["T5a"]:
        a = pd.read_csv(os.path.join(CD, f"{s}_t5.csv")); a = a[a.kind == "network"]
        h = pd.read_csv(os.path.join(g.OUT, "h7", f"{s}_h7.csv")); h = h[h.task == "modelA"]
        m = a.merge(h, on=["physio", "summary", "shift", "condition", "pair"], suffixes=("", "_h7"))
        assert len(m) == len(a), (s, len(m), len(a))
        diffs.append(np.abs(m[["z", "z_gs", "z_five"]].values - m[["z_h7", "z_gs_h7", "z_five_h7"]].values).max())
        assert (m.n_vol == m.n_vol_h7).all()
    out["t5_vs_h7_max_abs_diff"] = float(np.max(diffs)); out["t5_vs_h7_subjects"] = len(diffs)
    s, _, _ = L.participant_sets()
    d = pd.concat([pd.read_csv(os.path.join(CD, f"{x}_locfc.csv")) for x in s["all"]])
    d = d[(d.task == "total") & (d.physio == "base") & (d.summary == "median") & (d["shift"] == 4.0) & (d.pair == "SAL-SMN")]
    vc = d.groupby("condition").agg(n_subjects=("subject", "nunique"), vols_median=("n_vol", "median"),
                                    vols_min=("n_vol", "min"), vols_max=("n_vol", "max"),
                                    epochs_median=("n_epochs", "median"), epochs_min=("n_epochs", "min"))
    out["localizer_volumes_per_condition"] = vc.reset_index().to_dict(orient="records")
    nodes = pd.DataFrame([json.load(open(os.path.join(CD, f"{x}_t5_nodes.json"))) for x in S["T5a"]])
    out["t5_node_voxels_median"] = nodes.drop(columns="subject").median().to_dict()
    out["t5_overlap_dACC_sphere_n_subjects_gt0"] = int((nodes.overlap_dACC_sphere > 0).sum())
    json.dump(out, open(os.path.join(g.REPO, "Analysis", "localizer_stageC2_verification.json"), "w"), indent=1)
    print(json.dumps(out, indent=1, default=str))


# --------------------------------------------------------------------------- group
def holm(p):
    return H7.holm(p)


def fdr(p):
    return H7.fdr(p)


def paired(d, a, b, col, subs):
    x = d[(d.condition == a) & d.subject.isin(subs)].set_index("subject")[col]
    y = d[(d.condition == b) & d.subject.isin(subs)].set_index("subject")[col]
    k = x.index.intersection(y.index); diff = (x[k] - y[k]).values
    t, p = stats.ttest_1samp(diff, 0)
    return dict(N=len(k), mean_a=float(x[k].mean()), mean_b=float(y[k].mean()), mean_diff=float(diff.mean()),
                dz=float(diff.mean() / diff.std(ddof=1)), t=float(t), p=float(p))


def eligible(d, a, b, subs):
    """Participants with >= MIN_EPOCHS windows and >= MIN_VOLS volumes in both conditions."""
    e = d[(d.pair == d.pair.iloc[0])].set_index(["subject", "condition"])
    ok = []
    for s in subs:
        if all((s, c) in e.index and e.loc[(s, c), "n_epochs"] >= MIN_EPOCHS and e.loc[(s, c), "n_vol"] >= MIN_VOLS
               for c in (a, b)):
            ok.append(s)
    return ok


def group():
    s, _, _ = L.participant_sets()
    d = pd.concat([pd.read_csv(os.path.join(CD, f"{x}_locfc.csv")) for x in s["all"]], ignore_index=True)
    out, cm = [], []
    for task in ("total", "background"):
        for con, (a, b, sset) in CONTRASTS.items():
            for ptag, summ, shift in [("base", m, sh) for m in ("median", "mean") for sh in (4.0, 6.0)] + \
                    [("rvt", "median", 4.0), ("co2", "median", 4.0), ("nocsf", "median", 4.0)]:
                base_set = s[sset]
                if ptag in ("rvt", "co2"):
                    base_set = [x for x in base_set if x in s[ptag]]
                comps = [(ptag, base_set)] + ([("base", base_set)] if ptag in ("rvt", "co2") else [])
                for pt, subs in comps:
                    dd = d[(d.task == task) & (d.physio == pt) & (d.summary == summ) & (d["shift"] == shift)]
                    use = eligible(dd, a, b, subs)
                    label = pt if ptag not in ("rvt", "co2") or pt != "base" else f"base_in_{ptag}_set"
                    primary = pt == "base" and label == "base" and summ == "median" and shift == 4.0
                    for col in ("z", "z_gs", "z_five"):
                        res = []
                        for pair in NET_PAIRS:
                            r = paired(dd[dd.pair == pair], a, b, col, use); r.update(pair=pair); res.append(r)
                        res = pd.DataFrame(res)
                        prim = res.pair.isin(PRIMARY.get(con, []))
                        res["primary_pair"] = prim
                        res["p_holm_primary"] = np.nan
                        if prim.any():
                            res.loc[prim, "p_holm_primary"] = holm(res.loc[prim, "p"])
                        res["p_fdr10"] = fdr(res.p)
                        res = res.assign(analysis=task, contrast=con, cond_a=a, cond_b=b, physio=label, summary=summ,
                                         shift=shift, measure=col, primary_spec=primary)
                        out.append(res)
                    if summ == "median" and shift == 4.0:
                        sd = dd[dd.pair == "SAL-SMN"]
                        for n in NETS:
                            r = paired(sd, a, b, f"sd_{n}", use)
                            out.append(pd.DataFrame([dict(analysis=task, contrast=con, cond_a=a, cond_b=b, physio=label,
                                                          summary=summ, shift=shift, measure="sd", pair=n, **r,
                                                          primary_spec=primary)]))
    for task in ("total", "background"):
        dd = d[(d.task == task) & (d.physio == "base") & (d.summary == "median") & (d["shift"] == 4.0)]
        for (c, pair), x in dd.groupby(["condition", "pair"]):
            for col in ("z", "z_gs"):
                cm.append(dict(analysis=task, condition=c, pair=pair, measure=col, N=x.subject.nunique(),
                               mean=float(x[col].mean()), sd=float(x[col].std(ddof=1))))
        sdv = dd[dd.pair == "SAL-SMN"]
        for c, x in sdv.groupby("condition"):
            for n in NETS:
                cm.append(dict(analysis=task, condition=c, pair=n, measure="sd", N=x.subject.nunique(),
                               mean=float(x[f"sd_{n}"].mean()), sd=float(x[f"sd_{n}"].std(ddof=1))))
        vc = dd[dd.pair == "SAL-SMN"].groupby("condition")
        for c, x in vc:
            cm.append(dict(analysis=task, condition=c, pair="", measure="n_vol", N=x.subject.nunique(),
                           mean=float(x.n_vol.mean()), sd=float(x.n_vol.std(ddof=1)), median=float(x.n_vol.median()),
                           min=float(x.n_vol.min())))
    # co-activation (H7 check iv): network-mean 6 mm joint base effect of the same contrast
    from concurrent.futures import ThreadPoolExecutor
    mk = {n: NiftiMasker(mask_img=os.path.join(g.REPO, "ROIs", f"net_{n}.nii")).fit() for n in NETS}
    for con, (a, b, sset) in CONTRASTS.items():
        subs = s[sset]
        def nm(x):
            img = nib.load(L.fl_name(x, "joint", 6, "base", con, "effect"))
            return [float(mk[n].transform(img).mean()) for n in NETS]
        with ThreadPoolExecutor(8) as ex:
            V = np.array(list(ex.map(nm, subs)))
        for j, n in enumerate(NETS):
            v = V[:, j]; t, p = stats.ttest_1samp(v, 0)
            out.append(pd.DataFrame([dict(analysis="coactivation", contrast=con, cond_a=a, cond_b=b, physio="base",
                                          summary="mean", shift=np.nan, measure="coactivation_mean_effect", pair=n,
                                          N=len(v), mean_diff=float(v.mean()), dz=float(v.mean() / v.std(ddof=1)),
                                          t=float(t), p=float(p), primary_spec=True)]))
    res = pd.concat(out, ignore_index=True)
    res.to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_netfc.csv"), index=False)
    pd.DataFrame(cm).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_netfc_condition_means.csv"), index=False)
    print("wrote EXPLORATORY_LOC_C2_netfc.csv", flush=True)
    t5_group()


def t5_group():
    S = t5_sets()
    d = pd.concat([pd.read_csv(os.path.join(CD, f"{x}_t5.csv")) for x in S["T5a"]], ignore_index=True)
    out = []
    parts = dict(T5a=("paced", "free", S["T5a"]), T5b=("aware", "unaware", S["T5b"]))
    for part, (a, b, subs) in parts.items():
        for ptag, summ, shift in [("base", m, sh) for m in ("median", "mean") for sh in (4.0, 6.0)] + \
                [("rvt", "median", 4.0), ("co2", "median", 4.0)]:
            ss = subs if ptag == "base" else [x for x in subs if x in (S["S2"] if ptag == "rvt" else S["S1"])]
            comps = [(ptag, ss)] + ([("base", ss)] if ptag != "base" else [])
            for pt, use in comps:
                dd = d[(d.physio == pt) & (d.summary == summ) & (d["shift"] == shift)]
                label = pt if ptag == "base" else (ptag if pt != "base" else f"base_in_{ptag}_set")
                primary = label == "base" and summ == "median" and shift == 4.0
                for col in ("z", "z_gs", "z_five"):
                    for kind, plist in (("froi", T5_PAIRS), ("network", NET_PAIRS)):
                        res = []
                        for pair in plist:
                            r = paired(dd[dd.pair == pair], a, b, col, use); r.update(pair=pair); res.append(r)
                        res = pd.DataFrame(res)
                        res["p_fdr"] = fdr(res.p); res["fdr_family"] = len(plist)
                        out.append(res.assign(part=part, contrast=f"{a}-{b}", kind=kind, physio=label, summary=summ,
                                              shift=shift, measure=col, primary_spec=primary))
                if summ == "median" and shift == 4.0:
                    res = []
                    for nd in FR:
                        r = paired(dd[dd.pair == nd], a, b, "sd", use); r.update(pair=nd); res.append(r)
                    out.append(pd.DataFrame(res).assign(part=part, contrast=f"{a}-{b}", kind="sd", physio=label,
                                                        summary=summ, shift=shift, measure="sd", primary_spec=primary))
    res = pd.concat(out, ignore_index=True)
    res.to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_T5.csv"), index=False)
    # H7 network-level results in the same participants, from the registered H7 per-subject files
    h = pd.concat([pd.read_csv(os.path.join(g.OUT, "h7", f"{x}_h7.csv")) for x in S["T5a"]], ignore_index=True)
    h = h[(h.task == "modelA") & (h.physio == "base") & (h.summary == "median") & (h["shift"] == 4.0)]
    hh = []
    for part, (a, b, subs) in parts.items():
        for col in ("z", "z_gs"):
            r_ = []
            for pair in NET_PAIRS:
                r = paired(h[h.pair == pair], a, b, col, subs); r.update(pair=pair); r_.append(r)
            r_ = pd.DataFrame(r_); r_["p_fdr10"] = fdr(r_.p)
            hh.append(r_.assign(part=part, contrast=f"{a}-{b}", measure=col, source="H7 per-subject files"))
    pd.concat(hh, ignore_index=True).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_T5_H7_same_participants.csv"),
                                            index=False)
    print("wrote EXPLORATORY_LOC_C2_T5*.csv", flush=True)


# --------------------------------------------------------------------------- figure
def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
    from matplotlib.patches import Rectangle
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    res = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C2_netfc.csv"))
    res = res[res.primary_spec]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2})
    cons = ["C1", "C2", "C4", "C5"]
    titles = {"C1": "C1 Intero − Extero", "C2": "C2 Intero − Self", "C4": "C4 Match − Extero",
              "C5": "C5 Match − Intero"}
    rows_ = [("total", "z", "Network FC\n(nuisance only)"), ("background", "z", "Background FC\n(task removed)"),
             ("background", "z_gs", "Background FC,\nglobal signal partialled")]
    fig, axes = plt.subplots(3, 4, figsize=(14.0, 10.2), facecolor=SURF)
    fig.subplots_adjust(top=0.93, hspace=0.32)
    lim = 1.0; norm = Normalize(-lim, lim)
    for i, (task, meas, ylab) in enumerate(rows_):
        for j, con in enumerate(cons):
            ax = axes[i, j]; ax.set_facecolor(SURF)
            r = res[(res.analysis == task) & (res.contrast == con) & (res.measure == meas)].set_index("pair")
            for (ia, na), (ib, nb) in itertools.combinations(enumerate(NETS), 2):
                row = r.loc[f"{na}-{nb}"]
                star = "**" if row.p_fdr10 < .05 else ("*" if row.p < .05 else "")
                ax.add_patch(Rectangle((ia - 0.5, ib - 0.5), 1, 1, facecolor=cmap(norm(row.dz)), edgecolor=SURF, lw=2))
                ax.text(ia, ib, f"{row.dz:+.2f}{star}", ha="center", va="center", fontsize=7.8, color=INK)
            ax.set_xticks(range(4)); ax.set_yticks(range(1, 5))
            ax.set_xticklabels(NETS[:4]); ax.set_yticklabels(NETS[1:])
            ax.set_xlim(-0.5, 3.5); ax.set_ylim(4.5, 0.5); ax.set_aspect("equal")
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.tick_params(length=0)
            ax.set_title(f"{titles[con]}  (N {int(r.N.iloc[0])})", fontsize=9.5, color=INK, loc="left")
            if j == 0:
                ax.set_ylabel(ylab, fontsize=10, color=INK)
    im = ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.015, pad=0.02)
    cb.set_label("dz (positive = higher coupling in the first condition)", color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle("Localizer network connectivity, condition differences (exploratory). Cell = dz of the paired Fisher-z difference; * p < .05 uncorrected; ** BH-FDR < .05 over 10 pairs",
                 x=0.02, ha="left", fontsize=10.5, color=INK)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "loc_C2_network_fc.png"), dpi=160, facecolor=SURF, bbox_inches="tight")
    print("wrote Figures/localizer/loc_C2_network_fc.png")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "subject":
        subject(sys.argv[2])
    elif cmd == "t5subject":
        t5_subject(sys.argv[2])
    else:
        {"verify": verify, "group": group, "figure": figure}[cmd]()
