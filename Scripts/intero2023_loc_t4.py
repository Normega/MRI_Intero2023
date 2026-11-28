#!/usr/bin/env python
"""
intero2023_loc_t4.py -- localizer stage C3 (EXPLORATORY; OSF p5ztc plan §8 T4, §10.1):
per-epoch expression of the localizer C1 (Intero - Extero) pattern in the breathing runs.
Everything this script writes is exploratory.

Choices: Analysis/localizer_stageC3_implementation_notes.md §1 (committed before any value).

  python Scripts/intero2023_loc_t4.py lss <sub>     # LSS betas for every paced and free epoch (SAL u SMN, 4 mm)
  python Scripts/intero2023_loc_t4.py verify <sub>  # design checks for one participant (no index computed)
  python Scripts/intero2023_loc_t4.py index         # per-epoch indices -> EXPLORATORY_LOC_C3_T4_epochs.csv
  python Scripts/intero2023_loc_t4.py reliability   # split-half / cross-run reliability of the index
  (mixed models: Scripts/intero2023_loc_t4_mixed.R)

LSS (Mumford et al. 2012, the "by condition" form): for each target epoch, the registered swc7p
Model A design of that run (q.model_a: paced_aware / _unaware / _wrongdir / _nochange, free,
rating, subrating; SPM + derivative; nuisance sheet; drift none) with the target epoch moved to
its own regressor `target` (+ derivative). AR(1) via intero2023_firstlevel.fit; the target's
canonical beta (PSC) is kept. Data: breathing BOLD smoothed 4 mm (plan §3: pattern analyses),
PSC per run, voxels of CAREN SAL u SMN (13,189 voxels).

Index (per epoch): with template w (SAL u SMN voxels) and epoch beta b,
  dot   = b . w / (||w|| sqrt(V))     primary; the projection of b on the unit template, per voxel
  corr  = Pearson r(b, w) over voxels  secondary; spatial shape only (mean removed)
Templates: own = participant's joint 4 mm localizer C1 effect; loso = OLS t of the joint 4 mm C1
effect over the localizer Intero set (48) without the participant (as the stage-B fallback);
run1 / run2 = the participant's per-run 4 mm C1 effect (reliability only).
Sign: C1 = Intero - Extero is mostly NEGATIVE (stage A), so a positive dot index means the epoch
is lower where Intero < Extero, i.e. it looks more like Intero than Extero, including the global
deactivation; a positive corr means the same spatial shape, regardless of the epoch's mean level.
Outputs: $OUT/loc_c3/t4/<sub>_lss.npz (betas; outside the repo); Analysis/results/EXPLORATORY_LOC_C3_T4_*.csv
"""
import os, sys, json, time, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g
import intero2023_loc_common as L
from intero2023_firstlevel import fit
q = g.q
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
OD = os.path.join(g.OUT, "loc_c3", "t4")
NETS = ["net_SAL", "net_SMN"]
EPOCH_CONDS = ["paced_aware", "paced_unaware", "paced_wrongdir", "paced_nochange", "free"]
OUTCOME = {"paced_aware": "aware", "paced_unaware": "unaware", "paced_nochange": "nochange",
           "paced_wrongdir": "wrongdir", "free": "free"}
FWHM = 4


def union_mask():
    gm = nib.load(g.GROUP_MASK)
    m = np.zeros(gm.shape, bool)
    for n in NETS:
        m |= nib.load(os.path.join(g.REPO, "ROIs", f"{n}.nii")).get_fdata() > 0
    assert (m <= (gm.get_fdata() > 0)).all()
    return nib.Nifti1Image(m.astype(np.uint8), gm.affine), m


def net_labels(m):
    """For each union voxel (C order), which network it belongs to."""
    lab = np.empty(m.sum(), object)
    for n in NETS:
        v = nib.load(os.path.join(g.REPO, "ROIs", f"{n}.nii")).get_fdata() > 0
        lab[v[m]] = n[4:]
    return lab


def t4_set():
    s, _, _ = L.participant_sets()
    bsubs = set(g.run_list().subject)
    return [x for x in s["intero"] if x in bsubs], s


# --------------------------------------------------------------------------- LSS
def lss_designs(ev, cf, n):
    """Yields (row index in ev, design) for every paced / free epoch of the run."""
    base = ev[["onset", "duration", "condition"]].reset_index(drop=True)
    for i in np.flatnonzero(base.condition.isin(EPOCH_CONDS).values):
        e = base.copy()
        e.loc[i, "condition"] = "target"
        yield i, q.model_a(e, cf, n)


def lss(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    t0 = time.time()
    img, _ = union_mask()
    masker = NiftiMasker(mask_img=img, smoothing_fwhm=FWHM, standardize=False).fit()
    man = g.run_list(); runs = man[man.subject == sub]
    B, Vv, rows = [], [], []
    for r in runs.itertuples():
        n = int(r.n_vol)
        Y = masker.transform(g.bold_path(sub, r.bids_entity)); assert Y.shape[0] == n
        Y, _ = mean_scaling(Y, axis=0)
        ev, cf, _ = g.read_run_physio(sub, r.bids_entity, n, [])
        ev = ev.reset_index(drop=True)
        for i, dm in lss_designs(ev, cf, n):
            out, dof = fit(Y, dm, {"target": {"target": 1.0}})
            e, v = out["target"]
            B.append(e); Vv.append(v)
            rec = ev.loc[i]
            rows.append(dict(subject=sub, entity=r.bids_entity, run=int(r.bids_entity[-1]),
                             epoch_index=rec.epoch_index, condition=rec.condition, onset=rec.onset,
                             duration=rec.duration, confidence=rec.confidence, response=rec.response,
                             change_direction=rec.change_direction, salience=rec.salience, dof=dof,
                             n_reg=dm.shape[1]))
    os.makedirs(OD, exist_ok=True)
    np.savez_compressed(os.path.join(OD, f"{sub}_lss.npz"), beta=np.array(B, np.float32),
                        variance=np.array(Vv, np.float32))
    pd.DataFrame(rows).to_csv(os.path.join(OD, f"{sub}_lss_epochs.csv"), index=False)
    print(f"{sub}: {len(rows)} epochs in {time.time() - t0:.0f} s", flush=True)


def verify(sub):
    """Design checks (no index): LSS designs full rank; target column is the epoch's own boxcar;
    the LSS designs with the target merged back equal Model A."""
    man = g.run_list(); runs = man[man.subject == sub]
    res = dict(subject=sub, runs=[])
    for r in runs.itertuples():
        n = int(r.n_vol)
        ev, cf, _ = g.read_run_physio(sub, r.bids_entity, n, [])
        ev = ev.reset_index(drop=True)
        dmA = q.model_a(ev, cf, n)
        ranks, recon = [], []
        for i, dm in lss_designs(ev, cf, n):
            X = dm.values.astype(float); sd = X.std(0); sd[sd < 1e-12] = 1
            ranks.append(int(np.linalg.matrix_rank(X / sd)) == X.shape[1])
            c = ev.condition[i]
            merged = dm["target"] + (dm[c] if c in dm.columns else 0)
            recon.append(float(np.abs(merged.values - dmA[c].values).max()))
        res["runs"].append(dict(entity=r.bids_entity, n_epochs=len(ranks), all_full_rank=bool(all(ranks)),
                                max_abs_merge_diff=max(recon)))
    path = os.path.join(g.REPO, "Analysis", "localizer_stageC3_t4_design_check.json")
    json.dump(res, open(path, "w"), indent=1)
    print(json.dumps(res, indent=1))


# --------------------------------------------------------------------------- index
def loso_t(stack):
    n = stack.shape[0]; S1 = stack.sum(0); S2 = (stack ** 2).sum(0)
    out = np.empty_like(stack)
    for i in range(n):
        m = (S1 - stack[i]) / (n - 1)
        var = ((S2 - stack[i] ** 2) - (n - 1) * m ** 2) / (n - 2)
        out[i] = m / np.sqrt(np.maximum(var, 1e-20) / (n - 1))
    return out


def dot_index(B, w):
    return B @ w / (np.linalg.norm(w) * np.sqrt(len(w)))


def corr_index(B, w):
    Bc = B - B.mean(1, keepdims=True); wc = w - w.mean()
    return Bc @ wc / (np.linalg.norm(Bc, axis=1) * np.linalg.norm(wc))


def index():
    subs, s = t4_set()
    _, m = union_mask(); lab = net_labels(m)
    load = lambda p: nib.load(p).get_fdata(dtype=np.float32)[m].astype(np.float64)
    it = s["intero"]
    own = {x: load(L.fl_name(x, "joint", 4, "base", "C1", "effect")) for x in it}
    T = loso_t(np.vstack([own[x] for x in it]))
    loso = {x: T[i] for i, x in enumerate(it)}
    rows = []
    for x in subs:
        z = np.load(os.path.join(OD, f"{x}_lss.npz")); B = z["beta"].astype(np.float64)
        E = pd.read_csv(os.path.join(OD, f"{x}_lss_epochs.csv"))
        assert len(E) == B.shape[0]
        tmpl = {"own": own[x], "loso": loso[x]}
        for k in ("run1", "run2"):
            p = L.fl_name(x, k, 4, "base", "C1", "effect")
            if os.path.exists(p):
                tmpl[k] = load(p)
        for name, w in tmpl.items():
            E[f"dot_{name}"] = dot_index(B, w)
            E[f"corr_{name}"] = corr_index(B, w)
            if name in ("own", "loso"):
                for net in ("SAL", "SMN"):
                    sel = lab == net
                    E[f"dot_{name}_{net}"] = dot_index(B[:, sel], w[sel])
                    E[f"corr_{name}_{net}"] = corr_index(B[:, sel], w[sel])
        E["mean_beta"] = B.mean(1)
        E["template_mean_own"] = own[x].mean()
        E["template_rms_own"] = np.sqrt((own[x] ** 2).mean())
        E["outcome"] = E.condition.map(OUTCOME)
        E["epoch_type"] = np.where(E.condition == "free", "free", "paced")
        rows.append(E)
    D = pd.concat(rows, ignore_index=True)
    D["confidence"] = pd.to_numeric(D.confidence, errors="coerce")
    D.drop(columns=["onset", "duration"]).round(6).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_epochs.csv"), index=False)
    print(f"index: {D.subject.nunique()} participants, {len(D)} epochs "
          f"({(D.epoch_type == 'paced').sum()} paced)")


# --------------------------------------------------------------------------- reliability
def sb(r):
    return 2 * r / (1 + r)


def reliability():
    D = pd.read_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_epochs.csv"))
    P = D[D.epoch_type == "paced"].copy()
    out = []
    for kind in ("dot", "corr"):
        # (1) template split-half: the index from the run-1 vs the run-2 localizer C1 map,
        #     across paced epochs, within participant (participant means removed)
        d = P.dropna(subset=[f"{kind}_run1", f"{kind}_run2"])
        a = d[f"{kind}_run1"] - d.groupby("subject")[f"{kind}_run1"].transform("mean")
        b = d[f"{kind}_run2"] - d.groupby("subject")[f"{kind}_run2"].transform("mean")
        r = np.corrcoef(a, b)[0, 1]
        per = d.groupby("subject").apply(lambda u: np.corrcoef(u[f"{kind}_run1"], u[f"{kind}_run2"])[0, 1])
        out.append(dict(index=kind, measure="template split-half (loc run 1 vs run 2 C1 map), within-participant epoch r",
                        n_participants=d.subject.nunique(), n_epochs=len(d), r=r, spearman_brown=sb(r),
                        per_participant_r_median=float(per.median())))
        # between-participant part: participant mean index from each template
        mm = d.groupby("subject")[[f"{kind}_run1", f"{kind}_run2"]].mean()
        r = np.corrcoef(mm.iloc[:, 0], mm.iloc[:, 1])[0, 1]
        out.append(dict(index=kind, measure="template split-half, participant mean paced index (between-participant r)",
                        n_participants=len(mm), n_epochs=len(d), r=r, spearman_brown=sb(r)))
        # (2) breathing-run split-half with the own (two-run) template
        for tm in ("own", "loso"):
            col = f"{kind}_{tm}"
            mm = P.pivot_table(index="subject", columns="run", values=col, aggfunc="mean").dropna()
            r = np.corrcoef(mm[1], mm[2])[0, 1]
            out.append(dict(index=kind, template=tm, measure="breathing run 1 vs run 2, participant mean paced index",
                            n_participants=len(mm), r=r, spearman_brown=sb(r)))
            # outcome effect across breathing runs: aware - unaware per run
            Q = P[P.outcome.isin(["aware", "unaware"])].pivot_table(index=["subject", "run"], columns="outcome",
                                                                     values=col, aggfunc="mean").dropna()
            Q = (Q.aware - Q.unaware).unstack().dropna()
            if len(Q) >= 10:
                r = np.corrcoef(Q[1], Q[2])[0, 1]
                out.append(dict(index=kind, template=tm, measure="breathing run 1 vs run 2, participant aware - unaware",
                                n_participants=len(Q), r=r, spearman_brown=sb(r)))
            else:
                out.append(dict(index=kind, template=tm, measure="breathing run 1 vs run 2, participant aware - unaware",
                                n_participants=len(Q), note="fewer than 10 participants with both cells in both runs"))
    pd.DataFrame(out).round(4).to_csv(os.path.join(RD, "EXPLORATORY_LOC_C3_T4_reliability.csv"), index=False)
    print(pd.DataFrame(out).round(3).to_string())


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd in ("lss", "verify"):
        {"lss": lss, "verify": verify}[cmd](sys.argv[2])
    else:
        {"index": index, "reliability": reliability}[cmd]()
