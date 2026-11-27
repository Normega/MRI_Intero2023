#!/usr/bin/env python
"""
intero2023_loc_froi.py -- localizer stage B (OSF p5ztc; plan §7-§9): participant-specific
ROIs (fROIs), the reliability gate L4 with the LOSO group fallback, and the transfer tests
T1 (registered), T3b (registered), T3 (secondary), T2 and the C2 / C3 fROI rules
(exploratory), with the registered stage-B sensitivity analyses.
Choices: Analysis/localizer_stageB_implementation_notes.md (committed 03ef477, before any
stage-B statistic). Contrast mapping (swc7p numbering for breathing files):
  breathing C2joint = aware - unaware (H3) -> T1;  breathing C1 = paced - free (H1/H6) -> T3, T3b;
  breathing C3 = confidence (H2), E3changejoint = change - No Change (E3) -> T2.

  python Scripts/intero2023_loc_froi.py define     # fROIs, L4 gate, definition tables (localizer only)
  python Scripts/intero2023_loc_froi.py transfer   # extraction and T1 / T3b / T3 / T2 / exploratory / sensitivity

Reuses: intero2023_group (sh, reml_meta, subject_sets), intero2023_explore_awareness
(jzs_bf10, bf_plus0), intero2023_loc_group (jzs_reg_bf10), intero2023_loc_common (sets, names).
fROI masks: $OUT/froi/ (outside the repo). Tables: Analysis/results/LOC_*.csv.
"""
import os, sys, json, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
import intero2023_group as G
from intero2023_explore_awareness import jzs_bf10, bf_plus0
from intero2023_loc_group import jzs_reg_bf10
g = L.g
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
FD = os.path.join(g.OUT, "froi")
PARCELS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere"]
T1_PARCELS = PARCELS[:3]
SPHERE = "ieat_ACCsphere"
FRACS = (0.10, 0.05, 0.20)
MIN_VOX = 50
R_BF = 0.707
N_PERM, PERM_SEED, BOOT_SEED, N_BOOT = 10000, 20261126, 20261127, 10000
ALPHA = 0.05
# rule -> (localizer ranking map, direction). C1 is primary; C2, C3 exploratory; C5neg secondary (T3);
# C1bottom is the registered specificity control (Extero > Intero).
RULES = {"C1": ("C1", +1), "C2": ("C2", +1), "C3": ("C3", +1), "C5neg": ("C5", -1), "C1bottom": ("C1", -1)}
GATED = ("C1", "C2", "C3")


# --------------------------------------------------------------------------- helpers
def masks():
    gm = nib.load(g.GROUP_MASK).get_fdata() > 0
    out = {}
    for p in PARCELS:
        m = nib.load(G.ROI(p)).get_fdata() > 0
        assert (m <= gm).all(), p
        out[p] = m
    return out


PM = masks()


def load_parcels(path):
    """Raw voxel values (no smoothing) per parcel, in mask-voxel (C) order."""
    d = nib.load(path).get_fdata(dtype=np.float32)
    return {p: d[m].astype(np.float64) for p, m in PM.items()}


def pload(paths, jobs=8):
    with ThreadPoolExecutor(jobs) as ex:
        return list(ex.map(load_parcels, paths))


def froi_size(nvox, frac):
    return int(max(np.ceil(round(frac * nvox, 9)), MIN_VOX))


def select(vec, n, sign=+1):
    order = np.argsort(-sign * vec, kind="stable")
    sel = np.zeros(len(vec), bool); sel[order[:n]] = True
    return sel


def dice(a, b):
    return 2 * (a & b).sum() / (a.sum() + b.sum())


def fz(r):
    return np.arctanh(np.clip(r, -0.999999, 0.999999))


def sets():
    s, maia, u = L.participant_sets()
    bs, bmaia, man = G.subject_sets()
    h2 = json.load(open(os.path.join(RD, "H2_var-base_set-all_subjects.json")))["subjects"]
    e3 = json.load(open(os.path.join(RD, "EXPLORATORY_E3a_change_vs_nochange_subjects.json")))["subjects"]
    fr = s["froi"]
    S = dict(froi=fr, intero=s["intero"], v1=s["v1"],
             T1=[x for x in fr if x in bs["H3"]], T3=[x for x in fr if x in bs["all"]],
             T3b=[x for x in fr if x in bs["H6"]], T2conf=[x for x in fr if x in h2],
             T2chg=[x for x in fr if x in e3])
    S["T1nov1"] = [x for x in S["T1"] if x not in s["v1"]]
    S["T3bnov1"] = [x for x in S["T3b"] if x not in s["v1"]]
    return S, bmaia


# --------------------------------------------------------------------------- define + gate
def loso_t(stack):
    """Leave-one-out OLS one-sample t for each row of stack (n x v)."""
    n = stack.shape[0]; S1 = stack.sum(0); S2 = (stack ** 2).sum(0)
    out = np.empty_like(stack)
    for i in range(n):
        m = (S1 - stack[i]) / (n - 1)
        var = ((S2 - stack[i] ** 2) - (n - 1) * m ** 2) / (n - 2)
        out[i] = m / np.sqrt(np.maximum(var, 1e-20) / (n - 1))
    return out


def define():
    S, _ = sets()
    fr = S["froi"]; os.makedirs(FD, exist_ok=True)
    print(f"fROI-eligible N = {len(fr)}", flush=True)
    joint = {c: pload([L.fl_name(x, "joint", 4, "base", c, "t") for x in fr]) for c in ("C1", "C2", "C3", "C5")}
    jeff = {c: pload([L.fl_name(x, "joint", 4, "base", c, "effect") for x in fr]) for c in ("C1",)}
    run_t = {(k, c): pload([L.fl_name(x, k, 4, "base", c, "t") for x in fr]) for k in ("run1", "run2") for c in ("C1", "C2", "C3")}
    run_e = {(k, c): pload([L.fl_name(x, k, 4, "base", c, "effect") for x in fr]) for k in ("run1", "run2") for c in ("C1", "C2")}
    # LOSO group maps over the Intero set (48), target left out
    it = S["intero"]
    ieff = {c: pload([L.fl_name(x, "joint", 4, "base", c, "effect") for x in it]) for c in ("C1", "C2")}
    loso = {}
    for c in ("C1", "C2"):
        for p in PARCELS:
            T = loso_t(np.vstack([d[p] for d in ieff[c]]))
            for i, x in enumerate(it):
                loso[(c, p, x)] = T[i]
    froi, drows = {}, []
    for i, x in enumerate(fr):
        for p in PARCELS:
            nv = PM[p].sum()
            c1t, c1e = joint["C1"][i][p], jeff["C1"][i][p]
            maps = {r: (joint[c][i][p], sg) for r, (c, sg) in RULES.items()}
            maps["C1loso"] = (loso[("C1", p, x)], +1)
            maps["C2loso"] = (loso[("C2", p, x)], +1)
            maps["C3loso"] = (np.minimum(loso[("C1", p, x)], loso[("C2", p, x)]), +1)
            for rule, (mp, sg) in maps.items():
                fracs = FRACS if rule in ("C1", "C1loso") else (0.10,)
                for f in fracs:
                    sel = select(mp, froi_size(nv, f), sg)
                    froi[(x, rule, f, p)] = sel
                    drows.append(dict(subject=x, rule=rule, frac=f, parcel=p, n_vox=int(sel.sum()), parcel_vox=int(nv),
                                      froi_C1_t=float(c1t[sel].mean()), froi_C1_effect=float(c1e[sel].mean()),
                                      comp_C1_t=float(c1t[~sel].mean()), comp_C1_effect=float(c1e[~sel].mean()),
                                      parcel_C1_t=float(c1t.mean()), parcel_C1_effect=float(c1e.mean()),
                                      froi_rankmap_mean=float(mp[sel].mean())))
    # save fROI masks (derivatives only): one 4-D image per subject x rule (10%), a volume per parcel
    ref = nib.load(g.GROUP_MASK)
    for x in fr:
        od = os.path.join(FD, x); os.makedirs(od, exist_ok=True)
        for rule in list(RULES) + ["C1loso", "C2loso", "C3loso"]:
            vol = np.zeros(ref.shape + (len(PARCELS),), np.uint8)
            for k, p in enumerate(PARCELS):
                v = np.zeros(ref.shape, bool); v[PM[p]] = froi[(x, rule, 0.10, p)]; vol[..., k] = v
            nib.Nifti1Image(vol, ref.affine).to_filename(os.path.join(od, f"{x}_froi-{rule}_frac-10_parcels-5.nii.gz"))
    np.savez_compressed(os.path.join(FD, "froi_masks.npz"),
                        **{f"{x}|{r}|{f:.2f}|{p}": v for (x, r, f, p), v in froi.items()})
    json.dump(dict(parcels=PARCELS, subjects=fr), open(os.path.join(FD, "froi_index.json"), "w"), indent=1)
    D = pd.DataFrame(drows)
    D.round(6).to_csv(os.path.join(RD, "LOC_froi_definition.csv"), index=False)
    summ = []
    for (rule, f, p), d in D.groupby(["rule", "frac", "parcel"]):
        summ.append(dict(rule=rule, frac=f, parcel=p, N=len(d), n_vox=int(d.n_vox.iloc[0]),
                         froi_C1_t_median=d.froi_C1_t.median(), froi_C1_t_q25=d.froi_C1_t.quantile(.25),
                         froi_C1_t_q75=d.froi_C1_t.quantile(.75), froi_C1_effect_mean=d.froi_C1_effect.mean(),
                         froi_C1_effect_median=d.froi_C1_effect.median(),
                         share_froi_C1_effect_gt0=float((d.froi_C1_effect > 0).mean()),
                         share_froi_C1_t_gt0=float((d.froi_C1_t > 0).mean()),
                         parcel_C1_effect_mean=d.parcel_C1_effect.mean(), parcel_C1_t_median=d.parcel_C1_t.median(),
                         comp_C1_effect_mean=d.comp_C1_effect.mean()))
    pd.DataFrame(summ).round(6).to_csv(os.path.join(RD, "LOC_froi_definition_summary.csv"), index=False)
    gate(fr, run_t, run_e)
    print("define: done", flush=True)


def gate(fr, run_t, run_e):
    n = len(fr); rows, srows = [], []
    for rule in GATED:
        rank_c = rule
        eff_cs = ("C1", "C2") if rule == "C3" else (rule,)
        for p in PARCELS:
            nv = PM[p].sum(); k = froi_size(nv, 0.10)
            T1 = np.vstack([run_t[("run1", rank_c)][i][p] for i in range(n)])
            T2 = np.vstack([run_t[("run2", rank_c)][i][p] for i in range(n)])
            cross = {c: [] for c in eff_cs}; agree = {c: [] for c in eff_cs}; dc = []
            for i in range(n):
                a, b = select(T1[i], k), select(T2[i], k)
                dc.append(dice(a, b))
                for c in eff_cs:
                    e1, e2 = run_e[("run1", c)][i][p], run_e[("run2", c)][i][p]
                    d12, d21 = e2[a].mean(), e1[b].mean()       # defined in run j, estimated in run k
                    cross[c].append((d12 + d21) / 2); agree[c].append(np.sign(d12) == np.sign(d21))
            Z1 = (T1 - T1.mean(1, keepdims=True)) / T1.std(1, keepdims=True)
            Z2 = (T2 - T2.mean(1, keepdims=True)) / T2.std(1, keepdims=True)
            C = Z1 @ Z2.T / nv
            w = np.diag(C); zw = fz(w)
            zb = np.array([fz(np.delete(C[i], i)).mean() for i in range(n)])
            rb = np.tanh(zb)
            pt = stats.ttest_rel(zw, zb, alternative="greater")
            pw = stats.wilcoxon(zw - zb, alternative="greater")
            rec = dict(rule=rule, parcel=p, N=n, froi_vox=k, parcel_vox=int(nv))
            c1_ok = True
            for c in eff_cs:
                y = np.array(cross[c]); tt = stats.ttest_1samp(y, 0, alternative="greater")
                ok = bool(tt.pvalue < ALPHA and (y > 0).mean() >= 0.70)
                c1_ok &= ok
                rec.update({f"cross_{c}_mean": y.mean(), f"cross_{c}_dz": y.mean() / y.std(ddof=1),
                            f"cross_{c}_t": tt.statistic, f"cross_{c}_p_one": tt.pvalue,
                            f"cross_{c}_share_gt0": (y > 0).mean(),
                            f"cross_{c}_share_dir_agree": float(np.mean(agree[c])), f"crit1_{c}_pass": ok})
            dc = np.array(dc)
            rec.update(crit1_pass=c1_ok, dice_median=float(np.median(dc)), dice_q25=float(np.quantile(dc, .25)),
                       dice_q75=float(np.quantile(dc, .75)), dice_chance=k / nv, crit2_pass=bool(np.median(dc) >= 0.20),
                       r_within_median=float(np.median(w)), r_between_median=float(np.median(rb)),
                       z_within_minus_between_mean=float((zw - zb).mean()),
                       z_diff_dz=float((zw - zb).mean() / (zw - zb).std(ddof=1)),
                       crit3_t=float(pt.statistic), crit3_p_one=float(pt.pvalue), crit3_wilcoxon_p_one=float(pw.pvalue),
                       crit3_pass=bool(pt.pvalue < ALPHA))
            rec["L4_met"] = bool(rec["crit1_pass"] and rec["crit2_pass"] and rec["crit3_pass"])
            rec["fallback_used"] = not rec["crit3_pass"]
            rec["primary_froi"] = f"{rule}loso" if rec["fallback_used"] else rule
            rows.append(rec)
            for i, x in enumerate(fr):
                srows.append(dict(rule=rule, parcel=p, subject=x, dice=dc[i], r_within=w[i], r_between_mean=rb[i],
                                  **{f"cross_{c}": cross[c][i] for c in eff_cs}))
    pd.DataFrame(rows).to_csv(os.path.join(RD, "LOC_L4_gate.csv"), index=False)
    pd.DataFrame(srows).round(6).to_csv(os.path.join(RD, "LOC_L4_gate_subjects.csv"), index=False)


# --------------------------------------------------------------------------- extraction
def breath_path(x, sm, con, kind):
    d = os.path.join(g.OUT, "firstlevel", x)
    return os.path.join(d, f"{x}_var-base_con-{con}_{kind}.nii.gz") if sm == 6 else \
        os.path.join(d, f"{x}_sm-{sm}_var-base_con-{con}_{kind}.nii.gz")


def breath_dof(x, con):
    js = json.load(open(os.path.join(g.OUT, "firstlevel", x, f"{x}_breath_sm4_sm0.json")))["fits"]["sm4"]
    if con in ("C2joint", "E3changejoint"):
        return js["joint"]["dof"]
    key = "dof_A" if con == "C1" else "dof_B"
    return sum(r[key] for r in js["runs"] if con in r["contrasts"])


def load_froi():
    z = np.load(os.path.join(FD, "froi_masks.npz"))
    out = {}
    for k in z.files:
        x, r, f, p = k.split("|"); out[(x, r, float(f), p)] = z[k]
    return out


EXTRACT_PLAN = {   # smoothing -> contrasts; rules (with fractions) are taken from the saved fROIs
    4: ("C2joint", "C1", "C3", "E3changejoint"),
    0: ("C2joint", "C1"),
    6: ("C2joint", "C1"),
}


def extract(S):
    F = load_froi()
    subs = S["T3"]                      # every fROI-eligible participant with breathing data (41)
    rows = []
    for sm, cons in EXTRACT_PLAN.items():
        for con in cons:
            have = [x for x in subs if os.path.exists(breath_path(x, sm, con, "effect"))]
            E = dict(zip(have, pload([breath_path(x, sm, con, "effect") for x in have])))
            V = dict(zip(have, pload([breath_path(x, sm, con, "variance") for x in have])))
            for (x, r, f, p), sel in F.items():
                if x not in E or (sm != 4 and not (r in ("C1", "C1loso") or (f == 0.10 and r in ("C2", "C3", "C2loso", "C3loso")))):
                    continue
                e, v = E[x][p], V[x][p]
                rows.append(dict(subject=x, sm=sm, contrast=con, rule=r, frac=f, parcel=p,
                                 froi=e[sel].mean(), froi_var=v[sel].mean(), comp=e[~sel].mean(),
                                 comp_var=v[~sel].mean(), parcel_mean=e.mean(), parcel_var=v.mean()))
            print(f"extract sm{sm} {con}: {len(have)} participants", flush=True)
    X = pd.DataFrame(rows)
    X.to_csv(os.path.join(RD, "LOC_froi_extract.csv"), index=False, float_format="%.6g")
    return X


# --------------------------------------------------------------------------- statistics
def flame_batch(tag, items, subs, con):
    """FLAME 1 (flameo) on index values; one 'voxel' per item. items: list of (key, y, v)."""
    od = os.path.join(g.OUT, "group", f"LOC_stageB_flame_{tag}"); os.makedirs(od, exist_ok=True)
    n, m = len(subs), len(items)
    cope = np.zeros((m, 1, 1, n), np.float32); var = np.zeros_like(cope)
    for j, (_, y, v) in enumerate(items):
        cope[j, 0, 0] = y; var[j, 0, 0] = v
    A = np.eye(4)
    dof = np.array([breath_dof(x, con) for x in subs], np.float32)
    nib.Nifti1Image(cope, A).to_filename(f"{od}/cope.nii.gz")
    nib.Nifti1Image(var, A).to_filename(f"{od}/varcope.nii.gz")
    nib.Nifti1Image(np.broadcast_to(dof, (m, 1, 1, n)).astype(np.float32), A).to_filename(f"{od}/dof.nii.gz")
    nib.Nifti1Image(np.ones((m, 1, 1), np.uint8), A).to_filename(f"{od}/mask.nii.gz")
    open(f"{od}/design.mat", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    open(f"{od}/design.con", "w").write("/NumWaves 1\n/NumContrasts 1\n/Matrix\n1\n")
    open(f"{od}/design.grp", "w").write(f"/NumWaves 1\n/NumPoints {n}\n/Matrix\n" + "1\n" * n)
    G.sh(f"rm -rf {od}/stats; flameo --cope={od}/cope.nii.gz --vc={od}/varcope.nii.gz --dvc={od}/dof.nii.gz "
         f"--mask={od}/mask.nii.gz --ld={od}/stats --dm={od}/design.mat --tc={od}/design.con "
         f"--cs={od}/design.grp --runmode=flame1")
    rd = lambda f: nib.load(f"{od}/stats/{f}.nii.gz").get_fdata().ravel()
    z, t, tdof, re = rd("zstat1"), rd("tstat1"), rd("tdof_t1"), rd("mean_random_effects_var1")
    json.dump(dict(tag=tag, subjects=subs, dofs=dof.tolist(), keys=[str(k) for k, _, _ in items]),
              open(f"{od}/meta.json", "w"), indent=1)
    return {k: dict(flame_z=float(z[j]), flame_t=float(t[j]), flame_tdof=float(tdof[j]), revar=float(re[j]))
            for j, (k, _, _) in enumerate(items)}


def boot_idx(n):
    return np.random.default_rng(BOOT_SEED).integers(0, n, size=(N_BOOT, n))


def dz(y):
    return y.mean() / y.std(ddof=1)


def paired(yf, yc):
    d = yf - yc; B = boot_idx(len(d))
    bd = np.array([dz(d[b]) for b in B])
    bdd = np.array([dz(yf[b]) - dz(yc[b]) for b in B])
    tt = stats.ttest_rel(yf, yc, alternative="greater")
    return dict(diff_mean=d.mean(), diff_dz=dz(d), diff_dz_ci_lo=np.quantile(bd, .025), diff_dz_ci_hi=np.quantile(bd, .975),
                dzdiff=dz(yf) - dz(yc), dzdiff_ci_lo=np.quantile(bdd, .025), dzdiff_ci_hi=np.quantile(bdd, .975),
                paired_t=tt.statistic, paired_p_one=tt.pvalue)


def signflip_p(y, sign):
    """One-sided sign-flip permutation p of the one-sample t in direction `sign`."""
    rng = np.random.default_rng(PERM_SEED); n = len(y)
    t0 = sign * stats.ttest_1samp(y, 0).statistic
    fl = rng.choice([-1.0, 1.0], size=(N_PERM, n)) * y
    tp = sign * fl.mean(1) / (fl.std(1, ddof=1) / np.sqrt(n))
    return (1 + (tp >= t0).sum()) / (1 + N_PERM)


def maia_test(y, m):
    """OLS y ~ 1 + demeaned MAIA: slope t, one-sided (positive) permutation p (10,000, seed 20261126)."""
    x = m - m.mean(); n = len(y)
    def slope_t(xx):
        lr = stats.linregress(xx, y); return lr.slope, lr.slope / lr.stderr
    b, t0 = slope_t(x)
    rng = np.random.default_rng(PERM_SEED)
    tp = np.array([slope_t(rng.permutation(x))[1] for _ in range(N_PERM)])
    r = stats.pearsonr(x, y)[0]
    return dict(slope=b, slope_t=t0, p_param_one=stats.t.sf(t0, n - 2), p_perm_one=(1 + (tp >= t0).sum()) / (1 + N_PERM),
                r=r, bf10=jzs_reg_bf10(r ** 2, n))


def bf_cols(y, direction=+1):
    t = stats.ttest_1samp(y, 0).statistic
    b10 = jzs_bf10(float(t), len(y), R_BF)
    bdir = bf_plus0(direction * y, R_BF)[1]
    return dict(BF10=b10, BF01=1 / b10, BFdir=bdir, BF0dir=1 / bdir)


def get(X, subs, sm, con, rule, frac, parcel):
    d = X[(X.sm == sm) & (X.contrast == con) & (X.rule == rule) & (np.isclose(X.frac, frac)) & (X.parcel == parcel)]
    d = d.set_index("subject").reindex(subs)
    assert d.froi.notna().all(), (sm, con, rule, frac, parcel)
    return d


# --------------------------------------------------------------------------- transfer
def transfer():
    S, maia = sets()
    gt = pd.read_csv(os.path.join(RD, "LOC_L4_gate.csv"))
    prim = {(r.rule, r.parcel): r.primary_froi for r in gt.itertuples()}
    X = extract(S)
    e5 = pd.read_csv(os.path.join(RD, "EXPLORATORY_E5_bayes_factors.csv"))
    e5 = e5[(e5.target == "H3") & np.isclose(e5.prior_r, 0.707)].set_index("roi")
    h3 = pd.read_csv(os.path.join(RD, "H3_var-base_set-all_roi.csv")).set_index("roi")

    # ---------- T1 family (primary, sensitivity, exploratory rules): FLAME on indices
    specs = []   # (spec name, set name, sm, rule-of-record, frac, is_primary_spec)
    for rule in GATED:
        for alt in ("individual", "loso"):
            specs.append((f"{rule}_{alt}", "T1", 4, rule, alt, 0.10))
    for sm in (0, 6):
        specs.append((f"C1_sm{sm}", "T1", sm, "C1", "gate", 0.10))
    for f in (0.05, 0.20):
        specs.append((f"C1_frac{int(f * 100):02d}", "T1", 4, "C1", "gate", f))
    specs.append(("C1_nov1", "T1nov1", 4, "C1", "gate", 0.10))
    t1rows, wrows = [], []
    for setname in ("T1", "T1nov1"):
        subs = S[setname]; items, meta = [], []
        for name, sn, sm, rule, alt, f in specs:
            if sn != setname:
                continue
            for p in T1_PARCELS:
                r_used = prim[(rule, p)] if alt == "gate" else (rule if alt == "individual" else f"{rule}loso")
                d = get(X, subs, sm, "C2joint", r_used, f, p)
                for reg, yc, vc in (("froi", "froi", "froi_var"), ("complement", "comp", "comp_var"),
                                    ("parcel", "parcel_mean", "parcel_var")):
                    key = (name, p, reg)
                    items.append((key, d[yc].values, d[vc].values))
                meta.append((name, sm, rule, alt, r_used, f, p, d))
        fl = flame_batch(setname, items, subs, "C2joint")
        for name, sm, rule, alt, r_used, f, p, d in meta:
            gate_row = gt[(gt.rule == rule) & (gt.parcel == p)].iloc[0]
            is_primary_froi = r_used == prim[(rule, p)]
            for reg, yc, vc in (("froi", "froi", "froi_var"), ("complement", "comp", "comp_var"),
                                ("parcel", "parcel_mean", "parcel_var")):
                y, v = d[yc].values, d[vc].values
                F = fl[(name, p, reg)]
                tt = stats.ttest_1samp(y, 0)
                w = 1 / (v + F["revar"]); share = w / w.sum()
                meta_ = G.reml_meta(y, v)
                rec = dict(spec=name, set=setname, N=len(y), sm=sm, rule=rule, froi_used=r_used,
                           froi_is_gate_primary=is_primary_froi, L4_met=bool(gate_row.L4_met), frac=f, parcel=p, region=reg,
                           mean=y.mean(), dz=dz(y), t=tt.statistic, p_one_t=stats.t.sf(tt.statistic, len(y) - 1),
                           flame_z=F["flame_z"], flame_t=F["flame_t"], flame_tdof=F["flame_tdof"],
                           flame_p_one=float(stats.norm.sf(F["flame_z"])),
                           flame_p_one_bonf3=float(min(1.0, 3 * stats.norm.sf(F["flame_z"]))),
                           reml_mu=meta_["mu"], reml_z=meta_["z"], max_weight_share=float(share.max()),
                           **bf_cols(y, +1))
                if reg == "froi":
                    rec.update(paired(y, d["comp"].values))
                    rec["T1_froi_component"] = bool(rec["flame_p_one_bonf3"] < ALPHA)
                    rec["T1_paired_component"] = bool(rec["diff_dz_ci_lo"] > 0)
                    rec["T1_supported"] = rec["T1_froi_component"] and rec["T1_paired_component"]
                    if setname == "T1" and sm == 4 and f == 0.10 and rule == "C1" and is_primary_froi:
                        for s_, sh_ in zip(d.index, share):
                            wrows.append(dict(parcel=p, subject=s_, weight_share=sh_))
                if reg == "parcel":
                    rec.update(E5_N=int(e5.loc[p, "N"]), E5_dz=e5.loc[p, "dz"], E5_BF01=e5.loc[p, "BF01"],
                               E5_BF0plus=e5.loc[p, "BF0plus"], H3_N=int(h3.loc[p, "N"]), H3_peak_z=h3.loc[p, "peak_z"],
                               H3_p_fwe_bonf=h3.loc[p, "p_fwe_bonf"], H3_roi_dz=h3.loc[p, "roi_dz_unweighted"])
                t1rows.append(rec)
    T1 = pd.DataFrame(t1rows)
    T1["primary_spec"] = (T1.spec.isin(["C1_individual", "C1_loso"])) & T1.froi_is_gate_primary & (T1.set == "T1")
    T1.to_csv(os.path.join(RD, "LOC_T1.csv"), index=False)
    pd.DataFrame(wrows).to_csv(os.path.join(RD, "LOC_T1_weights.csv"), index=False)

    # ---------- T3b family: MAIA x paced - free in the sphere fROI
    t3b = []
    for setname in ("T3b", "T3bnov1"):
        subs = S[setname]; m = maia.loc[subs, "maia_ieat_total"].astype(float).values
        cases = []
        if setname == "T3b":
            for rule in GATED:
                for alt in ("individual", "loso"):
                    cases.append((f"{rule}_{alt}", 4, rule, rule if alt == "individual" else f"{rule}loso", 0.10))
            for sm in (0, 6):
                cases.append((f"C1_sm{sm}", sm, "C1", prim[("C1", SPHERE)], 0.10))
            for f in (0.05, 0.20):
                cases.append((f"C1_frac{int(f * 100):02d}", 4, "C1", prim[("C1", SPHERE)], f))
            cases.append(("C5neg_T3", 4, "C5neg", "C5neg", 0.10))
        else:
            cases.append(("C1_nov1", 4, "C1", prim[("C1", SPHERE)], 0.10))
        for name, sm, rule, r_used, f in cases:
            d = get(X, subs, sm, "C1", r_used, f, SPHERE)
            res = {reg: maia_test(d[c].values, m) for reg, c in (("froi", "froi"), ("complement", "comp"), ("parcel", "parcel_mean"))}
            dd = d["froi"].values - d["comp"].values
            B = boot_idx(len(dd)); x = m - m.mean()
            bs = np.array([stats.linregress(x[b], dd[b]).slope for b in B])
            for reg, rr in res.items():
                rec = dict(spec=name, set=setname, N=len(subs), sm=sm, rule=rule, froi_used=r_used, frac=f, region=reg, **rr)
                if reg == "froi":
                    lr = stats.linregress(x, dd)
                    rec.update(diff_slope=lr.slope, diff_slope_t=lr.slope / lr.stderr,
                               diff_slope_ci_lo=np.quantile(bs, .025), diff_slope_ci_hi=np.quantile(bs, .975),
                               T3b_supported=bool(rr["p_perm_one"] < ALPHA),
                               froi_is_gate_primary=(rule in GATED and r_used == prim[(rule, SPHERE)]))
                t3b.append(rec)
    T3B = pd.DataFrame(t3b); T3B.to_csv(os.path.join(RD, "LOC_T3b.csv"), index=False)

    # ---------- T3 (secondary): paced - free < 0 in the C5neg fROIs, x5
    subs = S["T3"]; t3 = []
    for p in PARCELS:
        d = get(X, subs, 4, "C1", "C5neg", 0.10, p)
        for reg, c in (("froi", "froi"), ("complement", "comp"), ("parcel", "parcel_mean")):
            y = d[c].values; tt = stats.ttest_1samp(y, 0)
            pp = signflip_p(y, -1)
            rec = dict(test="T3", N=len(y), parcel=p, region=reg, mean=y.mean(), dz=dz(y), t=tt.statistic,
                       p_perm_one_neg=pp, p_perm_bonf5=min(1.0, 5 * pp), **bf_cols(y, -1))
            if reg == "froi":
                pr = paired(-d["froi"].values, -d["comp"].values)     # more negative in the fROI = positive
                rec.update({f"neg_{k}": v for k, v in pr.items()})
            t3.append(rec)
    pd.DataFrame(t3).to_csv(os.path.join(RD, "LOC_T3.csv"), index=False)

    # ---------- T2 (exploratory): confidence (C3) and change - No Change (E3) in the C1 fROIs, two-sided
    t2 = []
    for con, setname in (("C3", "T2conf"), ("E3changejoint", "T2chg")):
        subs = S[setname]; items, meta = [], []
        for p in T1_PARCELS:
            d = get(X, subs, 4, con, prim[("C1", p)], 0.10, p)
            for reg, c, vc in (("froi", "froi", "froi_var"), ("complement", "comp", "comp_var"), ("parcel", "parcel_mean", "parcel_var")):
                items.append(((con, p, reg), d[c].values, d[vc].values))
            meta.append((p, d))
        fl = flame_batch(setname, items, subs, con)
        for p, d in meta:
            for reg, c in (("froi", "froi"), ("complement", "comp"), ("parcel", "parcel_mean")):
                y = d[c].values; F = fl[(con, p, reg)]; tt = stats.ttest_1samp(y, 0)
                p2 = float(2 * stats.norm.sf(abs(F["flame_z"])))
                rec = dict(test="T2", contrast=con, N=len(y), parcel=p, froi_used=prim[("C1", p)], region=reg,
                           mean=y.mean(), dz=dz(y), t=tt.statistic, flame_z=F["flame_z"], flame_p_two=p2,
                           flame_p_two_bonf3=min(1.0, 3 * p2), BF10=jzs_bf10(float(tt.statistic), len(y), R_BF))
                rec["BF01"] = 1 / rec["BF10"]
                if reg == "froi":
                    pr = paired(d["froi"].values, d["comp"].values)
                    rec.update(diff_dz=pr["diff_dz"], diff_dz_ci_lo=pr["diff_dz_ci_lo"], diff_dz_ci_hi=pr["diff_dz_ci_hi"])
                t2.append(rec)
    pd.DataFrame(t2).to_csv(os.path.join(RD, "LOC_T2.csv"), index=False)

    # ---------- compact sensitivity / exploratory-rule tables
    sens = []
    for r in T1[T1.region == "froi"].itertuples():
        par = T1[(T1.spec == r.spec) & (T1.set == r.set) & (T1.parcel == r.parcel) & (T1.region == "parcel")].iloc[0]
        sens.append(dict(test="T1", spec=r.spec, N=r.N, parcel=r.parcel, froi_used=r.froi_used,
                         froi_is_gate_primary=r.froi_is_gate_primary, froi_dz=r.dz, froi_flame_z=r.flame_z,
                         froi_p_bonf3=r.flame_p_one_bonf3, diff_dz=r.diff_dz, diff_dz_ci_lo=r.diff_dz_ci_lo,
                         diff_dz_ci_hi=r.diff_dz_ci_hi, froi_BF01=r.BF01, froi_BF0plus=r.BF0dir,
                         parcel_dz=par.dz, parcel_flame_z=par.flame_z, parcel_p_bonf3=par.flame_p_one_bonf3,
                         parcel_BF01=par.BF01, supported=r.T1_supported))
    for r in T3B[T3B.region == "froi"].itertuples():
        sens.append(dict(test="T3b", spec=r.spec, N=r.N, parcel=SPHERE, froi_used=r.froi_used,
                         froi_r=r.r, froi_slope_t=r.slope_t, froi_p_perm_one=r.p_perm_one, froi_BF10=r.bf10,
                         diff_slope=r.diff_slope, diff_slope_ci_lo=r.diff_slope_ci_lo, diff_slope_ci_hi=r.diff_slope_ci_hi,
                         supported=r.T3b_supported))
    SS = pd.DataFrame(sens)
    expl = SS[SS.spec.str.match(r"^C[23]_")]
    expl.to_csv(os.path.join(RD, "LOC_exploratory_rules.csv"), index=False)
    SS[~SS.spec.str.match(r"^C[23]_")].to_csv(os.path.join(RD, "LOC_stageB_sensitivity.csv"), index=False)
    json.dump({k: v for k, v in S.items()}, open(os.path.join(RD, "LOC_stageB_subjects.json"), "w"), indent=1)
    print("transfer: done", flush=True)


def posthoc():
    """POST HOC, descriptive (deviations log, after T1 was seen): the non-primary individual C1 fROI in
    AInsR under the T1 sensitivity specifications. Cannot change the T1 verdict."""
    S, _ = sets()
    X = pd.read_csv(os.path.join(RD, "LOC_froi_extract.csv"))
    p = "caren_AInsR"; rows = []
    for setname, cases in (("T1", [("individual_sm4", 4, 0.10), ("individual_sm0", 0, 0.10), ("individual_sm6", 6, 0.10),
                                   ("individual_frac05", 4, 0.05), ("individual_frac20", 4, 0.20)]),
                           ("T1nov1", [("individual_nov1", 4, 0.10)])):
        subs = S[setname]; items, meta = [], []
        for name, sm, f in cases:
            d = get(X, subs, sm, "C2joint", "C1", f, p)
            items.append((name, d.froi.values, d.froi_var.values)); meta.append((name, sm, f, d))
        fl = flame_batch(f"posthoc_{setname}", items, subs, "C2joint")
        for name, sm, f, d in meta:
            y = d.froi.values; F = fl[name]; pz = float(stats.norm.sf(F["flame_z"]))
            rec = dict(analysis="POST HOC descriptive", spec=name, N=len(y), parcel=p, sm=sm, frac=f, dz=dz(y),
                       flame_z=F["flame_z"], flame_p_one_bonf3=min(1.0, 3 * pz), **bf_cols(y, +1))
            rec.update(paired(y, d.comp.values)); rows.append(rec)
    pd.DataFrame(rows).to_csv(os.path.join(RD, "LOC_T1_posthoc_AInsR_individual.csv"), index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    {"define": define, "transfer": transfer, "posthoc": posthoc}[sys.argv[1]]()
