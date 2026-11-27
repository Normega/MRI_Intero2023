#!/usr/bin/env python
"""
intero2023_loc_behaviour.py -- localizer stage C1 (OSF p5ztc, plan §4): the EXPLORATORY
behavioural and physiological measures. Everything this script writes is exploratory.

Definitions: Preregistration/localizer_deviations.md, entry "Stage C1" (logged and
committed in a247228 before any stage-C value was computed). In brief:
  tracking PLV   presses -> belt breath phase in compliant Intero blocks; participant index =
                 |mean of block resultant vectors| (Indices 10); chance = 200 circular shifts
                 per block with the aggregate recomputed per shift
  pair score     the step-1 cross-validated, lead-centred pair rule (post hoc, disclosed)
  pacer PLV      localizer Extero / Match blocks vs the pacer, null = Rest/Self/Intero blocks vs a
                 virtual pacer (step-1 block PLVs); change-detection paced epochs vs the actual
                 pacer trace, null = free epochs vs a virtual pacer at the run's baseline rate,
                 plus an added length-matched mismatched-pacer null
  relations      Spearman (permutation p), Pearson, bootstrap CI, JZS BF10; BH-FDR over 18
  H1 moderation  breathing-run adherence as a group covariate on swc7p C1 (paced - free), the
                 swc7p H6 model (nilearn permutation) in each H1 mask, plus mask means

Reuses: intero2023_localizer_qc (physio_segments, acq_parts, acq_bold_blocks, grade_block,
cut_segment, bandpass, read_log, icc31), intero2023_explore_maia_conf (jzs_corr_bf10),
intero2023_group (ROI, ols_t, peak_mni, roi_mean, subject_sets), intero2023_glm_common (run_list).

  python Scripts/intero2023_loc_behaviour.py measures    # reads raw .acq (about 2 min); caches
  python Scripts/intero2023_loc_behaviour.py stats       # tables
  python Scripts/intero2023_loc_behaviour.py h1mod       # H1 moderation (imaging, permutation)
  python Scripts/intero2023_loc_behaviour.py figure
Tables: Analysis/results/EXPLORATORY_LOC_C1_*.csv. Cache: derivatives/localizer_qc_cache/.
No wallclock time is written to any table.
"""
import os, sys, glob, re, pickle, json, warnings
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_localizer_qc as qc
warnings.simplefilter("ignore")

REPO = qc.REPO
RD = os.path.join(REPO, "Analysis", "results")
FIG = os.path.join(REPO, "Figures", "localizer")
CACHE = qc.CACHE
PK = os.path.join(CACHE, "stageC1_measures.pkl")
PFX = "EXPLORATORY_LOC_C1_"
SKIP = qc.ENTRAIN_SKIP_S            # 3 s dropped at the start of every block / epoch
MIN_WIN_S = 8.0                     # change-detection epochs need >= 8 s of window
N_SHIFT, SHIFT_SEED = qc.N_CHANCE, 20261128
N_PERM, N_BOOT, PERM_SEED, BOOT_SEED = 10000, 10000, 20261126, 20261127
PARCELS = ["caren_AInsL", "caren_AInsR", "caren_dACC", "caren_PIns", "ieat_ACCsphere"]
H1_MASKS = ["caren_PIns", "caren_COper", "caren_SM"]


# =========================================================================== measures
def press_vector(ph, fs, pt, pk):
    """Block resultant vector of presses on belt phase (the step-1 press_lock): psi = phase at a
    '2' press, phase - pi at a '1' press, so psi = 0 on the press's own belt onset."""
    ii = np.clip((pt * fs).astype(int), 0, len(ph) - 1)
    return np.exp(1j * (ph[ii] - np.where(pk == 1, np.pi, 0.0))).mean()


def localizer_tracking(blocks, presses, seg, polarity, U):
    """Per compliant Intero block: observed resultant vector and N_SHIFT circularly shifted ones."""
    from scipy.signal import hilbert
    ok = U.set_index(["subject", "run"])
    rng = np.random.default_rng(SHIFT_SEED)
    rows = []
    for (sub, run), sg in sorted(seg.items()):
        if ok.loc[(sub, run), "belt_grade"] not in ("good", "degraded"):
            continue
        fs = sg["fs"]; pol = polarity.get(sub, 0)
        xf = qc.bandpass(sg["belt"], fs)
        ph = np.angle(hilbert(xf * (pol if pol else 1)))
        B = blocks[(blocks.subject == sub) & (blocks.run == run) & (blocks.trial_type == "intero")]
        for b in B.itertuples():
            if not (b.n_presses >= qc.MIN_PRESSES):
                continue
            P = presses[(presses.subject == sub) & (presses.run == run) & (presses.block_index == b.block_index)]
            pt, pk = P.t.values, P.key.values
            z = press_vector(ph, fs, pt, pk)
            zs = np.empty(N_SHIFT, complex)
            for i in range(N_SHIFT):
                sh = rng.uniform(5, b.duration - 5)
                zs[i] = press_vector(ph, fs, b.onset + np.mod(pt - b.onset + sh, b.duration), pk)
            rows.append(dict(subject=sub, run=run, block_index=b.block_index, z=z, z_shift=zs))
    return pd.DataFrame(rows)


def log_attempts_cd(pid):
    """Every change-detection attempt in a participant's logs: (run k, first-trigger wallclock,
    set of pacer-row times: 'ChangeTrial' in v2, 'Fixation' in v1). The wallclock is used only
    to find the BOLD block and is never written out."""
    out = []
    for f in sorted(glob.glob(os.path.join(qc.LOGDIR, "Log.Intero2023.P*.Session*.csv"))):
        if int(re.search(r"\.P(\d+)\.", f).group(1)) != pid:
            continue
        d = qc.read_log(f)
        if d is None:
            continue
        for k in (1, 2):
            L = d[d.Run == f"ChangeDetection{k}"]
            pr = list(L.index[L.Description == "PulseReceived"])
            for a, i0 in enumerate(pr):
                i1 = pr[a + 1] if a + 1 < len(pr) else L.index.max() + 1
                A = L.loc[i0:i1 - 1]
                out.append(dict(k=k, tw=pd.to_datetime(A.RunStartTime.iloc[0], errors="coerce"),
                                t=np.sort(A[A.Description.isin(["ChangeTrial", "Fixation"])].t.dropna().values)))
    return out


def piecewise_phase(tp):
    return lambda t: np.interp(t, tp, np.pi * np.arange(len(tp)))


def plv(phb, php):
    return float(np.abs(np.exp(1j * (phb - php)).mean()))


def breathing_runs():
    """Change-detection pacer PLV per epoch on the 90 swc7p runs."""
    from scipy.signal import hilbert
    import intero2023_glm_common as g
    man = g.run_list()
    rr = pd.read_csv(os.path.join(REPO, "Data", "resp_runs.csv"))
    parts = qc.acq_parts()
    runs, eps = [], []
    for sub in sorted(man.subject.unique()):
        pid = int(sub[4:])
        atts = log_attempts_cd(pid)
        bold = qc.acq_bold_blocks(parts[parts.pid == pid])
        for m in man[man.subject == sub].itertuples():
            ent = m.bids_entity; k = int(ent[-1])
            pt = pd.read_csv(os.path.join(qc.EVDIR, f"{sub}_{ent}_pacertrace.tsv"), sep="\t", na_values="n/a")
            ev = pd.read_csv(os.path.join(qc.EVDIR, f"{sub}_{ent}_events.tsv"), sep="\t", na_values="n/a")
            rec = dict(subject=sub, bids_entity=ent, n_vol=int(m.n_vol))
            g_ = rr[(rr.subject == sub) & (rr.bids_entity == ent)]
            rec["belt_grade"] = g_.belt_grade.iloc[0] if len(g_) else "no_belt_record"
            rec["adhere_r"] = g_.adhere_r.iloc[0] if len(g_) else np.nan
            # attempt of record = the one whose ChangeTrial rows contain the pacer-trace onsets
            best, score = None, 0.0
            for a in atts:
                if not len(a["t"]):
                    continue
                j = np.searchsorted(a["t"], pt.onset.values).clip(0, len(a["t"]) - 1)
                jm = (j - 1).clip(0)
                dmin = np.minimum(np.abs(a["t"][j] - pt.onset.values), np.abs(a["t"][jm] - pt.onset.values))
                sc = float((dmin < 1e-3).mean())
                if sc > score:
                    best, score = a, sc
            rec.update(log_match_frac=score, log_run_k=best["k"] if best else np.nan)
            if best is None or score < 0.99 or not bold or pd.isna(best["tw"]):
                rec["grade_time"] = "no_log_match" if (best is None or score < 0.99) else "no_bold_block"
                runs.append(rec); continue
            b, dt, margin, pdel, gr = qc.grade_block(bold, best["tw"], m.n_vol)
            rec.update(dt_s=round(dt, 3), margin_s=round(margin, 1), n_pulses=b["n_pulses"],
                       pulse_delta=pdel, grade_time=gr)
            rec["usable_phase"] = (rec["belt_grade"] in ("good", "degraded")) and gr in ("exact", "good") and pdel == 0
            runs.append(rec)
            if not rec["usable_phase"]:
                continue
            belt, fs = qc.cut_segment(b, m.n_vol)
            xf = qc.bandpass(belt, fs); t = np.arange(len(xf)) / fs
            ph = np.angle(hilbert(xf))
            hc_base = float(pt[pt.pulse_index == 2].period_s.median())
            paced = ev[ev.trial_type == "paced_breathing"]
            traces = {e: np.sort(d.onset.values) for e, d in pt.groupby("epoch_index")}
            for e in paced.itertuples():
                tp = traces.get(e.epoch_index)
                if tp is None or len(tp) < 3:
                    continue
                lo, hi = max(e.onset + SKIP, tp[0]), min(e.onset + e.duration, tp[-1])
                w = (t >= lo) & (t < hi)
                if hi - lo < MIN_WIN_S:
                    continue
                r_ = plv(ph[w], piecewise_phase(tp)(t[w]))
                # added length-matched null: the other paced epochs' traces, aligned at the first pulse
                mm = []
                for e2, tp2 in traces.items():
                    if e2 == e.epoch_index or len(tp2) < 3:
                        continue
                    tq = tp2 - tp2[0] + tp[0]
                    hi2 = min(hi, tq[-1]); w2 = (t >= lo) & (t < hi2)
                    if hi2 - lo >= MIN_WIN_S:
                        mm.append(plv(ph[w2], piecewise_phase(tq)(t[w2])))
                eps.append(dict(subject=sub, bids_entity=ent, kind="paced", epoch_index=e.epoch_index,
                                change_direction=e.change_direction, salience=e.salience, window_s=hi - lo,
                                plv=r_, plv_mismatched=np.mean(mm) if mm else np.nan, n_mismatched=len(mm)))
            for e in ev[ev.trial_type == "free_breathing"].itertuples():
                lo, hi = e.onset + SKIP, e.onset + e.duration
                w = (t >= lo) & (t < min(hi, t[-1]))
                if w.sum() / fs < MIN_WIN_S:
                    continue
                eps.append(dict(subject=sub, bids_entity=ent, kind="free", epoch_index=e.epoch_index,
                                window_s=w.sum() / fs, plv=plv(ph[w], np.pi * (t[w] - e.onset) / hc_base)))
            runs[-1]["baseline_halfcycle_s"] = hc_base
        print(f"  {sub}: {sum(r['subject'] == sub for r in runs)} runs", flush=True)
    return pd.DataFrame(runs), pd.DataFrame(eps)


def measures():
    x, blocks, presses, U, seg, POL, polarity = qc.physio_segments()
    TKZ = localizer_tracking(blocks, presses, seg, polarity, U)
    BRUN, BEP = breathing_runs()
    pickle.dump(dict(TKZ=TKZ, BRUN=BRUN, BEP=BEP, polarity=polarity), open(PK, "wb"))
    print(f"[measures] tracking blocks {len(TKZ)}; change-detection runs {len(BRUN)}, epochs {len(BEP)}")


# =========================================================================== statistics
def boot_ci_rho(x, y, rng):
    n = len(x); idx = rng.integers(0, n, (N_BOOT, n))
    rs = np.array([stats.spearmanr(x[i], y[i]).statistic for i in idx])
    return np.nanpercentile(rs, 2.5), np.nanpercentile(rs, 97.5)


def corr_row(xs, ys, **lab):
    from intero2023_explore_maia_conf import jzs_corr_bf10
    d = pd.concat([xs.rename("x"), ys.rename("y")], axis=1).dropna()
    n = len(d); rec = dict(lab, n=n)
    if n < 6:
        return rec
    xv, yv = d.x.values.astype(float), d.y.values.astype(float)
    rng = np.random.default_rng(PERM_SEED)
    rho = stats.spearmanr(xv, yv).statistic
    rx, ry = stats.rankdata(xv), stats.rankdata(yv)
    rx = (rx - rx.mean()) / rx.std(); ry = (ry - ry.mean()) / ry.std()
    perm = np.array([(rx * rng.permutation(ry)).mean() for _ in range(N_PERM)])
    lo, hi = boot_ci_rho(xv, yv, np.random.default_rng(BOOT_SEED))
    r, p_r = stats.pearsonr(xv, yv)
    rec.update(spearman_rho=rho, p_perm=(1 + (np.abs(perm) >= abs(rho) - 1e-12).sum()) / (N_PERM + 1),
               rho_ci_lo=lo, rho_ci_hi=hi, pearson_r=r, pearson_p=p_r, bf10_jzs=jzs_corr_bf10(r, n, 1.0))
    return rec


def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = np.empty(m); prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, p[o[i]] * m / (i + 1)); q[o[i]] = prev
    return q


def mwu_greater(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return stats.mannwhitneyu(a, b, alternative="greater").pvalue


def resultant_index(zs):
    return float(np.abs(np.mean(zs)))


def tracking_tables(M, usable_runs_set):
    """Tracking PLV (registered aggregate) and pair score, per run and per participant."""
    U = pd.read_csv(os.path.join(qc.OUT, "belt_usability.csv"))
    tr = pd.read_csv(os.path.join(qc.OUT, "tracking_accuracy_runs.csv"))
    _, _, _, TK, _ = pickle.load(open(os.path.join(CACHE, "physio_stage.pkl"), "rb"))
    Z = M["TKZ"]
    # check: block PLVs reproduce step 1 exactly
    chk = Z.merge(TK[["subject", "run", "block_index", "press_plv"]], on=["subject", "run", "block_index"])
    assert len(chk) == len(Z) and np.allclose(np.abs(chk.z.values), chk.press_plv.values, atol=1e-12)
    meas = tr[tr.usable_breath_timing & (tr.n_belt_breaths >= 4)][["subject", "run"]]
    meas["usable_run"] = [(s, r) in usable_runs_set for s, r in zip(meas.subject, meas.run)]
    rows = []
    for (sub, run), d in Z.groupby(["subject", "run"]):
        m = meas[(meas.subject == sub) & (meas.run == run)]
        if m.empty:
            continue
        obs = resultant_index(d.z.values)
        sh = np.abs(np.vstack(d.z_shift.values).mean(0))
        t_ = tr[(tr.subject == sub) & (tr.run == run)].iloc[0]
        rows.append(dict(subject=sub, run=run, usable_run=bool(m.usable_run.iloc[0]), n_blocks=len(d),
                         tracking_plv=obs, tracking_plv_chance_mean=sh.mean(),
                         tracking_plv_chance_p95=np.percentile(sh, 95),
                         tracking_plv_p_shift=(1 + (sh >= obs).sum()) / (N_SHIFT + 1),
                         tracking_plv_step1_meanblock=t_.press_plv,
                         pair_cv=t_.prop_correct_cv_posthoc, pair_cv_chance=t_.prop_correct_cv_chance))
    R = pd.DataFrame(rows)
    # participant level, over the usable runs (primary) and with the motion-failed run kept (check)
    P = {}
    for tag, keep in (("primary", R.usable_run), ("incl_motionfail", R.usable_run | True)):
        prow = []
        for sub, d in R[keep].groupby("subject"):
            Zs = Z.merge(d[["subject", "run"]], on=["subject", "run"])
            obs = resultant_index(Zs.z.values)
            sh = np.abs(np.vstack(Zs.z_shift.values).mean(0))
            tk = TK.merge(d[["subject", "run"]], on=["subject", "run"])
            cv = tk[tk.n_correct_cv.notna()]
            prow.append(dict(subject=sub, n_runs=len(d), n_blocks=len(Zs), tracking_plv=obs,
                             tracking_plv_chance_mean=sh.mean(), tracking_plv_chance_p95=np.percentile(sh, 95),
                             tracking_plv_p_shift=(1 + (sh >= obs).sum()) / (N_SHIFT + 1),
                             pair_cv=cv.n_correct_cv.sum() / cv.n_breaths.sum() if len(cv) else np.nan,
                             pair_cv_chance=cv.n_correct_cv_chance.sum() / cv.n_breaths.sum() if len(cv) else np.nan))
        P[tag] = pd.DataFrame(prow)
    return R, P


def loc_pacer_tables(usable_runs_set):
    """Localizer Extero entrainment and Match adherence (step-1 block PLVs), Inclusion 12 runs."""
    U = pd.read_csv(os.path.join(qc.OUT, "belt_usability.csv"))
    _, _, BR, _, _ = pickle.load(open(os.path.join(CACHE, "physio_stage.pkl"), "rb"))
    ok = U[U.usable_breath_timing][["subject", "run"]]
    ok = ok[[(s, r) in usable_runs_set for s, r in zip(ok.subject, ok.run)]]
    B = BR.merge(ok, on=["subject", "run"])
    B = B[B.plv.notna()]
    def index(d):
        null = d[d.condition.isin(["rest", "self", "intero"])].plv
        out = dict(null_plv=null.mean(), n_null=len(null))
        for c in ("extero", "match"):
            r = d[d.condition == c].plv
            out.update({f"{c}_plv": r.mean(), f"{c}_n": len(r), f"{c}_index": r.mean() - null.mean(),
                        f"{c}_p_mwu": mwu_greater(r, null)})
        return out
    R = pd.DataFrame([dict(subject=s, run=r, **index(d)) for (s, r), d in B.groupby(["subject", "run"])])
    P = pd.DataFrame([dict(subject=s, n_runs=d.run.nunique(), protocol_version=d.protocol_version.iloc[0], **index(d))
                      for s, d in B.groupby("subject")])
    return R, P, B


def cd_pacer_tables(M):
    BRUN, BEP = M["BRUN"], M["BEP"]
    def index(d):
        p, f = d[d.kind == "paced"], d[d.kind == "free"]
        return dict(paced_plv=p.plv.mean(), n_paced=len(p), free_plv=f.plv.mean(), n_free=len(f),
                    mismatched_plv=p.plv_mismatched.mean(),
                    cd_index=p.plv.mean() - f.plv.mean(), cd_index_mismatched=p.plv.mean() - p.plv_mismatched.mean(),
                    cd_p_mwu=mwu_greater(p.plv, f.plv),
                    cd_p_mismatched=(stats.wilcoxon(p.plv - p.plv_mismatched, alternative="greater").pvalue
                                     if p.plv_mismatched.notna().sum() >= 6 else np.nan))
    R = pd.DataFrame([dict(subject=s, bids_entity=e, **index(d)) for (s, e), d in BEP.groupby(["subject", "bids_entity"])])
    P = pd.DataFrame([dict(subject=s, n_runs=d.bids_entity.nunique(), **index(d)) for s, d in BEP.groupby("subject")])
    R = BRUN.drop(columns=["margin_s"], errors="ignore").merge(R, on=["subject", "bids_entity"], how="left")
    return R, P


def behaviour_breathing():
    """Detection rate, mean confidence (90 swc7p runs, model events), MAIA, adherence."""
    import intero2023_glm_common as g
    q = g.q
    man = g.run_list()
    ev = pd.concat([pd.read_csv(os.path.join(REPO, q.FL, f"{s}_{e}_desc-model_events.tsv"), sep="\t").assign(subject=s)
                    for s, e in zip(man.subject, man.bids_entity)], ignore_index=True)
    p = ev[ev.condition.str.startswith("paced")]
    ch = p[p.condition.isin(q.CHANGE_ONLY)]
    d = pd.DataFrame({"detection_rate": ch.groupby("subject").condition.apply(lambda c: (c == "paced_aware").mean()),
                      "n_change": ch.groupby("subject").size(),
                      "mean_confidence": p.groupby("subject").confidence.mean()})
    maia = pd.read_csv(os.path.join(REPO, "Data", "maia_scores.csv")).set_index("subject")
    d["maia"] = maia.maia_ieat_total.reindex(d.index)
    rr = pd.read_csv(os.path.join(REPO, "Data", "resp_runs.csv")).merge(man[["subject", "bids_entity"]])
    rr = rr[rr.belt_grade.isin(["good", "degraded"]) & rr.adhere_r.notna()]
    d["adhere_r"] = rr.groupby("subject").adhere_r.mean().reindex(d.index)
    d["adhere_n_runs"] = rr.groupby("subject").size().reindex(d.index).fillna(0).astype(int)
    return d


def c1_tables():
    """C1 (Intero - Extero) per participant: cross-run C1 in the individual C1 fROIs (L4 table),
    joint 4 mm C1 in the AInsR LOSO fROI, and joint 6 mm C1 mask means in the five parcels."""
    import nibabel as nib
    import intero2023_loc_common as L
    import intero2023_group as G
    gs = pd.read_csv(os.path.join(RD, "LOC_L4_gate_subjects.csv"))
    gs = gs[gs.rule == "C1"].pivot(index="subject", columns="parcel", values="cross_C1")
    out = gs.add_prefix("C1froi_cross_")
    s, _, _ = L.participant_sets()
    fm = np.load(os.path.join(L.OUT, "froi", "froi_masks.npz"))
    pm = nib.load(G.ROI("caren_AInsR")).get_fdata() > 0
    v = {}
    for x in s["froi"]:
        d = nib.load(L.fl_name(x, "joint", 4, "base", "C1", "effect")).get_fdata(dtype=np.float32)[pm]
        v[x] = float(d[fm[f"{x}|C1loso|0.10|caren_AInsR"]].mean())
    out["C1froi_loso_caren_AInsR"] = pd.Series(v)
    masks = {p: nib.load(G.ROI(p)).get_fdata() > 0 for p in PARCELS}
    mm = {}
    for x in s["intero"]:
        d = nib.load(L.fl_name(x, "joint", 6, "base", "C1", "effect")).get_fdata(dtype=np.float32)
        mm[x] = {f"C1mask_{p}": float(d[m].mean()) for p, m in masks.items()}
    mm = pd.DataFrame(mm).T
    # check against the stage-A L1 mask means
    l1 = pd.read_csv(os.path.join(RD, "LOC_L1_subjects.csv"))
    l1 = l1[l1.roi == "caren_AInsL"].set_index("subject").mask_mean
    assert np.allclose(mm.loc[l1.index, "C1mask_caren_AInsL"], l1, atol=1e-5)
    return out.join(mm, how="outer"), s


def run_stats():
    M = pickle.load(open(PK, "rb"))
    usable = qc.usable_runs(); usable = usable[usable.usable]
    uset = set(zip(usable.subject, usable.run))
    TR_R, TR_P = tracking_tables(M, uset)
    LP_R, LP_P, LP_B = loc_pacer_tables(uset)
    CD_R, CD_P = cd_pacer_tables(M)
    BEH = behaviour_breathing()
    C1, sets = c1_tables()

    # ---------------- participant table
    P = TR_P["primary"].set_index("subject").add_prefix("trk_")
    P = P.join(LP_P.set_index("subject").add_prefix("loc_"), how="outer")
    P = P.join(CD_P.set_index("subject").add_prefix("cd_"), how="outer")
    P = P.join(BEH, how="outer").join(C1, how="outer")
    P["in_intero_set"] = P.index.isin(sets["intero"]); P["in_froi_set"] = P.index.isin(sets["froi"])
    P.index.name = "subject"
    P.to_csv(os.path.join(RD, PFX + "participants.csv"))
    # run tables (no wallclock: only grades, pulse counts and the dt difference)
    TR_R.to_csv(os.path.join(RD, PFX + "tracking_runs.csv"), index=False)
    LP_R.to_csv(os.path.join(RD, PFX + "localizer_pacer_runs.csv"), index=False)
    CD_R.to_csv(os.path.join(RD, PFX + "changedetection_pacer_runs.csv"), index=False)

    # ---------------- above chance (the PI's question) and distributions
    rows = []
    def summ(measure, level, obs, null, above, pmask=None, extra=None):
        ok = obs.notna() & null.notna()
        o, n_ = obs[ok], null[ok]
        rec = dict(measure=measure, level=level, n=int(ok.sum()), median=o.median(), q25=o.quantile(.25),
                   q75=o.quantile(.75), min=o.min(), max=o.max(), null_median=n_.median(),
                   index_median=(o - n_).median(), n_obs_gt_null=int((o > n_).sum()),
                   n_above_chance_p05=int(above[ok].sum()), prop_above_chance_p05=float(above[ok].mean()),
                   group_wilcoxon_p=stats.wilcoxon(o - n_, alternative="greater").pvalue if ok.sum() >= 6 else np.nan,
                   group_dz=(o - n_).mean() / (o - n_).std(ddof=1))
        rec.update(extra or {})
        rows.append(rec)
    Tp, Tr = TR_P["primary"], TR_R[TR_R.usable_run]
    for lev, d in (("participant", Tp), ("run", Tr)):
        summ("tracking_plv", lev, d.tracking_plv, d.tracking_plv_chance_mean, d.tracking_plv_p_shift < .05,
             extra=dict(above_rule="observed > 95th pct of 200 circular shifts (p < .05)"))
        summ("pair_cv_posthoc", lev, d.pair_cv, d.pair_cv_chance, d.pair_cv > d.pair_cv_chance,
             extra=dict(above_rule="observed > shift mean (step-1 rule; no shift distribution kept)"))
    Ti = TR_P["incl_motionfail"]
    summ("tracking_plv_incl_motionfail_run", "participant", Ti.tracking_plv, Ti.tracking_plv_chance_mean,
         Ti.tracking_plv_p_shift < .05, extra=dict(above_rule="as tracking_plv; sub-108 run 1 kept"))
    for c in ("extero", "match"):
        for lev, d in (("participant", LP_P), ("run", LP_R)):
            summ(f"localizer_{c}_plv", lev, d[f"{c}_plv"], d.null_plv, d[f"{c}_p_mwu"] < .05,
                 extra=dict(above_rule="one-sided Mann-Whitney, real-pacer vs null blocks, p < .05"))
    for lev, d in (("participant", CD_P), ("run", CD_R[CD_R.paced_plv.notna()])):
        summ("changedetection_paced_plv", lev, d.paced_plv, d.free_plv, d.cd_p_mwu < .05,
             extra=dict(above_rule="one-sided Mann-Whitney, paced vs free epochs (registered null), p < .05"))
        summ("changedetection_paced_plv_vs_mismatched", lev, d.paced_plv, d.mismatched_plv, d.cd_p_mismatched < .05,
             extra=dict(above_rule="added null: other epochs' pacer traces; one-sided Wilcoxon over epochs, p < .05"))
    CH = pd.DataFrame(rows)
    # both tasks, the same participants
    both = P[P.loc_match_index.notna() & P.cd_cd_index.notna()]
    CH = pd.concat([CH, pd.DataFrame([dict(
        measure="breathing_along_both_tasks", level="participant", n=len(both),
        n_above_chance_p05=int(((both.loc_match_p_mwu < .05) & (both.cd_cd_p_mwu < .05)).sum()),
        prop_above_chance_p05=float(((both.loc_match_p_mwu < .05) & (both.cd_cd_p_mwu < .05)).mean()),
        n_obs_gt_null=int(((both.loc_match_index > 0) & (both.cd_cd_index > 0)).sum()),
        above_rule="Match and change-detection paced both above chance (p < .05)")])], ignore_index=True)
    CH.to_csv(os.path.join(RD, PFX + "above_chance.csv"), index=False)

    # ---------------- reliability
    rel = []
    def icc_row(measure, W):
        W = W.dropna()
        if len(W) < 5:
            return
        a, b = W.iloc[:, 0].values, W.iloc[:, 1].values
        rel.append(dict(measure=measure, n_both_runs=len(W), icc31=qc.icc31(a, b), pearson_r=stats.pearsonr(a, b)[0],
                        spearman_rho=stats.spearmanr(a, b).statistic))
    for m in ("tracking_plv", "tracking_plv_step1_meanblock", "pair_cv"):
        icc_row(m, Tr.pivot(index="subject", columns="run", values=m))
    icc_row("tracking_plv_incl_motionfail_run", TR_R.pivot(index="subject", columns="run", values="tracking_plv"))
    for c in ("extero", "match"):
        icc_row(f"localizer_{c}_index", LP_R.pivot(index="subject", columns="run", values=f"{c}_index"))
    icc_row("changedetection_index", CD_R.pivot(index="subject", columns="bids_entity", values="cd_index"))
    icc_row("changedetection_paced_plv", CD_R.pivot(index="subject", columns="bids_entity", values="paced_plv"))
    icc_row("adhere_r", CD_R[CD_R.belt_grade.isin(["good", "degraded"])].pivot(index="subject", columns="bids_entity", values="adhere_r"))
    REL = pd.DataFrame(rel)
    REL.to_csv(os.path.join(RD, PFX + "reliability.csv"), index=False)

    # ---------------- manipulation check (from step 1; not recomputed)
    b = pd.read_csv(os.path.join(qc.OUT, "breathing_by_condition.csv"))
    MC = b[(b.level.isin(["group", "paired", "entrainment"])) & (b.protocol_version.isin(["all", "v2"]))].dropna(axis=1, how="all")
    MC.to_csv(os.path.join(RD, PFX + "manipulation_check.csv"), index=False)

    # ---------------- relations
    trk = P.trk_tracking_plv
    rel_rows = []
    def add(family, x, y, xname, yname, subset=None, note=""):
        xx, yy = (x, y) if subset is None else (x[x.index.isin(subset)], y[y.index.isin(subset)])
        rel_rows.append(corr_row(xx, yy, family=family, x=xname, y=yname, note=note))
    add("FDR18", trk, P.maia, "tracking_plv", "maia_ieat_total")
    add("FDR18", trk, P.detection_rate, "tracking_plv", "breathing_detection_rate")
    add("FDR18", trk, P.mean_confidence, "tracking_plv", "breathing_mean_confidence")
    add("FDR18", trk, P.adhere_r, "tracking_plv", "breathing_adhere_r")
    for p in PARCELS:
        add("FDR18", trk, P[f"C1froi_cross_{p}"], "tracking_plv", f"C1_froi_crossrun_{p}", sets["froi"],
            note="individual C1 fROI, cross-run" + ("; gate-primary is LOSO here" if p == "caren_AInsR" else ""))
    for p in PARCELS:
        add("FDR18", trk, P[f"C1mask_{p}"], "tracking_plv", f"C1_mask_joint6mm_{p}", sets["intero"])
    add("FDR18", P.loc_match_index, P.adhere_r, "localizer_match_index", "breathing_adhere_r")
    # sensitivity / added (outside the FDR family)
    add("sensitivity", trk, P.cd_cd_index, "tracking_plv", "changedetection_index", note="PLV-scale adherence for adhere_r")
    add("sensitivity", P.loc_match_index, P.cd_cd_index, "localizer_match_index", "changedetection_index", note="same scale, both tasks")
    add("sensitivity", P.loc_extero_index, P.cd_cd_index, "localizer_extero_index", "changedetection_index")
    add("sensitivity", P.loc_match_index, P.cd_cd_index_mismatched, "localizer_match_index", "changedetection_index_mismatched")
    add("sensitivity", P.trk_pair_cv, P.maia, "pair_cv", "maia_ieat_total")
    add("sensitivity", P.trk_pair_cv, P.detection_rate, "pair_cv", "breathing_detection_rate")
    add("sensitivity", P.trk_pair_cv, P.mean_confidence, "pair_cv", "breathing_mean_confidence")
    add("sensitivity", P.trk_pair_cv, P.adhere_r, "pair_cv", "breathing_adhere_r")
    add("sensitivity", trk, P.C1froi_loso_caren_AInsR, "tracking_plv", "C1_froi_loso_caren_AInsR", sets["froi"],
        note="gate-primary AInsR fROI (LOSO), joint 4 mm")
    add("descriptive", P.cd_cd_index, P.adhere_r, "changedetection_index", "breathing_adhere_r",
        note="validity of the PLV-scale measure against the registered one")
    add("descriptive", P.cd_paced_plv, P.adhere_r, "changedetection_paced_plv", "breathing_adhere_r")
    add("descriptive", P.trk_tracking_plv, P.loc_match_index, "tracking_plv", "localizer_match_index")
    add("descriptive", P.trk_tracking_plv, P.loc_extero_index, "tracking_plv", "localizer_extero_index")
    # POST HOC (logged after the first tables were seen): the registered change-detection index is
    # dominated by its free-epoch null term, so the raw paced PLV and the length-matched index are added
    ph = "post hoc (deviation log): "
    add("posthoc", P.cd_cd_index_mismatched, P.adhere_r, "changedetection_index_mismatched", "breathing_adhere_r", note=ph + "length-matched null")
    add("posthoc", P.cd_cd_index, P.cd_free_plv, "changedetection_index", "changedetection_free_plv", note=ph + "index vs its own null term")
    add("posthoc", P.cd_free_plv, P.adhere_r, "changedetection_free_plv", "breathing_adhere_r", note=ph)
    add("posthoc", P.loc_null_plv, P.cd_free_plv, "localizer_null_plv", "changedetection_free_plv", note=ph + "null PLV as a participant trait")
    add("posthoc", P.loc_match_index, P.loc_null_plv, "localizer_match_index", "localizer_null_plv", note=ph + "index vs its own null term")
    add("posthoc", P.loc_match_plv, P.cd_paced_plv, "localizer_match_plv", "changedetection_paced_plv", note=ph + "raw PLV, both tasks")
    add("posthoc", P.loc_match_plv, P.adhere_r, "localizer_match_plv", "breathing_adhere_r", note=ph + "raw Match PLV")
    add("posthoc", trk, P.cd_paced_plv, "tracking_plv", "changedetection_paced_plv", note=ph + "raw paced PLV")
    RE = pd.DataFrame(rel_rows)
    # the three H1-moderation mask-mean rows join the family when h1mod has run
    hm = os.path.join(RD, PFX + "H1_moderation.csv")
    if os.path.exists(hm):
        H = pd.read_csv(hm)
        H = H[(H.kind == "mask_mean") & (H.moderator == "adhere_r")]
        RE = pd.concat([RE, H.assign(family="FDR18", x="breathing_adhere_r", y="H1_C1_" + H.roi)[
            [c for c in RE.columns if c in H.columns or c in ("family", "x", "y")]]], ignore_index=True)
    fam = RE.family == "FDR18"
    RE.loc[fam, "q_bh_fdr18"] = bh(RE.loc[fam, "p_perm"].values)
    RE["fdr_family_size"] = np.where(fam, fam.sum(), np.nan)
    RE.to_csv(os.path.join(RD, PFX + "relations.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(CH.round(3).to_string()); print(REL.round(3).to_string())
    print(RE[["family", "x", "y", "n", "spearman_rho", "p_perm", "rho_ci_lo", "rho_ci_hi", "pearson_r", "bf10_jzs", "q_bh_fdr18"]].round(3).to_string())


# =========================================================================== H1 moderation
def h1mod():
    import nibabel as nib
    from nilearn.glm.second_level import non_parametric_inference
    import intero2023_group as G
    P = pd.read_csv(os.path.join(RD, PFX + "participants.csv")).set_index("subject")
    sets, _, _ = G.subject_sets()
    rows = []
    for mod, col in (("adhere_r", "adhere_r"), ("changedetection_index", "cd_cd_index"),
                     ("changedetection_paced_plv_posthoc", "cd_paced_plv")):
        subs = [x for x in sets["all"] if pd.notna(P[col].get(x, np.nan)) and os.path.exists(G.fl(x, "base", "C1", "effect"))]
        imgs = [nib.load(G.fl(x, "base", "C1", "effect")) for x in subs]
        xv = P.loc[subs, col].astype(float).values
        dm = pd.DataFrame({"intercept": 1.0, "mod": xv - xv.mean()})
        for roi in H1_MASKS:
            res = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="mod", mask=G.ROI(roi),
                                           n_perm=G.N_PERM, two_sided_test=True, random_state=G.SEED, n_jobs=G.N_JOBS)
            rmask = nib.load(G.ROI(roi)).get_fdata() > 0
            lpv = res.get_fdata()[rmask]; pmin = float(10 ** -np.nanmax(lpv))
            tmap = G.ols_t(imgs, dm, "mod", G.ROI(roi))
            coord, _ = G.peak_mni(nib.Nifti1Image(np.abs(tmap.get_fdata()), tmap.affine), G.ROI(roi), +1)
            ijk = tuple(np.round(nib.affines.apply_affine(np.linalg.inv(tmap.affine), coord)).astype(int))
            rows.append(dict(kind="voxel_maxt", moderator=mod, roi=roi, n=len(subs), peak_t_signed=float(tmap.get_fdata()[ijk]),
                             peak_mni=coord, p_fwe_min=pmin, p_fwe_bonf3=min(1.0, 3 * pmin),
                             n_vox_p05_bonf3=int((lpv > -np.log10(0.05 / 3)).sum())))
            means = pd.Series(G.roi_mean(imgs, G.ROI(roi)), index=subs)
            X = np.column_stack([np.ones(len(subs)), xv - xv.mean()])
            beta, *_ = np.linalg.lstsq(X, means.values, rcond=None)
            rec = corr_row(P.loc[subs, col], means, kind="mask_mean", moderator=mod, roi=roi)
            rec.update(slope=beta[1], C1_mean=means.mean(), C1_dz=means.mean() / means.std(ddof=1))
            rows.append(rec)
            print(rows[-2]); print(rec, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RD, PFX + "H1_moderation.csv"), index=False)


# =========================================================================== figure
INK, INK2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    P = pd.read_csv(os.path.join(RD, PFX + "participants.csv")).set_index("subject")
    RE = pd.read_csv(os.path.join(RD, PFX + "relations.csv"))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(11.5, 7.4))
    gsp = fig.add_gridspec(2, 4, height_ratios=[1.05, 1], hspace=0.62, wspace=0.42, left=0.06, right=0.985, top=0.86, bottom=0.08)
    rng = np.random.default_rng(1)

    def strip(ax, cols, title):
        for i, (lab, obs, null, sig, col) in enumerate(cols):
            ok = obs.notna() & null.notna(); o, nl, s_ = obs[ok], null[ok], sig[ok]
            xj = i + rng.uniform(-0.13, 0.13, len(o))
            for a, b_, xx in zip(nl, o, xj):
                ax.plot([xx - 0.28, xx], [a, b_], color=GRID, lw=0.6, zorder=1)
            ax.scatter(xj - 0.28, nl, s=12, facecolor="none", edgecolor=INK2, lw=0.6, zorder=2)
            ax.scatter(xj[s_.values], o[s_.values], s=16, color=col, edgecolor=SURF, lw=0.5, zorder=3)
            ax.scatter(xj[~s_.values], o[~s_.values], s=16, facecolor=SURF, edgecolor=col, lw=0.9, zorder=3)
            ax.text(i, 1.03, f"{int(s_.sum())}/{len(o)} above chance", ha="center", va="bottom", color=INK, fontsize=7.6)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[0] for c in cols], fontsize=7.8, color=INK)
        ax.set_ylim(0, 1.12); ax.set_yticks([0, .25, .5, .75, 1]); ax.set_xlim(-0.75, len(cols) - 0.45)
        ax.set_ylabel("phase-locking value"); ax.set_title(title, loc="left", color=INK, fontsize=9.2, pad=16)
        ax.spines[["top", "right"]].set_visible(False); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)

    ax = fig.add_subplot(gsp[0, 0:2])
    strip(ax, [("Breath tracking\n(presses → belt)", P.trk_tracking_plv, P.trk_tracking_plv_chance_mean,
                P.trk_tracking_plv_p_shift < .05, BLUE),
               ("Localizer Match\n(belt → pacer)", P.loc_match_plv, P.loc_null_plv, P.loc_match_p_mwu < .05, ORANGE),
               ("Localizer Extero\n(unintended)", P.loc_extero_plv, P.loc_null_plv, P.loc_extero_p_mwu < .05, YELLOW),
               ("Change-detection\npaced epochs", P.cd_paced_plv, P.cd_free_plv, P.cd_cd_p_mwu < .05, AQUA)],
          "A  Each participant against their own chance level")
    ax.text(-0.72, -0.30, "Open grey circle = participant's chance / null; coloured = observed (filled: p < .05 vs own null).",
            transform=ax.transData, fontsize=7.2, color=INK2, clip_on=False)

    ax = fig.add_subplot(gsp[0, 2])
    d = P[["cd_paced_plv", "adhere_r"]].dropna()
    ax.scatter(d.adhere_r, d.cd_paced_plv, s=16, color=AQUA, edgecolor=SURF, lw=0.5)
    r = RE[(RE.x == "changedetection_paced_plv") & (RE.y == "breathing_adhere_r")].iloc[0]
    r2 = RE[(RE.x == "changedetection_index") & (RE.y == "breathing_adhere_r")].iloc[0]
    ax.set_xlabel("registered adherence (adhere_r)"); ax.set_ylabel("paced-epoch PLV")
    ax.set_title("B  Change-detection: PLV vs\n     registered adherence", loc="left", color=INK, fontsize=9.2)
    ax.text(0.03, 0.97, f"paced PLV: ρ = {r.spearman_rho:.2f}, n = {int(r.n)}\nPLV − free null: ρ = {r2.spearman_rho:.2f}",
            transform=ax.transAxes, va="top", color=INK, fontsize=7.4)
    ax.set_ylim(0.35, 1.02); ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(gsp[0, 3])
    d = P[["loc_match_plv", "cd_paced_plv"]].dropna()
    ax.scatter(d.loc_match_plv, d.cd_paced_plv, s=16, color=ORANGE, edgecolor=SURF, lw=0.5)
    r = RE[(RE.x == "localizer_match_plv") & (RE.y == "changedetection_paced_plv")].iloc[0]
    r2 = RE[(RE.x == "localizer_match_index") & (RE.y == "changedetection_index")].iloc[0]
    ax.set_xlabel("localizer Match PLV"); ax.set_ylabel("change-detection paced PLV")
    ax.set_title("C  Pacer coupling in both tasks\n     (same PLV scale)", loc="left", color=INK, fontsize=9.2)
    ax.text(0.03, 0.97, f"raw PLV: ρ = {r.spearman_rho:.2f}, n = {int(r.n)}\nPLV − null: ρ = {r2.spearman_rho:.2f}",
            transform=ax.transAxes, va="top", color=INK, fontsize=7.4)
    ax.set_xlim(0.3, 1.02); ax.set_ylim(0.35, 1.02); ax.spines[["top", "right"]].set_visible(False)

    sc = [("maia", "maia_ieat_total", "MAIA (IEAT total, 0–5)"),
          ("detection_rate", "breathing_detection_rate", "change-detection hit rate"),
          ("adhere_r", "breathing_adhere_r", "breathing-run adherence (adhere_r)"),
          ("C1mask_caren_PIns", "C1_mask_joint6mm_caren_PIns", "C1 Intero − Extero, PIns mask (%)")]
    for j, (col, yname, lab) in enumerate(sc):
        ax = fig.add_subplot(gsp[1, j])
        d = P[["trk_tracking_plv", col]].dropna()
        ax.scatter(d.trk_tracking_plv, d[col], s=16, color=BLUE, edgecolor=SURF, lw=0.5)
        r = RE[(RE.x == "tracking_plv") & (RE.y == yname)].iloc[0]
        ax.set_xlabel("breath-tracking PLV"); ax.set_ylabel(lab)
        ax.set_title(("D  " if j == 0 else "") + ["Tracking × MAIA", "Tracking × detection", "Tracking × adherence",
                                                   "Tracking × C1 (PIns)"][j], loc="left", color=INK, fontsize=9.2, pad=20)
        ax.text(0.0, 1.015, f"ρ {r.spearman_rho:+.2f} [{r.rho_ci_lo:+.2f}, {r.rho_ci_hi:+.2f}] · n {int(r.n)} · BF10 {r.bf10_jzs:.2f}",
                transform=ax.transAxes, va="bottom", ha="left", color=INK2, fontsize=7.0)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("EXPLORATORY (OSF p5ztc §4): breath tracking and pacer coupling in the localizer and change-detection runs",
                 x=0.06, y=0.975, ha="left", color=INK, fontsize=10.5)
    nq = int((RE.q_bh_fdr18 < .05).sum())
    fig.text(0.06, 0.935, f"Bracketed: bootstrap 95% CI on Spearman ρ. Relations surviving BH-FDR across the 18 listed relations: {nq}. "
             "B and C use raw paced PLV (post hoc; see deviation log).",
             color=INK2, fontsize=7.6)
    os.makedirs(FIG, exist_ok=True)
    out = os.path.join(FIG, "loc_C1_breath_tracking_adherence.png")
    fig.savefig(out, dpi=200); print(out)


if __name__ == "__main__":
    st = sys.argv[1:] or ["measures", "stats", "h1mod", "stats", "figure"]
    for s_ in st:
        dict(measures=measures, stats=run_stats, h1mod=h1mod, figure=figure)[s_]()
