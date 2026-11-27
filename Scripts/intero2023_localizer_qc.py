#!/usr/bin/env python
"""
intero2023_localizer_qc.py -- localizer step 1: behaviour, physiology and design
QC run BEFORE the localizer analysis plan (Preregistration/localizer_plan_draft.md)
is registered.

BLINDNESS. This script reads no localizer BOLD voxel. From localizer NIfTIs it
reads the header only (shape[3] = n volumes, pixdim[4] = TR); it never calls
get_fdata(). From fMRIPrep it reads the confounds TSV (framewise_displacement)
only. It relates no localizer behavioural or physiological measure to MAIA,
to breathing-run behaviour, to demographics or to any brain measure; the
breath-tracking measure is summarised by distribution, missingness and
run-1/run-2 reliability only.

Stages (all run by default; pass names to run a subset):
  events      log attempt of record -> Data/events/<sub>_task-localizer_run-<k>_events.tsv,
              Analysis/localizer_qc/events_qc.csv
  motion      registered run-level motion gate -> Analysis/localizer_qc/motion_gate.csv
  hpf         contrast efficiency under no drift / 128 s / 200 s cosine sets, VIFs
              -> Analysis/localizer_qc/hpf_efficiency.csv, hpf_efficiency_participant.csv,
                 hpf_sweep.csv, hpf_vif.csv
  compliance  Intero presses per block -> Analysis/localizer_qc/compliance.csv
  physio      raw Biopac belt per localizer run (needs bioread==2025.5.2 and the .acq
              copies in sourcedata/physio) -> belt_usability.csv, breathing_by_condition.csv,
              tracking_accuracy_runs.csv, tracking_accuracy_reliability.csv
              ('physio --reuse' re-summarises the cached per-block tables)

Usage: python Scripts/intero2023_localizer_qc.py [events motion hpf compliance physio]
Stages depend on the earlier ones' outputs. Results: Analysis/localizer_qc/README.md.

TIMING (from InteroPilot.4.51.py and verified on every log):
  RunTime is seconds since the first scanner trigger of the run.
  Every block: Instruction row (6 s image) -> [v2: Fixation row, 2 s '+'] -> task.
  v2 (P60+): task onset = the block's first 'Trial' row (Fixation + ~2.005 s),
     except MindWandering, which writes its Trial row at the END of its 29.5 s
     wait, so its onset = Fixation row + 2.0 s.
  v1 (P51-59): no Fixation row is logged. The circle tasks (Extero, Match,
     Intero) still start ~8.24 s after their Instruction row; Self and Rest start
     ~6.22 s after it (no fixation). Derived per run, never assumed:
       Extero/Match onset = first half-cycle row - pulsetime/2, where the pacer
         writes a row at the end of each grow/shrink phase, pulsetime/2 =
         0.8 x half-cycle and the half-cycle is the median row spacing;
       Intero onset = Instruction row + the run's median Extero/Match offset
         (checked against Completion - 30 s);
       Self onset = first word row - 4.5 s (word timeout 5 s - 0.5 s);
       Rest onset = Trial row - 29.5 s.
  Task offset = the block's Completion row (Extero/Match run whole pacer cycles,
  so in v1 they last 30.1-34.4 s).
"""
import os, sys, re, glob, json, warnings
import numpy as np, pandas as pd

warnings.simplefilter("ignore")
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ROOT = os.path.abspath(os.path.join(REPO, ".."))
LOGDIR = os.path.join(ROOT, "sourcedata", "logs")
PHYSDIR = os.path.join(ROOT, "sourcedata", "physio")
BIDS = os.path.join(ROOT, "bids")
FMRIPREP = os.path.join(ROOT, "derivatives", "fmriprep")
EVDIR = os.path.join(REPO, "Data", "events")
OUT = os.path.join(REPO, "Analysis", "localizer_qc")
CACHE = os.path.join(ROOT, "derivatives", "localizer_qc_cache")   # intermediates, outside the repo
TR = 2.0
PILOTS = {0, 1, 50, 999}
V2_FIRST = 60
COND = {"ActiveIntero": "intero", "ActiveExtero": "extero", "ActiveMatch": "match",
        "SelfRating": "self", "MindWandering": "rest"}
TASKS = ["intero", "extero", "match", "self", "rest"]
MIN_PRESSES = 4            # registered-in-draft compliance rule (plan §2)
END_TOL_S = 2.0            # last event may end this far past n_vol x TR
os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)


def secs(x):
    try:
        h, m, s = str(x).split(":")
        return int(h) * 3600 + int(m) * 60 + float(s)
    except Exception:
        return np.nan


def pversion(p):
    return "v2" if p >= V2_FIRST else "v1"


# ---------------------------------------------------------------------------
# BIDS inventory (headers only)
# ---------------------------------------------------------------------------
def bids_localizer_runs():
    import nibabel as nib
    rows = []
    for f in sorted(glob.glob(f"{BIDS}/sub-*/func/*task-localizer_run-*_bold.nii.gz")):
        sub = re.search(r"(sub-\d+)", f).group(1)
        run = int(re.search(r"run-(\d)", f).group(1))
        img = nib.load(f)                      # header only; no get_fdata()
        n_vol = int(img.shape[3])
        tr = float(img.header.get_zooms()[3])
        js = json.load(open(f.replace(".nii.gz", ".json")))
        rows.append(dict(subject=sub, pid=int(sub[4:]), run=run, n_vol=n_vol, tr=tr,
                         acq_time=js.get("AcquisitionTime"),
                         protocol=js.get("ProtocolName")))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------
def read_log(f):
    if os.path.getsize(f) == 0:
        return None
    try:
        d = pd.read_csv(f, dtype=str, index_col=False)
    except Exception:
        return None
    if "Run" not in d.columns or d.empty:
        return None
    d["t"] = d.RunTime.map(secs)
    return d


def log_attempts():
    """Every localizer attempt (segment from one PulseReceived row to the next)."""
    atts = []
    for f in sorted(glob.glob(f"{LOGDIR}/Log.Intero2023.P*.Session*.csv")):
        p = int(re.search(r"\.P(\d+)\.", f).group(1))
        if p in PILOTS:
            continue
        d = read_log(f)
        if d is None:
            continue
        for run in ("Localizer1", "Localizer2"):
            L = d[d.Run == run]
            if L.empty:
                continue
            pr = list(L.index[L.Description == "PulseReceived"])
            for a, i0 in enumerate(pr):
                i1 = pr[a + 1] if a + 1 < len(pr) else L.index.max() + 1
                A = L.loc[i0:i1 - 1]
                atts.append(dict(pid=p, run=int(run[-1]), file=os.path.basename(f),
                                 attempt=a, rows=A,
                                 trigger_wall=pd.to_datetime(A.RunStartTime.iloc[0], errors="coerce")))
    return atts


def parse_blocks(A, version):
    """Blocks of one attempt -> list of dicts with instruction/task timing."""
    ins = A[(A.Description == "Instruction") & (A.Block != "StartUp")]
    blocks = []
    for j, ix in enumerate(ins.index):
        nxt = ins.index[j + 1] if j + 1 < len(ins) else A.index.max() + 1
        B = A.loc[ix:nxt - 1]
        name = B.Block.iloc[0]
        t_ins = B.t.iloc[0]
        fx = B[B.Description == "Fixation"]
        tr = B[B.Description == "Trial"]
        comp = B[B.Description == "Completion"]
        b = dict(block_name=name, trial_type=COND.get(name, name), t_ins=t_ins,
                 t_fix=fx.t.iloc[0] if len(fx) else np.nan,
                 t_comp=comp.t.iloc[0] if len(comp) else np.nan,
                 onset=np.nan, anchor="", half_cycle=np.nan, presses=None)
        tnum = pd.to_numeric(tr.TrialNum, errors="coerce")
        if name in ("ActiveExtero", "ActiveMatch"):
            h = tr[tnum > 0].t.values
            if len(h) >= 3:
                b["half_cycle"] = float(np.median(np.diff(h)))
        if name == "ActiveIntero":
            pr = tr[(tnum > 0) & tr.Response.isin(["1", "2"])]
            b["presses"] = pr[["t", "Response"]].rename(columns={"Response": "key"}).assign(
                key=lambda x: x.key.astype(int)).reset_index(drop=True)
        if version == "v2":
            if name == "MindWandering":
                b["onset"] = b["t_fix"] + 2.0; b["anchor"] = "fixation+2"
            else:
                t0 = tr[tnum == 0].t
                if len(t0):
                    b["onset"] = t0.iloc[0]; b["anchor"] = "trial0_row"
        else:
            if name in ("ActiveExtero", "ActiveMatch") and np.isfinite(b["half_cycle"]):
                h1 = tr[tnum > 0].t.values[0]
                b["onset"] = h1 - 0.8 * b["half_cycle"]; b["anchor"] = "halfcycle1-0.8hc"
            elif name == "SelfRating":
                w1 = tr[tnum > 0].t
                if len(w1):
                    b["onset"] = w1.iloc[0] - 4.5; b["anchor"] = "word1-4.5"
            elif name == "MindWandering":
                if len(tr):
                    b["onset"] = tr.t.iloc[0] - 29.5; b["anchor"] = "trial-29.5"
        blocks.append(b)
    if version == "v1":
        offs = [b["onset"] - b["t_ins"] for b in blocks
                if b["block_name"] in ("ActiveExtero", "ActiveMatch") and np.isfinite(b["onset"])]
        off = float(np.median(offs)) if offs else np.nan
        for b in blocks:
            if b["block_name"] == "ActiveIntero":
                b["onset"] = b["t_ins"] + off; b["anchor"] = "instr+run_circle_offset"
                b["intero_check_s"] = (b["t_comp"] - 30.0) - b["onset"]
    return blocks


def stage_events():
    inv = bids_localizer_runs()
    atts = log_attempts()
    qc, all_blocks = [], []
    os.makedirs(EVDIR, exist_ok=True)
    for _, r in inv.iterrows():
        ver = pversion(r.pid)
        cands = [a for a in atts if a["pid"] == r.pid and a["run"] == r.run]
        acq = None
        if r.acq_time:
            hh, mm, ss = r.acq_time.split(":")
            acq = int(hh) * 3600 + int(mm) * 60 + float(ss)
        summ = []
        for a in cands:
            bl = parse_blocks(a["rows"], ver)
            tw = a["trigger_wall"]
            dt = (tw.hour * 3600 + tw.minute * 60 + tw.second + tw.microsecond / 1e6 - acq) \
                if (acq is not None and pd.notna(tw)) else np.nan
            summ.append(dict(a=a, blocks=bl, n_blocks=len(bl),
                             n_complete=sum(np.isfinite(b["t_comp"]) for b in bl), dt=dt))
        rec = dict(subject=r.subject, run=r.run, protocol_version=ver, n_vol=r.n_vol, tr=r.tr,
                   scan_dur_s=r.n_vol * r.tr, n_log_attempts=len(cands))
        if not summ:
            rec.update(selected_file="", flag="no_log")
            qc.append(rec); continue
        # attempt of record: most completed blocks, then closest to the scan's AcquisitionTime
        summ.sort(key=lambda s: (-s["n_complete"], abs(s["dt"]) if np.isfinite(s["dt"]) else 1e9))
        s = summ[0]
        bl = s["blocks"]
        rec.update(selected_file=s["a"]["file"], selected_attempt=s["a"]["attempt"],
                   trigger_minus_acqtime_s=s["dt"], n_blocks=s["n_blocks"],
                   n_blocks_complete=s["n_complete"],
                   other_attempts="; ".join(f'{x["a"]["file"]}#{x["a"]["attempt"]}:{x["n_complete"]}blk'
                                            for x in summ[1:]))
        ev = []
        for k, b in enumerate(bl, start=1):
            if not (np.isfinite(b["onset"]) and np.isfinite(b["t_comp"])):
                continue
            ev.append(dict(onset=b["t_ins"], duration=b["onset"] - b["t_ins"], trial_type="instruction",
                           block_index=k, protocol_version=ver, n_presses="n/a", compliant="n/a"))
            npress = len(b["presses"]) if b["presses"] is not None else None
            ev.append(dict(onset=b["onset"], duration=b["t_comp"] - b["onset"], trial_type=b["trial_type"],
                           block_index=k, protocol_version=ver,
                           n_presses=npress if npress is not None else "n/a",
                           compliant=(int(npress >= MIN_PRESSES) if npress is not None else "n/a")))
            all_blocks.append(dict(subject=r.subject, run=r.run, protocol_version=ver, block_index=k,
                                   trial_type=b["trial_type"], instr_onset=b["t_ins"], onset=b["onset"],
                                   duration=b["t_comp"] - b["onset"],
                                   instr_dur=b["onset"] - b["t_ins"], anchor=b["anchor"],
                                   half_cycle_s=b["half_cycle"], n_presses=npress,
                                   intero_v1_check_s=b.get("intero_check_s", np.nan)))
            if b["presses"] is not None and len(b["presses"]):
                pp = b["presses"].assign(subject=r.subject, run=r.run, block_index=k)
                PRESSES.append(pp)
        ev = pd.DataFrame(ev)
        task = ev[ev.trial_type != "instruction"]
        last_end = float((ev.onset + ev.duration).max()) if len(ev) else np.nan
        flags = []
        if len(task) != 10: flags.append(f"{len(task)}_task_blocks")
        if len(task) and ((task.duration < 29.5) | (task.duration > 35)).any(): flags.append("task_duration")
        if np.isfinite(last_end) and last_end > r.n_vol * r.tr + END_TOL_S: flags.append("overruns_scan")
        scan_end = r.n_vol * r.tr
        lost = float(((ev.onset + ev.duration).clip(upper=None) - scan_end).clip(lower=0).max()) if len(ev) else 0.0
        counts = task.trial_type.value_counts()
        rec.update(n_task_blocks=len(task),
                   **{f"n_{c}": int(counts.get(c, 0)) for c in TASKS},
                   task_dur_min=task.duration.min(), task_dur_max=task.duration.max(),
                   instr_dur_min=ev[ev.trial_type == "instruction"].duration.min(),
                   instr_dur_max=ev[ev.trial_type == "instruction"].duration.max(),
                   last_event_end_s=last_end, end_margin_s=scan_end - last_end,
                   s_lost_past_scan_end=lost,
                   flag=";".join(flags) if flags else "ok")
        qc.append(rec)
        fn = os.path.join(EVDIR, f"{r.subject}_task-localizer_run-{r.run}_events.tsv")
        ev.round(4).to_csv(fn, sep="\t", index=False, na_rep="n/a")
    qc = pd.DataFrame(qc)
    # Log-to-scan identity. The stimulus PC clock runs 2-3.5 min behind the
    # scanner's, so the absolute difference is uninformative; what identifies
    # the run is that it agrees with the participant's other localizer run.
    med = qc.groupby("subject").trigger_minus_acqtime_s.transform("median")
    nrun = qc.groupby("subject").run.transform("size")
    qc["clock_dev_from_subject_s"] = qc.trigger_minus_acqtime_s - med
    qc["log_scan_identity"] = np.where(nrun < 2, "single_run_unchecked",
                                np.where(qc.clock_dev_from_subject_s.abs() <= 15, "consistent", "MISMATCH"))
    qc["recommended_use"] = np.select(
        [qc.n_task_blocks < 10, qc.s_lost_past_scan_end > 0],
        ["exclude_partial_log", "model_clip_to_scan_end"], "model")
    # Log filenames carry the session date, and date + subject ID is an indirect
    # identifier (DEIDENT_REVIEW.md), so the committed table names each log file
    # only by its launch order within the participant; the dated names stay in
    # the cache outside the repo.
    qc.to_csv(os.path.join(CACHE, "events_qc_full.csv"), index=False)
    order = {}
    for f in sorted(os.listdir(LOGDIR)):
        m = re.search(r"\.P(\d+)\.Session(\w+?)_(\d{4}-\d{2}-\d{2}_\d+h[\d.]+)\.csv$", f)
        if m:
            order.setdefault(int(m.group(1)), []).append((m.group(3), m.group(2), f))
    lab = {}
    for p_, fs in order.items():
        for k, (_, sess, f) in enumerate(sorted(fs), start=1):
            lab[f] = f"P{p_}_{sess}_log{k}"
    pub = qc.copy()
    pub["selected_file"] = pub.selected_file.map(lambda f: lab.get(f, f))
    pub["other_attempts"] = pub.other_attempts.fillna("").map(
        lambda v: "; ".join(lab.get(x.split("#")[0], x.split("#")[0]) + "#" + x.split("#", 1)[1]
                            for x in v.split("; ") if x))
    pub = pub.rename(columns={"selected_file": "selected_log"})
    pub.to_csv(os.path.join(OUT, "events_qc.csv"), index=False)
    pd.DataFrame(all_blocks).to_csv(os.path.join(CACHE, "blocks.csv"), index=False)
    if PRESSES:
        pd.concat(PRESSES).to_csv(os.path.join(CACHE, "presses.csv"), index=False)
    print(f"[events] {len(qc)} scanned runs; flags:\n{qc.flag.value_counts().to_string()}")
    print(qc.log_scan_identity.value_counts().to_string()); print(qc.recommended_use.value_counts().to_string())
    return qc


# ---------------------------------------------------------------------------
# Motion gate (prereg §4, fixed 2026-08-14): a run fails if mean FD > 0.5 mm
# or > 25 % of its volumes have FD > 0.5 mm. Row-1 FD (NaN) is set to 0, as in
# Scripts/intero2023_covariates_linux.R.
# ---------------------------------------------------------------------------
FD_SPIKE, FD_GATE_FRAC = 0.5, 0.25


def stage_motion():
    ev = pd.read_csv(os.path.join(OUT, "events_qc.csv"))
    rows = []
    for _, r in ev.iterrows():
        f = f"{FMRIPREP}/{r.subject}/func/{r.subject}_task-localizer_run-{r.run}_desc-confounds_timeseries.tsv"
        if not os.path.exists(f):
            rows.append(dict(subject=r.subject, run=r.run, confounds_found=False)); continue
        fd = pd.read_csv(f, sep="\t", usecols=["framewise_displacement"]).framewise_displacement.fillna(0).values
        m, fr = fd.mean(), (fd > FD_SPIKE).mean()
        rows.append(dict(subject=r.subject, run=r.run, confounds_found=True, n_vol_confounds=len(fd),
                         n_vol_matches_header=len(fd) == r.n_vol, fd_mean=round(m, 4),
                         fd_frac_gt_0p5=round(fr, 4), n_spikes=int((fd > FD_SPIKE).sum()),
                         motion_gate_fail=bool(m > FD_SPIKE or fr > FD_GATE_FRAC)))
    mg = pd.DataFrame(rows)
    mg.to_csv(os.path.join(OUT, "motion_gate.csv"), index=False)
    x = mg.merge(ev[["subject", "run", "recommended_use"]], on=["subject", "run"])
    x["usable"] = (~x.motion_gate_fail) & (x.recommended_use != "exclude_partial_log")
    per = x.groupby("subject").usable.sum()
    print(f"[motion] {len(mg)} runs; gate fails {int(mg.motion_gate_fail.sum())}: "
          f"{', '.join(mg[mg.motion_gate_fail].subject + ' run-' + mg[mg.motion_gate_fail].run.astype(str))}")
    print(f"[motion] usable runs (gate + complete log) {int(x.usable.sum())}; participants >=1 run "
          f"{int((per >= 1).sum())}, 2 runs {int((per == 2).sum())}")
    return mg


def usable_runs():
    ev = pd.read_csv(os.path.join(CACHE, "events_qc_full.csv"))
    mg = pd.read_csv(os.path.join(OUT, "motion_gate.csv"))
    x = ev.merge(mg, on=["subject", "run"])
    x["usable"] = (~x.motion_gate_fail) & (x.recommended_use != "exclude_partial_log")
    return x


# ---------------------------------------------------------------------------
# High-pass check: design matrices only.
# ---------------------------------------------------------------------------
CONTRASTS = {"C1_intero_minus_extero": {"intero": 1, "extero": -1},
             "C2_intero_minus_self": {"intero": 1, "self": -1},
             "C4_match_minus_extero": {"match": 1, "extero": -1},
             "C5_match_minus_intero": {"match": 1, "intero": -1},
             "C6_self_minus_rest": {"self": 1, "rest": -1}}
DRIFTS = {"none": None, "cos128": 128.0, "cos200": 200.0}
SWEEP = [128.0, 160.0, 200.0, 256.0, 320.0, 400.0]
SWEEP_ROWS = []


def loc_design(ev, n_vol, hp, hrf="spm"):
    from nilearn.glm.first_level import make_first_level_design_matrix
    ft = TR * (np.arange(n_vol) + 0.5)            # as Scripts/intero2023_twomodel_qc.py
    e = ev[["onset", "duration", "trial_type"]].copy()
    over = (e.onset + e.duration) > n_vol * TR     # clip to the scan, as the breathing designs do
    e.loc[over, "duration"] = (n_vol * TR - e.loc[over, "onset"]).clip(lower=0)
    e = e[e.duration > 0]
    kw = dict(drift_model=None) if hp is None else dict(drift_model="cosine", high_pass=1.0 / hp)
    return make_first_level_design_matrix(ft, events=e, hrf_model=hrf, **kw)


def efficiency(dm, spec):
    if not set(spec) <= set(dm.columns):
        return np.nan
    X = dm.values.astype(float)
    c = np.zeros(X.shape[1])
    for k, w in spec.items():
        c[list(dm.columns).index(k)] = w
    return float(1.0 / (c @ np.linalg.pinv(X.T @ X) @ c))


def vifs(dm, cols):
    X = dm.values.astype(float); names = list(dm.columns); out = {}
    for k in cols:
        if k not in names:
            out[k] = np.nan; continue
        j = names.index(k); y = X[:, j]; Z = np.delete(X, j, axis=1)
        b, *_ = np.linalg.lstsq(Z, y, rcond=None)
        res = y - Z @ b; r2 = 1 - res.var() / y.var()
        out[k] = 1.0 / max(1e-12, 1 - r2)
    return out


def stage_hpf():
    x = usable_runs()
    rows, vrows = [], []
    for _, r in x.iterrows():
        ev = pd.read_csv(f"{EVDIR}/{r.subject}_task-localizer_run-{r.run}_events.tsv", sep="\t")
        for hrf in ("spm", "spm + derivative"):
            eff = {}
            for dname, hp in DRIFTS.items():
                dm = loc_design(ev, int(r.n_vol), hp, hrf)
                for cname, spec in CONTRASTS.items():
                    eff[(dname, cname)] = efficiency(dm, spec)
                if hrf == "spm":
                    v = vifs(dm, TASKS + ["instruction"])
                    vrows.append(dict(subject=r.subject, run=r.run, drift=dname, usable=r.usable,
                                      n_drift_cols=sum(c.startswith("drift") for c in dm.columns),
                                      **{f"vif_{k}": v[k] for k in v}))
            for cname in CONTRASTS:
                e0 = eff[("none", cname)]
                rows.append(dict(subject=r.subject, run=r.run, protocol_version=r.protocol_version,
                                 usable=r.usable, hrf=hrf, contrast=cname,
                                 eff_none=e0, eff_cos128=eff[("cos128", cname)], eff_cos200=eff[("cos200", cname)],
                                 loss_pct_128=100 * (1 - eff[("cos128", cname)] / e0),
                                 loss_pct_200=100 * (1 - eff[("cos200", cname)] / e0)))
        # implicit-rest variant: Rest not modelled, so the 5 task-ish regressors
        # no longer tile the run with the intercept (explains the VIF ~9-12 above)
        for dname, hp in DRIFTS.items():
            dm = loc_design(ev[ev.trial_type != "rest"], int(r.n_vol), hp, "spm")
            v = vifs(dm, [t for t in TASKS if t != "rest"] + ["instruction"])
            vrows.append(dict(subject=r.subject, run=r.run, drift=dname + "_implicit_rest", usable=r.usable,
                              n_drift_cols=sum(c.startswith("drift") for c in dm.columns),
                              **{f"vif_{k}": v[k] for k in v}))
        # cutoff sweep (SPM HRF): cosine sets and low-order polynomials
        dm0 = loc_design(ev, int(r.n_vol), None, "spm")
        e0 = {c: efficiency(dm0, sp) for c, sp in CONTRASTS.items()}
        for hp in SWEEP:
            dm = loc_design(ev, int(r.n_vol), hp, "spm")
            for c, sp in CONTRASTS.items():
                SWEEP_ROWS.append(dict(subject=r.subject, run=r.run, usable=r.usable, drift=f"cos{int(hp)}",
                                       contrast=c, loss_pct=100 * (1 - efficiency(dm, sp) / e0[c])))
        for order in (1, 2, 3):
            from nilearn.glm.first_level import make_first_level_design_matrix
            ft = TR * (np.arange(int(r.n_vol)) + 0.5)
            e = ev[["onset", "duration", "trial_type"]].copy()
            over = (e.onset + e.duration) > r.n_vol * TR
            e.loc[over, "duration"] = (r.n_vol * TR - e.loc[over, "onset"]).clip(lower=0)
            dm = make_first_level_design_matrix(ft, events=e[e.duration > 0], hrf_model="spm",
                                                drift_model="polynomial", drift_order=order)
            for c, sp in CONTRASTS.items():
                SWEEP_ROWS.append(dict(subject=r.subject, run=r.run, usable=r.usable, drift=f"poly{order}",
                                       contrast=c, loss_pct=100 * (1 - efficiency(dm, sp) / e0[c])))
    h = pd.DataFrame(rows); v = pd.DataFrame(vrows)
    sw = pd.DataFrame(SWEEP_ROWS)
    sws = sw[sw.usable].groupby(["drift", "contrast"]).loss_pct.agg(["median", "min", "max"]).round(1).reset_index()
    sws.to_csv(os.path.join(OUT, "hpf_sweep.csv"), index=False)
    print("[hpf] sweep (median % loss):\n" + sws.pivot(index="drift", columns="contrast", values="median").to_string())
    h.to_csv(os.path.join(OUT, "hpf_efficiency.csv"), index=False)
    v.to_csv(os.path.join(OUT, "hpf_vif.csv"), index=False)
    u = h[h.usable & (h.hrf == "spm")]
    # joint two-run model: task columns shared, so participant efficiency = sum over runs
    j = u.groupby(["subject", "contrast"])[["eff_none", "eff_cos128", "eff_cos200"]].sum().reset_index()
    j["loss_pct_128"] = 100 * (1 - j.eff_cos128 / j.eff_none); j["loss_pct_200"] = 100 * (1 - j.eff_cos200 / j.eff_none)
    j.to_csv(os.path.join(OUT, "hpf_efficiency_participant.csv"), index=False)
    print("[hpf] participant (joint-model) % loss:\n" +
          j.groupby("contrast")[["loss_pct_128", "loss_pct_200"]].agg(["median", "min", "max"]).round(1).to_string())
    summ = u.groupby("contrast")[["loss_pct_128", "loss_pct_200"]].agg(["median", "min", "max"]).round(1)
    print("[hpf] % efficiency loss vs no drift (usable runs, SPM HRF):\n" + summ.to_string())
    print("[hpf] VIF median (usable runs):\n" +
          v[v.usable].groupby("drift").median(numeric_only=True).filter(like="vif_").round(2).to_string())
    return h


# ---------------------------------------------------------------------------
# Compliance: Intero presses per block.
# ---------------------------------------------------------------------------
def stage_compliance():
    x = usable_runs()
    rows = []
    for _, r in x.iterrows():
        ev = pd.read_csv(f"{EVDIR}/{r.subject}_task-localizer_run-{r.run}_events.tsv", sep="\t")
        for _, e in ev[ev.trial_type == "intero"].iterrows():
            rows.append(dict(subject=r.subject, run=r.run, protocol_version=r.protocol_version,
                             block_index=e.block_index, n_presses=int(e.n_presses),
                             compliant=int(e.n_presses) >= MIN_PRESSES,
                             run_usable=r.usable, motion_gate_fail=r.motion_gate_fail,
                             recommended_use=r.recommended_use))
    c = pd.DataFrame(rows)
    pr = pd.read_csv(os.path.join(CACHE, "presses.csv"))
    ipi = pr.groupby(["subject", "run", "block_index"]).t.apply(
        lambda t: np.std(np.diff(t)) / np.mean(np.diff(t)) if len(t) > 2 else np.nan).rename("ipi_cv")
    alt = pr.groupby(["subject", "run", "block_index"]).key.apply(
        lambda k: np.mean(np.diff(k.values) != 0) if len(k) > 1 else np.nan).rename("alternation")
    c = c.merge(ipi.reset_index(), how="left").merge(alt.reset_index(), how="left")
    # participant-level
    def part(df):
        u = df[df.run_usable]
        nrun = u.run.nunique()
        nb, nc = len(u), int(u.compliant.sum())
        both = nrun == 2
        per_run_ok = all(u[u.run == k].compliant.any() for k in u.run.unique()) if nrun else False
        return pd.Series(dict(n_usable_runs=nrun, n_intero_blocks=nb, n_compliant=nc,
                              meets_3of4=both and nc >= 3,
                              meets_75pct=nb > 0 and nc / nb >= 0.75 and per_run_ok,
                              froi_eligible=both and nc >= 3 and per_run_ok))
    P = c.groupby("subject").apply(part).reset_index()
    c = c.merge(P, on="subject", how="left")
    c.to_csv(os.path.join(OUT, "compliance.csv"), index=False)
    allb = c
    ub = c[c.run_usable]
    print(f"[compliance] all scanned runs: {len(allb)} Intero blocks, median presses {allb.n_presses.median():.0f}, "
          f"<4 presses {int((~allb.compliant).sum())} ({100*(~allb.compliant).mean():.1f}%)")
    print(f"[compliance] usable runs: {len(ub)} blocks, <4 presses {int((~ub.compliant).sum())}; "
          f"participants with a usable run {len(P[P.n_usable_runs>0])}; meets >=3 of 4 (two runs) "
          f"{int(P.meets_3of4.sum())}; meets >=75% (one- or two-run) {int(P.meets_75pct.sum())}; "
          f"fROI-eligible {int(P.froi_eligible.sum())}")
    return c


# ---------------------------------------------------------------------------
# PHYSIO. Raw Biopac belt, segmented exactly as R/intero_physio.R does:
#   triggers detected at the native 2000 Hz (TTL > 2.5 V rising edges), blocks
#   split at gaps > 5 s, 5-pulse blocks = fieldmaps, >= 100 = BOLD; blocks placed
#   on the logs' clock via the header's earliest_marker_created_at (UTC ->
#   America/Toronto naive), matched to the log's first-trigger wallclock, and
#   graded by grade_time_match() (R/intero_qc.R): fail if dt > 15 s,
#   ambiguous if the runner-up is < 30 s worse, exact if dt <= 2 s and
#   pulses == volumes, degraded if |pulses - volumes| > 5, else good.
# Belt grade of record: R/intero_resp.R belt_grade() on the existing
#   Data/resp_calibration.csv filtered_iqr / band_ratio (fail < 1e-3 V IQR;
#   degraded band_ratio < .90; else good). Not re-derived here.
# ---------------------------------------------------------------------------
ACQ_RE = re.compile(r"^Intero_P([0-9]+)((?:[a-c]|_1|_A|_PartA)?)(?: [(][0-9]+[)])?[.]acq$")
FS_OUT = 50.0                      # belt analysis rate after decimation (2000 / 40)
# --- breath detection (fixed before any localizer measure was computed) ---
BAND = (0.05, 1.0)                 # Hz, zero-phase 2nd-order Butterworth band-pass
MIN_EXT_SEP_S = 0.8                # minimum time between a peak and the next trough
PROM_FRAC = 0.15                   # extremum prominence >= 15 % of the run's p95-p5 range
PERIOD_RANGE = (1.5, 20.0)         # s; breaths outside are discarded
# --- press-to-belt matching rule (pre-specified) ---
MATCH_WIN = (-1.0, 2.0)            # press - belt onset, s: a '1' matches an inhale onset
                                   # (trough) and a '2' an exhale onset (peak) inside it
POL_MIN_R = 0.30                   # |mean Match belt-pacer r| needed to set a participant's polarity
ENTRAIN_SKIP_S = 3.0               # first 3 s of each block dropped for rate/phase measures
N_CHANCE = 200


def acq_parts():
    import hashlib
    out, seen = [], set()
    for f in sorted(os.listdir(PHYSDIR)):
        m = ACQ_RE.match(f)
        if not m:
            continue
        path = os.path.join(PHYSDIR, f)
        h = hashlib.md5(open(path, "rb").read()).hexdigest()
        if h in seen:
            continue
        seen.add(h)
        out.append(dict(pid=int(m.group(1)), part=m.group(2) or "-", path=path))
    return pd.DataFrame(out)


def read_acq(path):
    import bioread
    from scipy.signal import resample_poly
    d = bioread.read_file(path)
    names = [c.name for c in d.channels]
    pick = lambda pat: next((i for i, n in enumerate(names) if re.search(pat, n, re.I)), None)
    ib, it = pick("DA100C"), pick("Scanner|Pulse|Trig")
    if ib is None or it is None:
        return None
    fs = float(d.channels[it].samples_per_second)
    trig = np.asarray(d.channels[it].data, float)
    belt = np.asarray(d.channels[ib].data, float)
    if len(trig) < fs * 30:
        return None
    hi = trig > 2.5
    edges = np.flatnonzero(hi & ~np.r_[False, hi[:-1]])
    t0 = getattr(d, "earliest_marker_created_at", None)
    if t0 is not None:
        t0 = pd.Timestamp(t0)
        t0 = (t0.tz_localize("UTC") if t0.tzinfo is None else t0).tz_convert("America/Toronto").tz_localize(None)
    q = int(round(fs / FS_OUT))
    belt_d = resample_poly(belt, 1, q) if len(belt) > 30 * q else belt
    blocks = []
    if len(edges):
        gaps = np.diff(edges) / fs
        brk = np.r_[0, np.flatnonzero(gaps > 5) + 1, len(edges)]
        for i in range(len(brk) - 1):
            ix = edges[brk[i]:brk[i + 1]]
            n = len(ix)
            blocks.append(dict(onset_s=ix[0] / fs, n_pulses=n,
                               kind="fieldmap" if n == 5 else ("bold" if n >= 100 else "fragment"),
                               abs_time=(t0 + pd.Timedelta(seconds=ix[0] / fs)) if t0 is not None else pd.NaT))
    return dict(fs=fs / q, belt=belt_d, blocks=pd.DataFrame(blocks), t0=t0)


def bandpass(x, fs):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(2, [BAND[0], BAND[1]], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x - np.median(x))


def detect_breaths(x, fs):
    """Peaks and troughs of the band-passed belt, alternated. Polarity-free."""
    from scipy.signal import find_peaks
    rng = np.percentile(x, 95) - np.percentile(x, 5)
    kw = dict(prominence=PROM_FRAC * rng, distance=max(1, int(MIN_EXT_SEP_S * fs)))
    pk, _ = find_peaks(x, **kw); tr, _ = find_peaks(-x, **kw)
    ext = sorted([(i, 1) for i in pk] + [(i, -1) for i in tr])
    alt = []
    for i, k in ext:                     # enforce alternation: keep the more extreme of a run
        if alt and alt[-1][1] == k:
            if (k == 1 and x[i] > x[alt[-1][0]]) or (k == -1 and x[i] < x[alt[-1][0]]):
                alt[-1] = (i, k)
        else:
            alt.append((i, k))
    return np.array([a[0] for a in alt], int), np.array([a[1] for a in alt], int)


def breath_table(x, fs, polarity):
    """Breaths as inhale-onset -> exhale-onset -> next inhale-onset.

    polarity +1: belt rises on inhalation, so inhale onset = trough, exhale onset = peak.
    """
    idx, kind = detect_breaths(x, fs)
    kind = kind * polarity                     # now +1 = exhale onset (end of inhale), -1 = inhale onset
    t = idx / fs
    rows = []
    for j in range(len(idx) - 2):
        if kind[j] == -1 and kind[j + 1] == 1 and kind[j + 2] == -1:
            per = t[j + 2] - t[j]
            if PERIOD_RANGE[0] <= per <= PERIOD_RANGE[1]:
                amp = abs(x[idx[j + 1]] - 0.5 * (x[idx[j]] + x[idx[j + 2]]))
                rows.append(dict(t_in=t[j], t_ex=t[j + 1], t_end=t[j + 2], period=per, amp=amp))
    return pd.DataFrame(rows, columns=["t_in", "t_ex", "t_end", "period", "amp"])


def pacer_radius(t, onset, hc):
    """Pacer circle radius trajectory (InteroPilot pulse()): grow 0.8 hc, hold 0.2 hc,
    shrink 0.8 hc, hold 0.2 hc; full cycle 2 hc. Returns radius in [0, 1] and phase."""
    u = np.mod(t - onset, 2 * hc) / hc         # 0..2
    r = np.where(u < 0.8, u / 0.8, np.where(u < 1.0, 1.0, np.where(u < 1.8, 1 - (u - 1.0) / 0.8, 0.0)))
    return r, 2 * np.pi * np.mod(t - onset, 2 * hc) / (2 * hc)


def icc31(a, b):
    X = np.c_[a, b]; n, k = X.shape
    gm = X.mean(); rm = X.mean(1); cm = X.mean(0)
    msr = k * ((rm - gm) ** 2).sum() / (n - 1)
    mse = ((X - rm[:, None] - cm[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse)


def acq_bold_blocks(pp):
    """Every BOLD trigger block in a participant's .acq parts (rows of acq_parts())."""
    bold = []
    for p in pp.itertuples():
        a = read_acq(p.path)
        if a is None or a["blocks"].empty:
            continue
        for b in a["blocks"][a["blocks"].kind == "bold"].itertuples():
            bold.append(dict(part=p.part, a=a, onset_s=b.onset_s, n_pulses=b.n_pulses, abs_time=b.abs_time))
    return bold


def grade_block(bold, tw, n_vol):
    """Pick the BOLD block nearest the log's first-trigger wallclock `tw` and grade it with
    grade_time_match() (R/intero_qc.R). Returns (block, dt_s, margin_s, pulse_delta, grade)."""
    d = np.array([abs((b["abs_time"] - tw).total_seconds()) if pd.notna(b["abs_time"]) else np.inf
                  for b in bold])
    o = np.argsort(d); b = bold[o[0]]
    dt = d[o[0]]; margin = d[o[1]] - dt if len(o) > 1 else np.inf
    pdel = abs(b["n_pulses"] - n_vol)
    if dt > 15: gr = "fail"
    elif margin < 30: gr = "ambiguous"
    elif dt <= 2 and pdel == 0: gr = "exact"
    elif pdel > 5: gr = "degraded"
    else: gr = "good"
    return b, dt, margin, pdel, gr


def cut_segment(b, n_vol):
    """Belt samples from the block's first pulse to n_vol x TR + 5 s."""
    a = b["a"]; fs = a["fs"]
    i0 = int(round(b["onset_s"] * fs)); i1 = min(len(a["belt"]), i0 + int((n_vol * TR + 5) * fs))
    return a["belt"][i0:i1], fs


def stage_physio():
    x, blocks, presses, U, seg, POL, polarity = physio_segments()
    return physio_measures(x, blocks, presses, U, seg, POL, polarity)


def physio_segments():
    """Localizer belt segments, usability table and per-participant belt polarity."""
    x = usable_runs()
    blocks = pd.read_csv(os.path.join(CACHE, "blocks.csv"))
    presses = pd.read_csv(os.path.join(CACHE, "presses.csv"))
    # log first-trigger wallclock of the attempt of record
    trig = {}
    for _, r in x.iterrows():
        d = read_log(os.path.join(LOGDIR, r.selected_file))
        L = d[(d.Run == f"Localizer{r.run}") & (d.Description == "PulseReceived")]
        trig[(r.subject, r.run)] = pd.to_datetime(L.RunStartTime.iloc[int(r.selected_attempt)])
    # belt grade of record
    cal = pd.read_csv(os.path.join(REPO, "Data", "resp_calibration.csv"))
    cal = cal[cal.task == "localizer"].copy()
    cal["belt_grade"] = np.where(cal.filtered_iqr < 1e-3, "fail", np.where(cal.band_ratio < 0.90, "degraded", "good"))
    cal["run"] = cal.bids_entity.str[-1].astype(int)
    cov = pd.read_csv(os.path.join(REPO, "Data", "physio_coverage.csv"))
    cov = cov[cov.bids_entity.str.contains("localizer", na=False)].copy()
    cov["run"] = cov.bids_entity.str[-1].astype(int)
    parts = acq_parts()
    seg = {}                                   # (subject, run) -> dict(belt, fs, grade, ...)
    urows = []
    for sub in sorted(x.subject.unique()):
        pid = int(sub[4:])
        pp = parts[parts.pid == pid]
        bold = acq_bold_blocks(pp)
        for _, r in x[x.subject == sub].iterrows():
            rec = dict(subject=sub, run=r.run, protocol_version=r.protocol_version, n_vol=r.n_vol,
                       events_usable=r.recommended_use != "exclude_partial_log",
                       motion_gate_fail=r.motion_gate_fail, n_acq_parts=len(pp))
            g = cal[(cal.subject == sub) & (cal.run == r.run)]
            rec["belt_grade"] = g.belt_grade.iloc[0] if len(g) else "no_belt_record"
            rec["band_ratio"] = g.band_ratio.iloc[0] if len(g) else np.nan
            c = cov[(cov.subject == sub) & (cov.run == r.run)]
            rec["grade_time_recorded"] = c.grade_time.iloc[0] if len(c) else "no_record"
            tw = trig[(sub, r.run)]
            if bold and pd.notna(tw):
                b, dt, margin, pdel, gr = grade_block(bold, tw, r.n_vol)
                rec.update(acq_part=b["part"], dt_s=round(dt, 3), margin_s=round(margin, 1),
                           n_pulses=b["n_pulses"], pulse_delta=pdel, grade_time=gr)
                if gr != "fail":
                    belt, fs = cut_segment(b, r.n_vol)
                    seg[(sub, r.run)] = dict(belt=belt, fs=fs, grade_time=gr)
            else:
                rec.update(grade_time="no_acq" if not len(pp) else "no_bold_block")
            urows.append(rec)
    U = pd.DataFrame(urows)

    # ---- polarity, from Match blocks only (instructed: inhale while the circle grows) ----
    ok_belt = lambda k: k in seg and U.set_index(["subject", "run"]).loc[k, "belt_grade"] in ("good", "degraded")
    pol_rows = []
    for (sub, run), sg in seg.items():
        if not ok_belt((sub, run)):
            continue
        xf = bandpass(sg["belt"], sg["fs"]); t = np.arange(len(xf)) / sg["fs"]
        for b in blocks[(blocks.subject == sub) & (blocks.run == run) & (blocks.trial_type == "match")].itertuples():
            w = (t >= b.onset + ENTRAIN_SKIP_S) & (t < b.onset + b.duration)
            if w.sum() < 10 * sg["fs"]:
                continue
            best = None
            for lag in np.arange(0, 2.51, 0.25):
                rr, _ = pacer_radius(t[w] - lag, b.onset, b.half_cycle_s)
                cc = np.corrcoef(xf[w], rr)[0, 1]
                if best is None or abs(cc) > abs(best[1]):
                    best = (lag, cc)
            pol_rows.append(dict(subject=sub, run=run, block_index=b.block_index, lag=best[0], r=best[1]))
    POL = pd.DataFrame(pol_rows)
    prun = POL.groupby(["subject", "run"]).r.mean()
    # Belt polarity is NOT constant across the cohort: the Match-block belt-pacer
    # correlation is strongly negative in every session up to sub-068 and strongly
    # positive from sub-069 on (a hardware/wiring change between two consecutive
    # session dates in early 2024). It is therefore set PER PARTICIPANT from that
    # participant's Match blocks (instructed: inhale as the circle grows), which
    # uses neither Intero presses nor any outcome. |mean r| < POL_MIN_R = unresolved.
    psub = POL.groupby("subject").r.mean()
    polarity = {k: (int(np.sign(v)) if abs(v) >= POL_MIN_R else 0) for k, v in psub.items()}
    U = U.merge(prun.rename("match_pacer_r").reset_index(), how="left")
    U["polarity"] = U.subject.map(polarity).fillna(0).astype(int)
    print(f"[physio] polarity from {len(POL)} Match blocks in {len(prun)} runs / {len(psub)} participants: "
          f"+1 {sum(v == 1 for v in polarity.values())}, -1 {sum(v == -1 for v in polarity.values())}, "
          f"unresolved {sum(v == 0 for v in polarity.values())}; median |r| {prun.abs().median():.2f}; "
          f"runs whose sign disagrees with their participant's: "
          f"{int(sum(np.sign(v) != polarity.get(k[0], 0) for k, v in prun.items() if polarity.get(k[0], 0)))}")
    return x, blocks, presses, U, seg, POL, polarity


def physio_measures(x, blocks, presses, U, seg, POL, polarity):
    from scipy.signal import hilbert
    ok_belt = lambda k: k in seg and U.set_index(["subject", "run"]).loc[k, "belt_grade"] in ("good", "degraded")
    # ---- per-block breathing measures ----
    regdir = os.path.join(REPO, "Data", "regressors")
    brows, trows = [], []
    rng = np.random.default_rng(20261127)
    for (sub, run), sg in seg.items():
        if not ok_belt((sub, run)):
            continue
        fs = sg["fs"]; xf = bandpass(sg["belt"], fs); t = np.arange(len(xf)) / fs
        pol = polarity.get(sub, 0)
        bt = breath_table(xf, fs, pol if pol else 1)   # rate/depth are polarity-invariant
        if bt.empty:
            continue
        amp_run = bt.amp.median()
        ph = np.angle(hilbert(xf * (pol if pol else 1)))   # 0 at the belt maximum = end of inhalation
        rf = os.path.join(regdir, f"{sub}_task-localizer_run-{run}_resp.tsv")
        reg = pd.read_csv(rf, sep="\t", na_values="n/a") if os.path.exists(rf) else None
        B = blocks[(blocks.subject == sub) & (blocks.run == run)]
        hc_run = B.half_cycle_s.median()
        for b in B.itertuples():
            on, off = b.onset + ENTRAIN_SKIP_S, b.onset + b.duration
            inb = bt[(bt.t_in >= on) & (bt.t_end <= off)]
            row = dict(subject=sub, run=run, protocol_version=b.protocol_version, block_index=b.block_index,
                       condition=b.trial_type, grade_time=sg["grade_time"], n_breaths=len(inb),
                       rate_brpm=60 / inb.period.median() if len(inb) >= 2 else np.nan,
                       period_cv=inb.period.std() / inb.period.mean() if len(inb) >= 3 else np.nan,
                       depth_rel=inb.amp.median() / amp_run if len(inb) >= 2 else np.nan)
            if reg is not None and len(reg):
                vt = TR * (np.arange(len(reg)) + 0.5)
                wv = (vt >= on) & (vt < off)
                row["rvt_z_mean"] = reg.rvt_z[wv].mean(); row["resp_rate_reg_brpm"] = reg.resp_rate_brpm[wv].mean()
            # entrainment to the pacer (Extero/Match) or to a virtual pacer of the run's rate (others: null)
            hc = b.half_cycle_s if np.isfinite(b.half_cycle_s) else hc_run
            if np.isfinite(hc):
                w = (t >= on) & (t < off)
                rr, php = pacer_radius(t[w], b.onset, hc)
                z = np.exp(1j * (ph[w] - php))
                row.update(pacer_rate_brpm=60 / (2 * hc), plv=float(np.abs(z.mean())),
                           phase_lag_rad=float(np.angle(z.mean())) if pol else np.nan,
                           rate_over_pacer=row["rate_brpm"] / (60 / (2 * hc)),
                           pacer_is_virtual=b.trial_type not in ("extero", "match"))
            brows.append(row)
            # ---- breath tracking (Intero, compliant blocks only) ----
            if b.trial_type == "intero" and b.n_presses >= MIN_PRESSES:
                P = presses[(presses.subject == sub) & (presses.run == run) & (presses.block_index == b.block_index)]
                pt, pkk = P.t.values, P.key.values
                # secondary, window-free and polarity-invariant: phase-locking of presses
                # to the belt cycle. psi = phase at a '2' press, phase - pi at a '1' press;
                # both are 0 when the press falls exactly on its belt onset.
                def press_lock(tt):
                    ii = np.clip((tt * fs).astype(int), 0, len(ph) - 1)
                    psi = ph[ii] - np.where(pkk == 1, np.pi, 0.0)
                    return np.exp(1j * psi).mean()
                zl = press_lock(pt)
                rec = dict(subject=sub, run=run, block_index=b.block_index, polarity=pol,
                           press_plv=float(np.abs(zl)), press_phase_rad=float(np.angle(zl)) if pol else np.nan)
                ch, chl = [], []
                for _ in range(N_CHANCE):
                    sh = rng.uniform(5, b.duration - 5)
                    tt = b.onset + np.mod(pt - b.onset + sh, b.duration)
                    chl.append(np.abs(press_lock(tt)))
                    if pol:
                        o = np.argsort(tt)
                        ch.append(track(bt, tt[o], pkk[o], b.onset, b.onset + b.duration)["n_correct"])
                rec["press_plv_chance"] = float(np.mean(chl))
                if pol:
                    res = track(bt, pt, pkk, b.onset, b.onset + b.duration)
                    resf = track(breath_table(xf, fs, -pol), pt, pkk, b.onset, b.onset + b.duration)
                    rec.update(**res, n_correct_flipped=resf["n_correct"], n_breaths_flipped=resf["n_breaths"],
                               n_correct_chance=float(np.mean(ch)))
                rec["_ctx"] = (bt, pt, pkk, b.onset, b.duration, np.angle(zl), inb.period.median())
                trows.append(rec)
        # POST HOC (added after the pre-specified window was seen to sit at chance):
        # the same pair rule with the window re-centred on the participant's own
        # constant press lead, estimated from the OTHER Intero block of the run
        # (cross-validated, so no block centres its own window).
        mine = [r_ for r_ in trows if r_["subject"] == sub and r_["run"] == run and "_ctx" in r_]
        for i, r_ in enumerate(mine):
            other = [o for j, o in enumerate(mine) if j != i]
            if not pol or not other or not np.isfinite(r_["_ctx"][6]):
                continue
            ob = other[0]["_ctx"]
            if not np.isfinite(ob[6]):
                continue
            lead = ob[5] / (2 * np.pi) * ob[6]          # s; negative = press before the belt onset
            win = (lead - 1.5, lead + 1.5)
            bt_, pt_, pk_, on_, du_ = r_["_ctx"][:5]
            rc = track(bt_, pt_, pk_, on_, on_ + du_, win)
            chc = []
            for _ in range(N_CHANCE):
                sh = rng.uniform(5, du_ - 5)
                tt = on_ + np.mod(pt_ - on_ + sh, du_); o = np.argsort(tt)
                chc.append(track(bt_, tt[o], pk_[o], on_, on_ + du_, win)["n_correct"])
            r_.update(cv_lead_s=lead, n_correct_cv=rc["n_correct"], n_correct_cv_chance=float(np.mean(chc)))
    for r_ in trows:
        r_.pop("_ctx", None)
    BR = pd.DataFrame(brows); TK = pd.DataFrame(trows)
    return U, POL, BR, TK, polarity


def track(bt, pt, pk, on, off, win=MATCH_WIN):
    """Match presses to belt onsets inside `win`. Greedy, earliest press first."""
    bb = bt[(bt.t_in >= on) & (bt.t_ex < off)]
    used = np.zeros(len(pt), bool)
    def match(t_on, key):
        lo, hi = t_on + win[0], t_on + win[1]
        cand = np.flatnonzero(~used & (pk == key) & (pt >= lo) & (pt <= hi))
        if not len(cand):
            return None
        j = cand[0]; used[j] = True
        return j
    n_corr, lags = 0, []
    for b in bb.itertuples():
        j1 = match(b.t_in, 1); j2 = match(b.t_ex, 2)
        if j1 is not None: lags.append(pt[j1] - b.t_in)
        if j2 is not None: lags.append(pt[j2] - b.t_ex)
        if j1 is not None and j2 is not None and pt[j1] < pt[j2]:
            n_corr += 1
    inwin = (pt >= on) & (pt < off)
    return dict(n_breaths=len(bb), n_correct=n_corr, n_presses=int(inwin.sum()),
                n_false=int((inwin & ~used).sum()), lags=";".join(f"{v:.3f}" for v in lags))


def summarise_physio(U, POL, BR, TK, polarity):
    from scipy import stats
    U = U.copy()
    U["belt_ok"] = U.belt_grade.isin(["good", "degraded"])
    U["usable_rate_depth"] = U.belt_ok & U.grade_time.isin(["exact", "good"]) & U.events_usable
    # Breath-level timing needs the acq block's first pulse to BE the first volume:
    # pulses == volumes and the block identified (exact or good; "good" here is a
    # 2-15 s header-clock offset, which does not move the first pulse).
    U["timing_recoverable"] = U.grade_time.isin(["exact", "good"]) & (U.pulse_delta == 0)
    U["usable_breath_timing"] = U.belt_ok & U.timing_recoverable & U.events_usable & (U.polarity != 0)
    U["usable_breath_timing_strict_exact"] = U.usable_breath_timing & (U.grade_time == "exact")
    U.to_csv(os.path.join(OUT, "belt_usability.csv"), index=False)
    POL.to_csv(os.path.join(CACHE, "polarity_match_blocks.csv"), index=False)
    print(f"[physio] belt grade: {U.belt_grade.value_counts().to_dict()}")
    print(f"[physio] timing grade (re-derived): {U.grade_time.value_counts().to_dict()}; "
          f"recorded: {U.grade_time_recorded.value_counts().to_dict()}")
    print(f"[physio] usable for rate/depth: {int(U.usable_rate_depth.sum())} runs "
          f"({U[U.usable_rate_depth].subject.nunique()} participants); breath-level timing: "
          f"{int(U.usable_breath_timing.sum())} runs ({U[U.usable_breath_timing].subject.nunique()} participants); "
          f"strict exact-only {int(U.usable_breath_timing_strict_exact.sum())} runs")

    # ---------------- breathing by condition ----------------
    key = U.set_index(["subject", "run"])
    BR = BR.merge(U[["subject", "run", "usable_rate_depth"]], on=["subject", "run"])
    BRu = BR[BR.usable_rate_depth]
    per_run = BRu.groupby(["subject", "run", "protocol_version", "condition"]).agg(
        n_blocks=("block_index", "size"), n_breaths=("n_breaths", "sum"),
        rate_brpm=("rate_brpm", "mean"), period_cv=("period_cv", "mean"), depth_rel=("depth_rel", "mean"),
        rvt_z_mean=("rvt_z_mean", "mean"), pacer_rate_brpm=("pacer_rate_brpm", "mean"),
        rate_over_pacer=("rate_over_pacer", "mean"), plv=("plv", "mean")).reset_index()
    per_run.insert(0, "level", "run")
    # participant means (over runs) then group summary per version x condition
    per_sub = per_run.groupby(["subject", "protocol_version", "condition"]).mean(numeric_only=True).reset_index()
    g = []
    for (ver, cond), d in list(per_sub.groupby(["protocol_version", "condition"])) + \
                          [(("all", c), d) for c, d in per_sub.groupby("condition")]:
        rec = dict(level="group", protocol_version=ver, condition=cond, n_participants=d.subject.nunique())
        for m in ["rate_brpm", "period_cv", "depth_rel", "rvt_z_mean", "rate_over_pacer", "plv"]:
            rec[f"{m}_median"] = d[m].median(); rec[f"{m}_iqr"] = d[m].quantile(.75) - d[m].quantile(.25)
        g.append(rec)
    G = pd.DataFrame(g)
    # paired contrasts within participant: Intero vs Extero (and others) per version
    pc = []
    W = per_sub.pivot_table(index=["subject", "protocol_version"], columns="condition",
                            values=["rate_brpm", "depth_rel", "rvt_z_mean", "plv", "period_cv"])
    for ver in ["v1", "v2", "all"]:
        Wv = W if ver == "all" else W.xs(ver, level="protocol_version", drop_level=False)
        for m in ["rate_brpm", "depth_rel", "rvt_z_mean", "plv", "period_cv"]:
            for a, b in [("intero", "extero"), ("intero", "rest"), ("match", "intero"), ("match", "extero"),
                         ("extero", "rest"), ("self", "rest")]:
                if (m, a) not in Wv or (m, b) not in Wv:
                    continue
                dd = (Wv[(m, a)] - Wv[(m, b)]).dropna()
                if len(dd) < 3:
                    continue
                w = stats.wilcoxon(dd) if len(dd) >= 6 else None
                pc.append(dict(level="paired", protocol_version=ver, condition=f"{a}-{b}", measure=m,
                               n_participants=len(dd), diff_median=dd.median(), diff_mean=dd.mean(),
                               dz=dd.mean() / dd.std() if dd.std() > 0 else np.nan,
                               wilcoxon_p=w.pvalue if w is not None else np.nan))
    PC = pd.DataFrame(pc)
    # entrainment: PLV real vs the virtual-pacer null in the same participant
    ent = []
    for ver in ["v1", "v2", "all"]:
        d = per_sub if ver == "all" else per_sub[per_sub.protocol_version == ver]
        null = d[d.condition.isin(["rest", "self", "intero"])].groupby("subject").plv.mean()
        for c in ["extero", "match"]:
            real = d[d.condition == c].set_index("subject").plv
            j = pd.concat([real, null], axis=1, keys=["real", "null"]).dropna()
            ro = d[d.condition == c].rate_over_pacer.dropna()
            ent.append(dict(level="entrainment", protocol_version=ver, condition=c, n_participants=len(j),
                            plv_median=j.real.median(), plv_null_median=j.null.median(),
                            plv_minus_null_median=(j.real - j.null).median(),
                            plv_gt_null_wilcoxon_p=stats.wilcoxon(j.real - j.null).pvalue if len(j) >= 6 else np.nan,
                            n_plv_gt_null=int((j.real > j.null).sum()),
                            rate_over_pacer_median=ro.median(),
                            n_within_10pct_of_pacer=int((ro.sub(1).abs() <= 0.10).sum()), n_rate=len(ro)))
    EN = pd.DataFrame(ent)
    out = pd.concat([per_run, G, PC, EN], ignore_index=True, sort=False)
    out.to_csv(os.path.join(OUT, "breathing_by_condition.csv"), index=False)
    pd.set_option("display.width", 250)
    print("[breathing] group medians (participant means):\n" +
          G[["protocol_version", "condition", "n_participants", "rate_brpm_median", "depth_rel_median",
             "rvt_z_mean_median", "rate_over_pacer_median", "plv_median"]].round(3).to_string(index=False))
    print("[breathing] paired:\n" + PC[PC.condition.isin(["intero-extero", "match-intero", "intero-rest"])]
          .round(4).to_string(index=False))
    print("[breathing] entrainment:\n" + EN.round(3).to_string(index=False))

    # ---------------- tracking accuracy ----------------
    TK = TK.merge(U[["subject", "run", "usable_breath_timing", "protocol_version"]], on=["subject", "run"])
    for c in ["n_breaths", "n_correct", "n_presses", "n_false", "lags", "n_correct_flipped",
              "n_breaths_flipped", "n_correct_chance", "n_correct_cv", "n_correct_cv_chance", "cv_lead_s"]:
        if c not in TK:
            TK[c] = np.nan
    TK = TK.merge(U[["subject", "run", "grade_time"]], on=["subject", "run"])
    runs = []
    for (sub, run), d in TK.groupby(["subject", "run"]):
        lags = np.array([float(v) for s_ in d.lags.dropna() for v in str(s_).split(";") if v != ""])
        nb, nc = d.n_breaths.sum(), d.n_correct.sum()
        zz = (d.press_plv * np.exp(1j * d.press_phase_rad.fillna(0))).mean()
        runs.append(dict(subject=sub, run=run, protocol_version=d.protocol_version.iloc[0],
                         polarity=int(d.polarity.iloc[0]), grade_time=d.grade_time.iloc[0],
                         usable_breath_timing=bool(d.usable_breath_timing.iloc[0]),
                         press_plv=d.press_plv.mean(), press_plv_chance=d.press_plv_chance.mean(),
                         press_phase_rad=float(np.angle(zz)) if d.polarity.iloc[0] else np.nan,
                         n_blocks=len(d), n_belt_breaths=nb, n_presses=d.n_presses.sum(),
                         prop_correct=nc / nb if nb else np.nan,
                         prop_correct_cv_posthoc=(d.n_correct_cv.sum() / d.n_breaths[d.n_correct_cv.notna()].sum()
                                                  if d.n_correct_cv.notna().any() else np.nan),
                         prop_correct_cv_chance=(d.n_correct_cv_chance.sum() / d.n_breaths[d.n_correct_cv.notna()].sum()
                                                 if d.n_correct_cv.notna().any() else np.nan),
                         press_lead_s=d.cv_lead_s.mean(),
                         prop_correct_chance=d.n_correct_chance.sum() / nb if nb else np.nan,
                         prop_correct_flipped_polarity=d.n_correct_flipped.sum() / d.n_breaths_flipped.sum()
                         if d.n_breaths_flipped.sum() else np.nan,
                         lag_median_s=np.median(lags) if len(lags) else np.nan,
                         lag_sd_s=np.std(lags, ddof=1) if len(lags) > 1 else np.nan,
                         false_press_prop=d.n_false.sum() / d.n_presses.sum() if d.n_presses.sum() else np.nan,
                         false_press_per_min=d.n_false.sum() / (len(d) * 0.5)))
    R = pd.DataFrame(runs)
    R.to_csv(os.path.join(OUT, "tracking_accuracy_runs.csv"), index=False)
    Ru = R[R.usable_breath_timing & (R.n_belt_breaths >= 4)]
    rel = []
    for m in ["prop_correct", "lag_median_s", "lag_sd_s", "false_press_prop", "press_plv",
              "press_phase_rad", "prop_correct_cv_posthoc", "press_lead_s"]:
        W = Ru.pivot_table(index="subject", columns="run", values=m).dropna()
        dist = Ru[m]
        rec = dict(measure=m, n_runs=int(dist.notna().sum()), n_participants=Ru.subject.nunique(),
                   median=dist.median(), q25=dist.quantile(.25), q75=dist.quantile(.75),
                   min=dist.min(), max=dist.max(), n_both_runs=len(W))
        if len(W) >= 5:
            rec.update(icc31=icc31(W[1].values, W[2].values), pearson_r=stats.pearsonr(W[1], W[2])[0],
                       spearman_rho=stats.spearmanr(W[1], W[2])[0])
        rel.append(rec)
    REL = pd.DataFrame(rel)
    allpairs = U[U.events_usable]
    miss = dict(measure="missingness",
                n_runs=len(allpairs), n_participants=allpairs.subject.nunique(),
                n_runs_measured=len(Ru), n_participants_measured=Ru.subject.nunique(),
                n_runs_lost_belt=int((~allpairs.belt_ok).sum()),
                n_runs_lost_timing=int((allpairs.belt_ok & ~allpairs.timing_recoverable).sum()),
                n_runs_lost_polarity=int((allpairs.belt_ok & allpairs.timing_recoverable & (allpairs.polarity == 0)).sum()),
                chance_prop_correct_median=Ru.prop_correct_chance.median(),
                chance_press_plv_median=Ru.press_plv_chance.median(),
                chance_prop_correct_cv_median=Ru.prop_correct_cv_chance.median(),
                n_runs_cv_above_chance=int((Ru.prop_correct_cv_posthoc > Ru.prop_correct_cv_chance).sum()),
                n_runs_press_plv_above_chance=int((Ru.press_plv > Ru.press_plv_chance).sum()),
                flipped_polarity_prop_correct_median=Ru.prop_correct_flipped_polarity.median(),
                n_runs_above_chance=int((Ru.prop_correct > Ru.prop_correct_chance).sum()))
    REL = pd.concat([REL, pd.DataFrame([miss])], ignore_index=True, sort=False)
    REL.to_csv(os.path.join(OUT, "tracking_accuracy_reliability.csv"), index=False)
    print("[tracking]\n" + REL.round(3).to_string(index=False))


PRESSES = []


if __name__ == "__main__":
    stages = sys.argv[1:] or ["events", "motion", "hpf", "compliance", "physio"]
    if "events" in stages:
        stage_events()
    if "motion" in stages:
        stage_motion()
    if "hpf" in stages:
        stage_hpf()
    if "compliance" in stages:
        stage_compliance()
    if "physio" in stages:
        import pickle
        pk = os.path.join(CACHE, "physio_stage.pkl")
        if "--reuse" in stages and os.path.exists(pk):
            res = pickle.load(open(pk, "rb"))
        else:
            res = stage_physio(); pickle.dump(res, open(pk, "wb"))
        summarise_physio(*res)
