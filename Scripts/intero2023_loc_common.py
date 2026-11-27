"""
intero2023_loc_common.py -- shared paths, run list, participant sets and design
builders for the registered LOCALIZER analysis (OSF p5ztc; plan of record
Preregistration/intero2023_localizer_prereg.md v1.0, frozen at d29136b).

Analyst choices are in Analysis/localizer_glm_implementation_notes.md.
Everything here is built on the breathing-run pipeline (intero2023_glm_common,
intero2023_twomodel_qc, intero2023_h3_joint_qc), imported, never modified.

Imaging outputs go OUTSIDE the repository: $INTERO_GLM_OUT/firstlevel/<sub>/ with
`task-localizer` in every name.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_glm_common as g          # noqa: E402  paths, save_float32, TR, bold_path
from nilearn.glm.first_level import make_first_level_design_matrix   # noqa: E402

REPO, ROOT, TR = g.REPO, g.ROOT, g.TR
OUT = g.OUT
EVDIR = os.path.join(REPO, "Data", "events")
REGDIR = os.path.join(REPO, "Data", "regressors")
QCDIR = os.path.join(REPO, "Analysis", "localizer_qc")
SHEETS = os.path.join(REPO, "FirstLevel", "localizer")
MANIFEST = os.path.join(SHEETS, "localizer_manifest.csv")
FD_SPIKE, FD_GATE_FRAC, PHYSIO_MIN_COV = 0.5, 0.25, 0.5
MIN_PRESSES = 4

# Model regressors. Rest is the implicit baseline (plan §3), so it is not modelled.
TASK = ["intero", "extero", "match", "self", "instruction", "intero_nc"]
CONTRASTS = {
    "C1": {"intero": 1.0, "extero": -1.0},
    "C2": {"intero": 1.0, "self": -1.0},
    "C4": {"match": 1.0, "extero": -1.0},
    "C5": {"match": 1.0, "intero": -1.0},
    "C6": {"self": 1.0},
    "intero": {"intero": 1.0},
    "extero": {"extero": 1.0},
    "match": {"match": 1.0},
}
HP_PRIMARY = 200.0
PHYSIO_RVT = ["rvt_rrf_z"]
PHYSIO_CO2 = ["petco2_hrf_lag00", "petco2_hrf_lag10", "petco2_hrf_lag20", "petco2_hrf_lag30"]
PHYSIO_ALL = ["rvt_rrf_z", "rvt_z", "resp_rate_c"] + PHYSIO_CO2


def entity(run):
    return f"task-localizer_run-{int(run)}"


def events_path(sub, run):
    return os.path.join(EVDIR, f"{sub}_{entity(run)}_events.tsv")


def sheet_path(sub, run):
    return os.path.join(SHEETS, f"{sub}_{entity(run)}_desc-model_confounds.tsv")


def fl_dir(sub):
    return os.path.join(OUT, "firstlevel", sub)


def fl_name(sub, model, sm, var, con, kind):
    """model: 'joint' or 'run1'/'run2'; sm: 4 or 6; kind: effect | variance | t."""
    return os.path.join(fl_dir(sub), f"{sub}_task-localizer_loc-{model}_sm-{sm}_var-{var}_con-{con}_{kind}.nii.gz")


# --------------------------------------------------------------------------- runs / sets
def run_list():
    """The registered run set (p5ztc, Data inclusion 4-5): 102 runs, 54 participants."""
    m = pd.read_csv(MANIFEST)
    u = m[m.usable].sort_values(["subject", "run"]).reset_index(drop=True)
    assert len(u) == 102 and u.subject.nunique() == 54, (len(u), u.subject.nunique())
    return u


def participant_sets():
    """Participant-level sets (p5ztc Data inclusion 8-13)."""
    u = run_list()
    c = pd.read_csv(os.path.join(QCDIR, "compliance.csv")).groupby("subject").first()
    allsub = sorted(u.subject.unique())
    s = dict(all=allsub)
    s["intero"] = sorted(c[c.meets_75pct.astype(bool)].index.intersection(allsub))
    s["froi"] = sorted(c[c.froi_eligible.astype(bool)].index.intersection(allsub))
    maia = pd.read_csv(os.path.join(REPO, "Data", "maia_scores.csv")).set_index("subject")
    s["L3"] = [x for x in s["intero"] if x in maia.index and pd.notna(maia.loc[x, "maia_ieat_total"])]
    gb = u.groupby("subject")
    s["rvt"] = sorted(gb.has_rvt.all().pipe(lambda x: x[x]).index)          # RVT on every used run
    s["rvt_tworun"] = sorted(set(s["rvt"]) & set(gb.size().pipe(lambda x: x[x == 2]).index))
    s["co2"] = sorted(gb.has_petco2.all().pipe(lambda x: x[x]).index)       # PetCO2 basis on every used run
    s["v1"] = sorted(gb.protocol_version.first().pipe(lambda x: x[x == "v1"]).index)
    s["v2"] = sorted(set(allsub) - set(s["v1"]))
    return s, maia, u


# --------------------------------------------------------------------------- design
def read_events(sub, run, n_vol):
    """Registered event coding: non-compliant Intero blocks (< 4 presses) -> intero_nc;
    Rest dropped (implicit baseline); events clipped to the scan end."""
    ev = pd.read_csv(events_path(sub, run), sep="\t")
    ev = ev[["onset", "duration", "trial_type", "n_presses", "compliant"]].copy()
    it = ev.trial_type == "intero"
    nc = it & (pd.to_numeric(ev.n_presses, errors="coerce") < MIN_PRESSES)
    assert (nc == (it & (pd.to_numeric(ev.compliant, errors="coerce") == 0))).all(), (sub, run)
    ev.loc[nc, "trial_type"] = "intero_nc"
    ev = ev[ev.trial_type != "rest"]
    scan_end = n_vol * TR
    over = (ev.onset + ev.duration) > scan_end
    ev.loc[over, "duration"] = (scan_end - ev.loc[over, "onset"]).clip(lower=0)
    n_clip = int(over.sum())
    ev = ev[ev.duration > 0]
    return ev[["onset", "duration", "trial_type"]].reset_index(drop=True), n_clip


def read_confounds(sub, run, n_vol, physio=()):
    """The swc7p nuisance sheet with fMRIPrep's 128 s cosine_* columns dropped (p5ztc
    Transformations 5) and the requested physio columns kept."""
    cf = pd.read_csv(sheet_path(sub, run), sep="\t")
    assert len(cf) == n_vol, (sub, run, len(cf), n_vol)
    missing = [c for c in physio if c not in cf.columns]
    if missing:
        raise KeyError(f"{sub} run-{run}: physio columns missing {missing}")
    drop = [c for c in cf.columns if c.startswith("cosine")]
    drop += [c for c in PHYSIO_ALL if c in cf.columns and c not in physio]
    return cf.drop(columns=drop)


def design(ev, cf, n_vol, hp=HP_PRIMARY):
    """One run: SPM canonical HRF, boxcars over the exact event durations, DCT high-pass
    at `hp` s (nilearn cosine set, as the step-1 efficiency QC), nuisance as add_regs.
    Frame times TR x (k + 0.5), as the breathing designs."""
    ft = TR * (np.arange(n_vol) + 0.5)
    return make_first_level_design_matrix(
        ft, events=ev, hrf_model="spm", drift_model="cosine", high_pass=1.0 / hp,
        add_regs=cf.values, add_reg_names=list(cf.columns))


def build_joint_design(run_designs, task=TASK):
    """intero2023_h3_joint_qc.build_joint_design with the localizer task columns: task
    regressors shared across runs; nuisance, drift and intercept run-specific."""
    tcols = [c for c in task if any(c in dm.columns for dm in run_designs)]
    parts = []
    for k, dm in enumerate(run_designs):
        nuis = dm.drop(columns=[c for c in dm.columns if c in tcols]).add_prefix(f"run{k + 1}_")
        parts.append((dm.reindex(columns=tcols, fill_value=0.0).reset_index(drop=True),
                      nuis.reset_index(drop=True)))
    ncols = [c for _, n in parts for c in n.columns]
    return pd.concat([pd.concat([t, n.reindex(columns=ncols, fill_value=0.0)], axis=1)
                      for t, n in parts], ignore_index=True)


def spec_vector(dm, spec):
    return np.array([spec.get(c, 0.0) for c in dm.columns])


def estimable(dm, spec):
    """Contrast estimable: all its columns present and c in the row space of X."""
    if not set(spec) <= set(dm.columns):
        return False
    X = dm.values.astype(float)
    sd = X.std(axis=0); sd[sd < 1e-12] = 1.0
    Xs = X / sd; c = spec_vector(dm, spec) / sd
    u, sv, vt = np.linalg.svd(Xs, full_matrices=False)
    rank = int((sv > max(Xs.shape) * np.finfo(float).eps * sv.max()).sum())
    Vr = vt[:rank]
    return bool(np.linalg.norm(c - Vr.T @ (Vr @ c)) <= 1e-8 * max(1.0, np.linalg.norm(c)))


def run_designs(sub, runs, physio=(), hp=HP_PRIMARY):
    """Per-run designs for one participant's used runs (runs: manifest rows)."""
    out = []
    for r in runs.itertuples():
        ev, _ = read_events(sub, r.run, int(r.n_vol))
        cf = read_confounds(sub, r.run, int(r.n_vol), physio)
        out.append((int(r.run), design(ev, cf, int(r.n_vol), hp)))
    return out
