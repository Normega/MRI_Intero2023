"""
intero2023_glm_common.py -- shared paths, run list and design builders for the
registered GLM analysis (prereg v0.8, OSF swc7p as updated 2026-11-26).

Every design is built by the same functions the pre-registration design QC used
(Scripts/intero2023_twomodel_qc.py and Scripts/intero2023_h3_joint_qc.py), so
the fitted designs are the designs that were checked.

Imaging outputs go OUTSIDE the repository (the repo is public; images are not
published): $INTERO_GLM_OUT, default <intero2024>/derivatives/intero2023_glm.
"""
import os, sys, hashlib
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import intero2023_twomodel_qc as q          # noqa: E402  read_run, model_a/b/c, design
import intero2023_h3_joint_qc as hj          # noqa: E402  build_joint_design

REPO = os.path.dirname(HERE)
ROOT = os.path.dirname(REPO)
FMRIPREP = os.environ.get("INTERO_FMRIPREP", os.path.join(ROOT, "derivatives", "fmriprep"))
OUT = os.environ.get("INTERO_GLM_OUT", os.path.join(ROOT, "derivatives", "intero2023_glm"))
TR = 2.0
SPACE = "space-MNI152NLin2009cAsym_res-2"
GROUP_MASK = os.path.join(REPO, "ROIs", "group_brain_mask.nii")
PHYSIO_RVT = ["rvt_rrf_z"]
PHYSIO_CO2 = ["petco2_hrf_lag00", "petco2_hrf_lag10", "petco2_hrf_lag20", "petco2_hrf_lag30"]
ALL_PHYSIO = q.PHYSIO


def bold_path(sub, ent):
    return os.path.join(FMRIPREP, sub, "func", f"{sub}_{ent}_{SPACE}_desc-preproc_bold.nii.gz")


def mask_path(sub, ent):
    return os.path.join(FMRIPREP, sub, "func", f"{sub}_{ent}_{SPACE}_desc-brain_mask.nii.gz")


def fmriprep_confounds_path(sub, ent):
    return os.path.join(FMRIPREP, sub, "func", f"{sub}_{ent}_desc-confounds_timeseries.tsv")


def run_list():
    """The registered run list: both manifests, motion gate passed, sub-056 out (prereg §4)."""
    os.chdir(REPO)
    man = pd.concat([pd.read_csv(f"{q.FL}/firstlevel_manifest.csv"),
                     pd.read_csv(f"{q.FL}/firstlevel_manifest_candidates.csv")], ignore_index=True)
    man = man[(man.motion_gate_fail.astype(str).str.upper() == "FALSE") & (man.subject != "sub-056")]
    man = man.drop_duplicates(["subject", "bids_entity"]).sort_values(["subject", "bids_entity"])
    assert len(man) == 90 and man.subject.nunique() == 47, (len(man), man.subject.nunique())
    return man.reset_index(drop=True)


def read_run_physio(sub, ent, n_vol, physio=()):
    """As q.read_run, but keeps the requested physio columns (sensitivity models)."""
    ev = pd.read_csv(f"{q.FL}/{sub}_{ent}_desc-model_events.tsv", sep="\t")
    cf = pd.read_csv(f"{q.FL}/{sub}_{ent}_desc-model_confounds.tsv", sep="\t")
    assert len(cf) == n_vol
    keep_physio = [c for c in physio if c in cf.columns]
    missing = [c for c in physio if c not in cf.columns]
    if missing:
        raise KeyError(f"{sub} {ent}: physio columns missing {missing}")
    cf = cf.drop(columns=[c for c in ALL_PHYSIO if c in cf.columns and c not in keep_physio])
    ev = ev[ev.condition.isin(q.CONDITIONS)].copy()
    scan_end = n_vol * TR
    over = (ev.onset + ev.duration) > scan_end
    ev.loc[over, "duration"] = (scan_end - ev.loc[over, "onset"]).clip(lower=0)
    ev = ev[ev.duration > 0]
    return ev, cf, int(over.sum())


def c1_spec(dm):
    """C1 (H1): unweighted mean of the paced conditions present in the run, minus free."""
    present = [k for k in q.PACED_SIDE if k in dm.columns]
    spec = {k: 1.0 / len(present) for k in present}
    spec["free"] = -1.0
    return spec


def spec_vector(dm, spec):
    return np.array([spec.get(c, 0.0) for c in dm.columns])


def save_float32(img, path):
    """Write an image as float32. nilearn's inverse_transform copies the header of the
    (uint8) mask, so a plain to_filename() would quantise the data to 8 bits."""
    import nibabel as nib
    out = nib.Nifti1Image(np.asarray(img.get_fdata(dtype=np.float32)), img.affine)
    out.set_data_dtype(np.float32)
    out.to_filename(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def conf_varies(ev):
    """Model B is estimable in a run only if paced confidence varies (prereg §4)."""
    c = ev[ev.condition.isin(q.PACED_MOD)].confidence
    c = pd.to_numeric(c, errors="coerce").dropna()
    return len(c) >= 2 and c.std() > 0
