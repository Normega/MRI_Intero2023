#!/usr/bin/env python
"""
intero2023_firstlevel.py -- registered first-level models (prereg v0.8 §7-§8).

One participant per call:  python Scripts/intero2023_firstlevel.py sub-057

Per run (the same computation as nilearn FirstLevelModel, done once per run so
every model and variant reuses the same smoothed data; checked on sub-082 run-1
to reproduce FirstLevelModel's effect and variance exactly):
  * NiftiMasker on ROIs/group_brain_mask.nii, 6 mm FWHM smoothing;
  * percent-signal-change scaling per voxel (nilearn signal_scaling=0);
  * run_glm with the AR(1) noise model on the design built by the SAME
    functions the pre-registration design QC used (twomodel_qc.model_a/b/c);
  * 'spm + derivative' HRF, drift_model=None (cosines are in the confounds).
Models and contrasts:
  A -> C1 paced_gt_free (paced conditions present in the run, minus free),
       C2 aware_gt_unaware (per run, where both cells exist; robustness form)
  B -> C3 confidence_mod (only runs whose paced confidence varies)
  C -> C5 salience_high_gt_low (where both cells exist)
Variants: base; +rvt (rvt_rrf_z); +co2 (4-lag PetCO2 basis); +both, each only
on runs that have the regressor (manifest has_rvt / has_petco2).
Participant level: precision-weighted fixed effects over runs (inverse-variance,
voxelwise) for C1/C2/C3/C5.
H3 joint model (prereg §8, v0.7): Model A designs of all of the participant's
runs joined by build_joint_design() (shared task columns, run-specific nuisance
and intercept), each run PSC-scaled before concatenation, AR(1); C2 effect,
variance and residual dof. Variants as above when every run has the regressor.

Outputs (outside the repo): $OUT/firstlevel/<sub>/
  <sub>_run-<k>_var-<v>_con-<c>_{effect,variance}.nii.gz  (per run)
  <sub>_var-<v>_con-<c>_{effect,variance}.nii.gz          (FFX over runs)
  <sub>_var-<v>_con-C2joint_{effect,variance}.nii.gz      (H3 joint model)
  <sub>_firstlevel.json                                   (runs, dofs, columns)
No statistic is summarised or printed here.
"""
import os, sys, json, time, warnings
import numpy as np, pandas as pd, nibabel as nib
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level import run_glm
from nilearn.glm.first_level.first_level import mean_scaling
from nilearn.glm.contrasts import compute_contrast
import nilearn
import intero2023_glm_common as g
q, hj = g.q, g.hj
warnings.simplefilter("ignore")

VARIANTS = {"base": [], "rvt": g.PHYSIO_RVT, "co2": g.PHYSIO_CO2, "both": g.PHYSIO_RVT + g.PHYSIO_CO2}


def fit(Y, dm, specs):
    """run_glm + contrasts; returns {name: (effect, variance)}, dof."""
    X = dm.values.astype(float)
    labels, results = run_glm(Y, X, noise_model="ar1")
    dof = int(X.shape[0] - np.linalg.matrix_rank(X))
    out = {}
    for name, spec in specs.items():
        if spec is None or not set(spec) <= set(dm.columns):
            continue
        con = compute_contrast(labels, results, g.spec_vector(dm, spec), stat_type="t")
        out[name] = (np.asarray(con.effect_size()).ravel().astype(np.float32),
                     np.asarray(con.effect_variance()).ravel().astype(np.float32))
    return out, dof


def ffx(pairs):
    w = [1.0 / np.maximum(v, 1e-12) for _, v in pairs]
    W = np.sum(w, axis=0)
    return (np.sum([wi * e for wi, (e, _) in zip(w, pairs)], axis=0) / W).astype(np.float32), \
           (1.0 / W).astype(np.float32)


def main(sub):
    t0 = time.time()
    man = g.run_list()
    runs = man[man.subject == sub]
    assert len(runs), sub
    od = os.path.join(g.OUT, "firstlevel", sub); os.makedirs(od, exist_ok=True)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    save = lambda arr, name: masker.inverse_transform(arr).to_filename(os.path.join(od, name))
    log = dict(subject=sub, nilearn=nilearn.__version__, runs=[], joint={})
    C2 = {"paced_aware": 1.0, "paced_unaware": -1.0}
    C5 = {"paced_highsal": 1.0, "paced_lowsal": -1.0}
    per = {}            # (variant, contrast) -> list of (effect, variance)
    joint_parts = {v: [] for v in VARIANTS}   # variant -> list of (Ypsc, dmA)
    for r in runs.itertuples():
        k = r.bids_entity.split("_")[-1]
        has = dict(base=True, rvt=str(r.has_rvt).upper() == "TRUE",
                   co2=str(r.has_petco2).upper() == "TRUE")
        has["both"] = has["rvt"] and has["co2"]
        Y = masker.transform(g.bold_path(sub, r.bids_entity))   # float32, as FirstLevelModel
        assert Y.shape[0] == int(r.n_vol)
        Y, _ = mean_scaling(Y, axis=0)
        rl = dict(entity=r.bids_entity, n_vol=int(r.n_vol), rvt_source=r.rvt_source, variants={})
        for v, cols in VARIANTS.items():
            if not has[v]:
                continue
            ev, cf, _ = g.read_run_physio(sub, r.bids_entity, int(r.n_vol), cols)
            n = int(r.n_vol)
            dmA = q.model_a(ev, cf, n)
            outA, dofA = fit(Y, dmA, {"C1": g.c1_spec(dmA), "C2": C2})
            joint_parts[v].append((Y, dmA))
            vrec = dict(dof_A=dofA, n_reg_A=dmA.shape[1], contrasts=sorted(outA))
            if g.conf_varies(ev):
                dmB = q.model_b(ev, cf, n, q.PACED_MOD)
                outB, vrec["dof_B"] = fit(Y, dmB, {"C3": {"conf_pm": 1.0}})
                outA.update(outB)
            dmC = q.model_c(ev, cf, n)
            if dmC is not None:
                outC, vrec["dof_C"] = fit(Y, dmC, {"C5": C5})
                outA.update(outC)
            vrec["contrasts"] = sorted(outA)
            for c, (e, var) in outA.items():
                save(e, f"{sub}_{k}_var-{v}_con-{c}_effect.nii.gz")
                save(var, f"{sub}_{k}_var-{v}_con-{c}_variance.nii.gz")
                per.setdefault((v, c), []).append((e, var))
            rl["variants"][v] = vrec
        log["runs"].append(rl)
    # participant-level FFX
    n_runs = len(runs)
    log["ffx"] = {}
    for (v, c), pairs in per.items():
        e, var = ffx(pairs)
        save(e, f"{sub}_var-{v}_con-{c}_effect.nii.gz")
        save(var, f"{sub}_var-{v}_con-{c}_variance.nii.gz")
        log["ffx"][f"{v}:{c}"] = len(pairs)
    # H3 joint model: only when every run has the variant's regressors
    for v, parts in joint_parts.items():
        if len(parts) != n_runs:
            continue
        J = hj.build_joint_design([dm for _, dm in parts])
        Yj = np.vstack([Yr for Yr, _ in parts])
        out, dof = fit(Yj, J, {"C2joint": C2})
        if "C2joint" in out:
            e, var = out["C2joint"]
            save(e, f"{sub}_var-{v}_con-C2joint_effect.nii.gz")
            save(var, f"{sub}_var-{v}_con-C2joint_variance.nii.gz")
            log["joint"][v] = dict(dof=dof, n_reg=J.shape[1], n_vol=int(Yj.shape[0]))
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(os.path.join(od, f"{sub}_firstlevel.json"), "w"), indent=1)
    print(f"{sub} done in {log['seconds']} s")


if __name__ == "__main__":
    main(sys.argv[1])
