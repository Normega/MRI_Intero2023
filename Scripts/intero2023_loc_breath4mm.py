#!/usr/bin/env python
"""
intero2023_loc_breath4mm.py -- localizer stage B (OSF p5ztc, plan §3): the breathing-run
contrasts used in transfer, refitted at 4 mm (primary) and unsmoothed (sensitivity),
with no other change to the registered swc7p model.

  python Scripts/intero2023_loc_breath4mm.py <sub>          # 4 mm and unsmoothed fits
  python Scripts/intero2023_loc_breath4mm.py check <sub>    # same code at 6 mm vs the swc7p maps

Everything except `smoothing_fwhm` is intero2023_firstlevel.main / the E3a joint model:
  g.run_list(); g.read_run_physio(sub, ent, n_vol, []) (base nuisance set);
  q.model_a / q.model_b; masker on ROIs/group_brain_mask.nii; mean_scaling per run;
  intero2023_firstlevel.fit (run_glm AR(1) + compute_contrast); ffx over runs;
  hj.build_joint_design for the joint model.
Contrasts (swc7p numbering, see Analysis/localizer_stageB_implementation_notes.md):
  C1  paced - free (Model A per run, FFX)            -> T3, T3b
  C3  confidence modulator (Model B per run, FFX)    -> T2 (exploratory)
  C2joint  aware - unaware (joint Model A)            -> T1
  E3changejoint  change - No Change (joint Model A)   -> T2 (exploratory)
Outputs (outside the repo): $OUT/firstlevel/<sub>/<sub>_sm-<4|0>_var-base_con-<c>_{effect,variance}.nii.gz
and <sub>_breath_sm4_sm0.json. No statistic is summarised here.
"""
import os, sys, json, time, warnings
import numpy as np, nibabel as nib
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
import nilearn
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_glm_common as g
from intero2023_firstlevel import fit, ffx
from intero2023_explore_awareness import CHANGE
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
C2 = {"paced_aware": 1.0, "paced_unaware": -1.0}
KERNELS = {4: 4, 0: None}          # label -> smoothing_fwhm


def fit_subject(sub, fwhm):
    """Returns {contrast: (effect, variance)}, log, masker."""
    man = g.run_list(); runs = man[man.subject == sub]
    assert len(runs), sub
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=fwhm, standardize=False).fit()
    per, Ys, dms, log = {}, [], [], dict(runs=[])
    for r in runs.itertuples():
        n = int(r.n_vol)
        Y = masker.transform(g.bold_path(sub, r.bids_entity))
        assert Y.shape[0] == n
        Y, _ = mean_scaling(Y, axis=0)
        ev, cf, n_clip = g.read_run_physio(sub, r.bids_entity, n, [])
        dmA = q.model_a(ev, cf, n)
        outA, dofA = fit(Y, dmA, {"C1": g.c1_spec(dmA)})
        rec = dict(entity=r.bids_entity, n_vol=n, dof_A=dofA, n_clipped=n_clip)
        if g.conf_varies(ev):
            dmB = q.model_b(ev, cf, n, q.PACED_MOD)
            outB, rec["dof_B"] = fit(Y, dmB, {"C3": {"conf_pm": 1.0}})
            outA.update(outB)
        rec["contrasts"] = sorted(outA)
        for c, ev_ in outA.items():
            per.setdefault(c, []).append(ev_)
        Ys.append(Y); dms.append(dmA); log["runs"].append(rec)
    out = {c: ffx(pairs) for c, pairs in per.items()}
    J = hj.build_joint_design(dms)
    specs = {"C2joint": C2}
    present = [c for c in CHANGE if c in J.columns]
    if "paced_nochange" in J.columns and present:
        sp = {c: 1.0 / len(present) for c in present}; sp["paced_nochange"] = -1.0
        specs["E3changejoint"] = sp
    outJ, dofJ = fit(np.vstack(Ys), J, specs)
    out.update(outJ)
    log["joint"] = dict(dof=dofJ, n_reg=int(J.shape[1]), n_vol=int(sum(y.shape[0] for y in Ys)),
                        contrasts=sorted(outJ))
    log["ffx_runs"] = {c: len(p) for c, p in per.items()}
    return out, log, masker


def main(sub):
    t0 = time.time()
    od = os.path.join(g.OUT, "firstlevel", sub)
    full = dict(subject=sub, nilearn=nilearn.__version__, fits={})
    for lab, fwhm in KERNELS.items():
        out, log, masker = fit_subject(sub, fwhm)
        for c, (e, v) in out.items():
            if not (v > 0).all():
                raise ValueError(f"{sub} sm{lab} {c}: {(v <= 0).sum()} non-positive variances")
            g.save_float32(masker.inverse_transform(e), os.path.join(od, f"{sub}_sm-{lab}_var-base_con-{c}_effect.nii.gz"))
            g.save_float32(masker.inverse_transform(v), os.path.join(od, f"{sub}_sm-{lab}_var-base_con-{c}_variance.nii.gz"))
        full["fits"][f"sm{lab}"] = log
    full["seconds"] = round(time.time() - t0, 1)
    json.dump(full, open(os.path.join(od, f"{sub}_breath_sm4_sm0.json"), "w"), indent=1)
    print(f"{sub} done in {full['seconds']} s", flush=True)


def check(sub):
    """The same code at 6 mm must reproduce the swc7p first-level maps."""
    out, log, masker = fit_subject(sub, 6)
    od = os.path.join(g.OUT, "firstlevel", sub)
    res = dict(subject=sub, fwhm=6, contrasts={})
    raw = NiftiMasker(mask_img=g.GROUP_MASK).fit()          # no smoothing when reading saved maps
    js = json.load(open(os.path.join(od, f"{sub}_firstlevel.json")))
    for c, (e, v) in out.items():
        e0 = raw.transform(nib.load(os.path.join(od, f"{sub}_var-base_con-{c}_effect.nii.gz"))).ravel()
        v0 = raw.transform(nib.load(os.path.join(od, f"{sub}_var-base_con-{c}_variance.nii.gz"))).ravel()
        res["contrasts"][c] = dict(max_abs_diff_effect=float(np.abs(e - e0).max()),
                                   max_abs_diff_variance=float(np.abs(v - v0).max()),
                                   max_abs_effect=float(np.abs(e0).max()))
    res["joint_dof"] = log["joint"]["dof"]; res["swc7p_joint_dof"] = js["joint"]["base"]["dof"]
    path = os.path.join(g.REPO, "Analysis", "localizer_stageB_refit_check.json")
    json.dump(res, open(path, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check(sys.argv[2])
    else:
        main(sys.argv[1])
