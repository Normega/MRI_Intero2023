#!/usr/bin/env python
"""
intero2023_loc_firstlevel.py -- registered localizer first levels (OSF p5ztc, plan §3).

One participant per call:  python Scripts/intero2023_loc_firstlevel.py sub-060

Computation = intero2023_firstlevel.fit (nilearn run_glm, AR(1) + compute_contrast),
on group-mask data smoothed at 6 mm and at 4 mm, PSC per run (mean_scaling).
Designs: intero2023_loc_common.design (SPM canonical HRF, DCT high-pass, swc7p nuisance
sheet without fMRIPrep cosine_*), joined over runs by intero2023_loc_common.build_joint_design.
Choices: Analysis/localizer_glm_implementation_notes.md.

Fits:
  6 mm  joint  base (200 s), hp128, hp320, rvt (+rvt_rrf_z, if on every used run),
               co2 (+PetCO2 basis, if on every used run)
  6 mm  per run  base
  4 mm  joint + per run  base       (stage B)
Saved per fit and contrast (C1 C2 C4 C5 C6 intero extero match): effect, variance, t;
C3 = min(t_C1, t_C2) as t only. Outputs (outside the repo):
  $OUT/firstlevel/<sub>/<sub>_task-localizer_loc-<joint|run1|run2>_sm-<6|4>_var-<v>_con-<c>_<kind>.nii.gz
  $OUT/firstlevel/<sub>/<sub>_task-localizer_firstlevel.json
No statistic is summarised or printed here.
"""
import os, sys, json, time, warnings
import numpy as np
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
import nilearn
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intero2023_loc_common as L
from intero2023_firstlevel import fit
g = L.g
warnings.simplefilter("ignore")


def variants_for(runs):
    v = {"base": ((), 200.0), "hp128": ((), 128.0), "hp320": ((), 320.0)}
    if runs.has_rvt.all():
        v["rvt"] = (tuple(L.PHYSIO_RVT), 200.0)
    if runs.has_petco2.all():
        v["co2"] = (tuple(L.PHYSIO_CO2), 200.0)
    return v


def main(sub):
    t0 = time.time()
    runs = L.run_list()
    runs = runs[runs.subject == sub]
    assert len(runs), sub
    os.makedirs(L.fl_dir(sub), exist_ok=True)
    log = dict(subject=sub, nilearn=nilearn.__version__, runs=[], fits={})
    for r in runs.itertuples():
        _, n_clip = L.read_events(sub, r.run, int(r.n_vol))
        log["runs"].append(dict(run=int(r.run), n_vol=int(r.n_vol), n_events_clipped=n_clip,
                                has_rvt=bool(r.has_rvt), has_petco2=bool(r.has_petco2)))
    V = variants_for(runs)
    for sm in (6, 4):
        masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=sm, standardize=False).fit()

        def save(arr, model, var, con, kind):
            if kind == "variance" and not (arr > 0).all():
                raise ValueError(f"{sub} {model} sm{sm} {var} {con}: {(arr <= 0).sum()} non-positive variances")
            g.save_float32(masker.inverse_transform(arr), L.fl_name(sub, model, sm, var, con, kind))

        Y = {}
        for r in runs.itertuples():
            y = masker.transform(g.bold_path(sub, L.entity(r.run)))     # float32, as FirstLevelModel
            assert y.shape[0] == int(r.n_vol)
            Y[int(r.run)], _ = mean_scaling(y, axis=0)
        for var, (physio, hp) in V.items():
            if sm == 4 and var != "base":
                continue
            dms = L.run_designs(sub, runs, physio, hp)
            models = [("joint", np.vstack([Y[k] for k, _ in dms]), L.build_joint_design([dm for _, dm in dms]))]
            if var == "base":
                models += [(f"run{k}", Y[k], dm) for k, dm in dms]
            for model, Ym, dm in models:
                specs = {c: s for c, s in L.CONTRASTS.items() if L.estimable(dm, s)}
                out, dof = fit(Ym, dm, specs)
                tmaps = {}
                for c, (e, v) in out.items():
                    t = (e / np.sqrt(v)).astype(np.float32)
                    save(e, model, var, c, "effect"); save(v, model, var, c, "variance"); save(t, model, var, c, "t")
                    tmaps[c] = t
                if "C1" in tmaps and "C2" in tmaps:
                    save(np.minimum(tmaps["C1"], tmaps["C2"]), model, var, "C3", "t")
                    out["C3"] = None
                log["fits"][f"{model}:sm{sm}:{var}"] = dict(
                    dof=dof, n_vol=int(Ym.shape[0]), n_reg=int(dm.shape[1]), high_pass_s=hp,
                    physio=list(physio), contrasts=sorted(out),
                    task_columns=[c for c in L.TASK if c in dm.columns])
    log["seconds"] = round(time.time() - t0, 1)
    json.dump(log, open(os.path.join(L.fl_dir(sub), f"{sub}_task-localizer_firstlevel.json"), "w"), indent=1)
    print(f"{sub} done in {log['seconds']} s")


if __name__ == "__main__":
    main(sys.argv[1])
