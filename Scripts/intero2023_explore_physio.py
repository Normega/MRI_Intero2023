#!/usr/bin/env python
"""
intero2023_explore_physio.py -- EXPLORATORY E2 (deviations.md, 2026-11-26): whole-brain maps
of the physiological regressors themselves, from the registered +RVT / +PetCO2 Model A fits.

Stage 1 (per participant):  python Scripts/intero2023_explore_physio.py subject sub-057
  Model A + rvt_rrf_z  -> contrast RVT   (beta of RVT x RRF, % signal per SD of RVT)
  Model A + 4 CO2 lags -> contrast CO2sum (sum of the 4 lag betas: net PetCO2 response)
                          and CO2lag00/10/20/30 (kept, not tested)
  Same masker, PSC scaling, AR(1) and precision-weighted FFX over runs as the registered first level.
Stage 2:  python Scripts/intero2023_explore_physio.py group
  One-sample sign-flip permutation, 10,000, two-sided, cluster-forming p < .001 per tail,
  cluster FWE: RVT in S2b (41) and S2 (30), CO2sum in S1 (25).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
from nilearn.glm.second_level import non_parametric_inference
import intero2023_glm_common as g
from intero2023_firstlevel import fit, ffx
q = g.q
warnings.simplefilter("ignore")
LAGS = g.PHYSIO_CO2


def subject(sub):
    man = g.run_list(); runs = man[man.subject == sub]
    od = os.path.join(g.OUT, "firstlevel", sub)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    per = {}
    for r in runs.itertuples():
        has_rvt = str(r.has_rvt).upper() == "TRUE"; has_co2 = str(r.has_petco2).upper() == "TRUE"
        if not (has_rvt or has_co2):
            continue
        Y, _ = mean_scaling(masker.transform(g.bold_path(sub, r.bids_entity)), axis=0)
        n = int(r.n_vol)
        if has_rvt:
            ev, cf, _ = g.read_run_physio(sub, r.bids_entity, n, g.PHYSIO_RVT)
            out, _ = fit(Y, q.model_a(ev, cf, n), {"RVT": {"rvt_rrf_z": 1.0}})
            for c, ev_ in out.items(): per.setdefault(c, []).append(ev_)
        if has_co2:
            ev, cf, _ = g.read_run_physio(sub, r.bids_entity, n, LAGS)
            specs = {"CO2sum": {l: 1.0 for l in LAGS}}
            specs.update({f"CO2{l[-5:]}": {l: 1.0} for l in LAGS})
            out, _ = fit(Y, q.model_a(ev, cf, n), specs)
            for c, ev_ in out.items(): per.setdefault(c, []).append(ev_)
    nruns = len(runs); log = {}
    for c, pairs in per.items():
        e, v = ffx(pairs)
        g.save_float32(masker.inverse_transform(e), os.path.join(od, f"{sub}_var-physio_con-{c}_effect.nii.gz"))
        g.save_float32(masker.inverse_transform(v), os.path.join(od, f"{sub}_var-physio_con-{c}_variance.nii.gz"))
        log[c] = dict(n_runs=len(pairs), all_runs=len(pairs) == nruns)
    json.dump(log, open(os.path.join(od, f"{sub}_explore_physio.json"), "w"), indent=1)
    print(sub, log)


def group():
    import intero2023_group as G
    sets, _, _ = G.subject_sets()
    rd = os.path.join(g.REPO, "Analysis", "results")
    for con, subset in (("RVT", "S2b"), ("RVT", "S2"), ("CO2sum", "S1")):
        subs = sets[subset]
        files = [os.path.join(g.OUT, "firstlevel", s, f"{s}_var-physio_con-{con}_effect.nii.gz") for s in subs]
        files = [f for f in files if os.path.exists(f)]
        assert len(files) == len(subs), (con, subset, len(files), len(subs))
        label = f"EXPLORATORY_E2_{con}_set-{subset}"
        od = os.path.join(g.OUT, "group", label); os.makedirs(od, exist_ok=True)
        dm = pd.DataFrame({"intercept": np.ones(len(files))})
        wb = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept",
                                      mask=g.GROUP_MASK, n_perm=G.N_PERM, two_sided_test=True,
                                      threshold=0.002, random_state=G.SEED, n_jobs=G.N_JOBS)
        for k in ("t", "logp_max_size"):
            g.save_float32(wb[k], os.path.join(od, f"wholebrain_{k}.nii.gz"))
        cl = G.cluster_table(wb["t"], wb["logp_max_size"], len(files) - 1)
        cl.to_csv(os.path.join(rd, f"{label}_wholebrain.csv"), index=False)
        print(label, "N", len(files), "clusters p<.05:", int((cl.p_fwe < 0.05).sum()) if len(cl) else 0)


if __name__ == "__main__":
    subject(sys.argv[2]) if sys.argv[1] == "subject" else group()
