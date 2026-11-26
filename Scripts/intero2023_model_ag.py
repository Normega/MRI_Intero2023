#!/usr/bin/env python
"""
intero2023_model_ag.py -- EXPLORATORY graded-awareness model, Model A-g (prereg v0.8 §10).

One participant per call:  python Scripts/intero2023_model_ag.py sub-057
Joint two-run form exactly as the H3 joint model (shared task columns, run-specific
nuisance and intercept, PSC per run, AR(1)), base nuisance only. Three versions:
  AG      aware + unaware change epochs as ONE regressor (paced_change) with
          modulator g = +confidence (aware) / -confidence (unaware), mean-centred
          within run (paced_gmod); wrong-direction, No Change, free, rating and
          sub-rating as in Model A. Contrast: paced_gmod.
  AGrank  as AG with g rank-transformed within run (average ranks), then centred.
  AGaw    "confidence on aware trials only": Model A (aware and unaware separate)
          plus a confidence modulator on the aware epochs only, centred within the
          run's aware trials (paced_awconf). Contrast: paced_awconf. This reading of
          the §10 variant was fixed before any result was computed.
Outputs $OUT/firstlevel/<sub>/<sub>_var-base_con-<AG|AGrank|AGaw>joint_{effect,variance}.nii.gz
and <sub>_model_ag.json (dofs). Group level: intero2023_group.py AG|AGrank|AGaw (FLAME, H3 masks).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd
from scipy.stats import rankdata
from nilearn.maskers import NiftiMasker
from nilearn.glm.first_level.first_level import mean_scaling
import intero2023_glm_common as g
from intero2023_firstlevel import fit
q, hj = g.q, g.hj
warnings.simplefilter("ignore")


def events_ag(ev, kind):
    ev = ev.copy()
    ev["confidence"] = pd.to_numeric(ev.confidence, errors="coerce")
    base_types = ["paced_wrongdir", "paced_nochange", "free", "rating", "subrating"]
    if kind == "AGaw":
        base_types = ["paced_aware", "paced_unaware"] + base_types
    rows = [ev[ev.condition.isin(base_types)][["onset", "duration", "condition"]]
            .rename(columns={"condition": "trial_type"}).assign(modulation=1.0)]
    if kind in ("AG", "AGrank"):
        ch = ev[ev.condition.isin(["paced_aware", "paced_unaware"]) & ev.confidence.notna()]
        gv = np.where(ch.condition == "paced_aware", 1.0, -1.0) * ch.confidence.values
        if kind == "AGrank":
            gv = rankdata(gv)
        rows.append(ch[["onset", "duration"]].assign(trial_type="paced_change", modulation=1.0))
        if len(ch) >= 2 and np.std(gv) > 0:
            rows.append(ch[["onset", "duration"]].assign(trial_type="paced_gmod", modulation=gv - gv.mean()))
    else:
        aw = ev[(ev.condition == "paced_aware") & ev.confidence.notna()]
        if len(aw) >= 2 and aw.confidence.std() > 0:
            rows.append(aw[["onset", "duration"]].assign(trial_type="paced_awconf",
                                                           modulation=aw.confidence - aw.confidence.mean()))
    return pd.concat(rows, ignore_index=True)


def main(sub):
    man = g.run_list()
    runs = man[man.subject == sub]
    od = os.path.join(g.OUT, "firstlevel", sub)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    Ys, evs = [], []
    for r in runs.itertuples():
        Y = masker.transform(g.bold_path(sub, r.bids_entity))
        Ys.append(mean_scaling(Y, axis=0)[0])
        ev, cf, _ = g.read_run_physio(sub, r.bids_entity, int(r.n_vol), [])
        evs.append((ev, cf, int(r.n_vol)))
    log = {}
    for kind, col in (("AG", "paced_gmod"), ("AGrank", "paced_gmod"), ("AGaw", "paced_awconf")):
        dms = [q.design(events_ag(ev, kind), cf, n) for ev, cf, n in evs]
        J = hj.build_joint_design(dms)
        if col not in J.columns:
            log[kind] = "not estimable"; continue
        out, dof = fit(np.vstack(Ys), J, {kind: {col: 1.0}})
        e, v = out[kind]
        g.save_float32(masker.inverse_transform(e), os.path.join(od, f"{sub}_var-base_con-{kind}joint_effect.nii.gz"))
        g.save_float32(masker.inverse_transform(v), os.path.join(od, f"{sub}_var-base_con-{kind}joint_variance.nii.gz"))
        log[kind] = dict(dof=dof, n_reg=J.shape[1])
    json.dump(log, open(os.path.join(od, f"{sub}_model_ag.json"), "w"), indent=1)
    print(sub, {k: (v if isinstance(v, str) else "ok") for k, v in log.items()})


if __name__ == "__main__":
    main(sys.argv[1])
