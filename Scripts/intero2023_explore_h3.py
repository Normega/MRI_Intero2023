#!/usr/bin/env python
"""
intero2023_explore_h3.py -- EXPLORATORY (deviations.md, 2026-11-26, E1): aware vs unaware,
two-sided, in the three H3 masks, on the registered H3 FLAME fit (base, N = 45).

Per mask and tail: peak z, GRF voxel FWE p (3-D EC term, as the registered test), Bonferroni
x6 (3 masks x 2 tails); REML meta-analysis of ROI means, two-sided. Nothing is refitted.
Output: Analysis/results/EXPLORATORY_E1_H3_twosided.csv
"""
import os, json
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
import intero2023_glm_common as g
import intero2023_group as G

label = "H3_var-base_set-all"
od = os.path.join(g.OUT, "group", label)
meta = json.load(open(os.path.join(od, "flame_meta.json")))
z = nib.load(os.path.join(od, "stats", "zstat1.nii.gz")).get_fdata()
cope = nib.load(os.path.join(od, "cope.nii.gz")).get_fdata()
var = nib.load(os.path.join(od, "varcope.nii.gz")).get_fdata()
ref = nib.load(g.GROUP_MASK)
rows = []
for roi in G.SETS["H3"]:
    m = nib.load(G.ROI(roi)).get_fdata() > 0
    R = m.sum() / meta["smoothest"]["RESELS"]
    y, v = cope[m].mean(axis=0), var[m].mean(axis=0)
    mt = G.reml_meta(y, v); p2 = 2 * stats.norm.sf(abs(mt["z"]))
    for tail, sgn in (("aware>unaware", 1), ("aware<unaware", -1)):
        zz = np.where(m, sgn * z, -np.inf); ijk = np.unravel_index(np.argmax(zz), z.shape)
        zmax = float(zz[ijk]); p = G.grf_p(zmax, R) if zmax > 1 else 1.0
        rows.append(dict(analysis="EXPLORATORY E1", roi=roi, tail=tail, N=len(meta["subjects"]),
                         peak_z=round(zmax, 3),
                         peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(ref.affine, ijk)],
                         p_fwe_tail=p, p_fwe_bonf6=min(1.0, 6 * p),
                         z_thresh_bonf6=round(float(G.sh(f"ptoz {0.05 / 6} -g {R}").strip()), 3),
                         reml_mu=mt["mu"], reml_z=mt["z"], reml_p_two_sided=p2))
d = pd.DataFrame(rows)
out = os.path.join(g.REPO, "Analysis", "results", "EXPLORATORY_E1_H3_twosided.csv")
d.to_csv(out, index=False)
print(d.drop(columns=["analysis"]).to_string(index=False))
