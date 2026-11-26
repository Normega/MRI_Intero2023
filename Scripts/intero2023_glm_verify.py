#!/usr/bin/env python
"""
intero2023_glm_verify.py -- handoff §3 verification before any model is fitted.

Reads image HEADERS only (shape, affine, number of volumes); no voxel intensity.
  1. every run: MNI preproc BOLD + brain mask + fMRIPrep confounds exist; volume
     count == manifest n_vol == sheet rows == fMRIPrep confound rows; grid ==
     group mask grid;
  2. events: end <= scan + 10 s except the three registered clipped runs;
     clipped-event count; all seven condition labels present;
  3. ROI masks: SHA-256 of every file matches ROIs/roi_manifest.csv;
  4. design QC: 0 rank-deficient designs in Models A/B/C (Analysis/design_qc_twomodel.csv);
  5. registered per-hypothesis Ns.
Writes Analysis/glm_verification.csv and prints a summary. Exits non-zero on any failure.
"""
import sys
import numpy as np, pandas as pd, nibabel as nib
import intero2023_glm_common as g
q = g.q

fail = []
man = g.run_list()
gm = nib.load(g.GROUP_MASK)
rows = []
CLIPPED = {("sub-052", "task-breathing_run-1"), ("sub-053", "task-breathing_run-1"),
           ("sub-053", "task-breathing_run-2")}
labels = set()
for _, r in man.iterrows():
    s, e, n = r.subject, r.bids_entity, int(r.n_vol)
    rec = dict(subject=s, entity=e, n_vol_manifest=n)
    for k, p in (("bold", g.bold_path(s, e)), ("mask", g.mask_path(s, e)),
                 ("conf", g.fmriprep_confounds_path(s, e))):
        try:
            open(p, "rb").close(); rec[f"has_{k}"] = True
        except OSError:
            rec[f"has_{k}"] = False; fail.append(f"{s} {e}: missing {k}")
    if rec["has_bold"]:
        img = nib.load(g.bold_path(s, e))
        rec["n_vol_bold"] = img.shape[3]
        rec["grid_ok"] = img.shape[:3] == gm.shape[:3] and np.allclose(img.affine, gm.affine, atol=1e-4)
        rec["tr_header"] = float(img.header.get_zooms()[3])
        if rec["n_vol_bold"] != n: fail.append(f"{s} {e}: BOLD {rec['n_vol_bold']} vols vs manifest {n}")
        if not rec["grid_ok"]: fail.append(f"{s} {e}: grid differs from group mask")
        if abs(rec["tr_header"] - g.TR) > 1e-3: fail.append(f"{s} {e}: TR {rec['tr_header']}")
    sheet = pd.read_csv(f"{q.FL}/{s}_{e}_desc-model_confounds.tsv", sep="\t")
    rec["n_rows_sheet"] = len(sheet)
    if len(sheet) != n: fail.append(f"{s} {e}: sheet {len(sheet)} rows vs {n}")
    if rec["has_conf"]:
        rec["n_rows_fmriprep"] = len(pd.read_csv(g.fmriprep_confounds_path(s, e), sep="\t"))
        if rec["n_rows_fmriprep"] != n: fail.append(f"{s} {e}: fMRIPrep confounds {rec['n_rows_fmriprep']} rows")
    ev = pd.read_csv(f"{q.FL}/{s}_{e}_desc-model_events.tsv", sep="\t")
    end = float((ev.onset + ev.duration).max())
    rec["events_end_s"] = round(end, 2); rec["scan_dur_s"] = n * g.TR
    rec["registered_clip"] = (s, e) in CLIPPED
    if end > n * g.TR + 10 and not rec["registered_clip"]:
        fail.append(f"{s} {e}: events end {end:.1f} s > scan {n*g.TR:.0f} + 10 s (not a registered clip)")
    _, _, n_trunc = q.read_run(s, e, n)
    rec["n_clipped_events"] = n_trunc
    labels |= set(ev.condition)
    rows.append(rec)
d = pd.DataFrame(rows)
missing_labels = set(q.CONDITIONS) - labels
if missing_labels: fail.append(f"condition labels absent: {missing_labels}")

# ROI checksums
roi = pd.read_csv(f"{g.REPO}/ROIs/roi_manifest.csv")
bad = [rr.file for rr in roi.itertuples() if g.sha256(f"{g.REPO}/{rr.file}") != rr.sha256_file]
if bad: fail.append(f"ROI checksum mismatch: {bad}")

# Design QC
qc = pd.read_csv(f"{g.REPO}/Analysis/design_qc_twomodel.csv")
rd = {m: int(qc[f"{m}_rank_deficient"].fillna(0).sum()) for m in ("A", "Ball", "C")}
if any(rd.values()): fail.append(f"rank-deficient designs: {rd}")

# Ns
sets = pd.read_csv(f"{g.REPO}/Data/h3_sets.csv")
maia = pd.read_csv(f"{g.REPO}/Data/maia_scores.csv")
subs = sorted(man.subject.unique())
n_h3 = int(sets.in_H3.sum()); n_s7 = int(sets.in_S7.sum())
n_h6 = int(maia[maia.subject.isin(subs)].maia_ieat_total.notna().sum())
print(f"runs {len(d)}, participants {len(subs)}; one-run participants "
      f"{(man.groupby('subject').size() == 1).sum()}")
for s in ("sub-065", "sub-076", "sub-082", "sub-104"):
    print(f"  {s}: {list(man[man.subject == s].bids_entity)}")
print(f"clipped events: {int(d.n_clipped_events.sum())} in runs "
      f"{[(r.subject, r.entity) for r in d[d.n_clipped_events > 0].itertuples()]}")
print(f"ROI masks checked: {len(roi)}, mismatches {len(bad)}")
print(f"rank-deficient designs: {rd}")
print(f"N H3 {n_h3}, S7/H7b {n_s7}, H6 (MAIA) {n_h6}")
d.to_csv(f"{g.REPO}/Analysis/glm_verification.csv", index=False)
if n_h3 != 45 or n_s7 != 23 or n_h6 != 46: fail.append("registered Ns not reproduced")
if fail:
    print("FAIL:"); [print("  ", f) for f in fail]; sys.exit(1)
print("ALL CHECKS PASSED")
