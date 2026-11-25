#!/usr/bin/env fslpython
"""Build the registered ROI masks (prereg v0.5, section 8) on the analysis grid.

Run with FSL's Python (nibabel + numpy):
    $FSLDIR/bin/fslpython Scripts/intero2023_build_rois.py

Inputs (no BOLD intensities are read):
  * CAREN atlas (Doucet, Lee & Frangou, 2019, Hum Brain Mapp 40:4577) --
    network masks and network x AAL-3 region files, MNI152 2 mm.
    Default Atlas/caren/ (gitignored: redistribution terms not established).
  * Harvard-Oxford cortical maxprob-thr50 2 mm (FSL) -- atlas-sensitivity set.
  * fMRIPrep `desc-brain_mask` images (space-MNI152NLin2009cAsym_res-2) for the
    90 analysed runs, to fix the target grid and the group coverage mask.

Target grid: fMRIPrep MNI152NLin2009cAsym res-2 (97 x 115 x 97, 2 mm).
Atlas -> grid by sform, nearest neighbour (the procedure used in MindLock,
deviations C6). The MNI152 (FSL/ICBM) vs 2009cAsym template difference is
small and declared in the prereg.

Outputs (ROIs/):
  <set>_<roi>.nii        uint8 tested mask = atlas ROI AND group brain mask
                         (gitignored; checksums are the registered record)
  group_brain_mask.nii   intersection of the 90 run brain masks
  roi_manifest.csv       definition, hypotheses, voxels, coverage, centroid,
                         sha256 of file and of data array (committed)
  roi_overlap.csv        pairwise shared voxels (committed)
  roi_sources.csv        sha256 of every atlas file used (committed)
"""
import csv
import hashlib
import os
import sys

import nibabel as nib
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAREN = os.environ.get("CAREN_DIR", os.path.join(REPO, "Atlas", "caren"))
FMRIPREP = os.environ.get(
    "FMRIPREP_DIR", "/media/normega/FB_Cluster_Data/intero2024/derivatives/fmriprep")
FSLDIR = os.environ.get("FSLDIR", "/home/normega/fsl")
HO = os.path.join(FSLDIR, "data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz")
OUT = os.path.join(REPO, "ROIs")
SPACE = "space-MNI152NLin2009cAsym_res-2"

# CAREN network numbering (verified against CAREN_byRegions_AAL.xlsx in MindLock):
# RSN01 SAL, RSN02 ECN, RSN03 SMN, RSN04 DMN, RSN05 VIS.  Region files are
# CAREN_RSNnn_AAL_RegKK.nii with KK the AAL-3 index (Atlas/caren/CAREN_regions.csv).
def caren_regions(net, idx):
    return [os.path.join(CAREN, net, f"CAREN_{net}_AAL_Reg{i:02d}.nii") for i in idx]

# Harvard-Oxford maxprob labels are the XML index + 1.
HO_LAB = {"insula": 2, "cing_ant": 29, "central_operc": 42, "precentral": 7,
          "postcentral": 17, "sma": 26}

ROIS = [
    # set, roi, label, hypotheses, kind, spec
    ("caren", "AInsL", "Left anterior insula (salience-network insula)", "H2;H3;H4",
     "caren", caren_regions("RSN01", [33])),
    ("caren", "AInsR", "Right anterior insula (salience-network insula)", "H2;H3;H4",
     "caren", caren_regions("RSN01", [34])),
    ("caren", "dACC", "Dorsal ACC / anterior MCC (salience-network cingulate)", "H2;H3;H4;H6-secondary",
     "caren", caren_regions("RSN01", [37, 38, 155, 156])),
    ("caren", "PIns", "Posterior insula, bilateral (somatomotor-network insula)", "H1",
     "caren", caren_regions("RSN03", [33, 34])),
    ("caren", "COper", "Central operculum, bilateral (somatomotor-network Rolandic operculum)", "H1",
     "caren", caren_regions("RSN03", [13, 14])),
    ("caren", "SM", "Sensorimotor cortex, bilateral (somatomotor-network pre/postcentral)", "H1",
     "caren", caren_regions("RSN03", [1, 2, 61, 62])),
    ("caren", "SMA", "Supplementary motor area, bilateral (somatomotor-network SMA)", "descriptive",
     "caren", caren_regions("RSN03", [15, 16])),
    ("ieat", "ACCsphere", "IEAT ACC sphere, 10 mm at MNI (-4, +12, +36)", "H6",
     "sphere", ((-4.0, 12.0, 36.0), 10.0)),
    # Harvard-Oxford atlas-sensitivity set (reported, never used for inference)
    ("ho", "AInsL", "HO insula, x < 0, y > 0", "atlas-sensitivity",
     "ho", ("insula", "L", "ant")),
    ("ho", "AInsR", "HO insula, x > 0, y > 0", "atlas-sensitivity",
     "ho", ("insula", "R", "ant")),
    ("ho", "dACC", "HO cingulate gyrus, anterior division", "atlas-sensitivity",
     "ho", ("cing_ant", None, None)),
    ("ho", "PIns", "HO insula, y <= 0, bilateral", "atlas-sensitivity",
     "ho", ("insula", None, "post")),
    ("ho", "COper", "HO central opercular cortex, bilateral", "atlas-sensitivity",
     "ho", ("central_operc", None, None)),
    ("ho", "SM", "HO pre- and postcentral gyri, bilateral", "atlas-sensitivity",
     "ho", (("precentral", "postcentral"), None, None)),
]

NETWORKS = [("net", n, f"CAREN {n} network", "H7 (proposed)", "caren",
             [os.path.join(CAREN, f)]) for n, f in
            [("SAL", "CAREN1_SAL.nii"), ("ECN", "CAREN2_ECN.nii"), ("SMN", "CAREN3_SMN.nii"),
             ("DMN", "CAREN4_DMN.nii"), ("VIS", "CAREN5_VIS.nii")]]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_list():
    rows = []
    for fn in ("firstlevel_manifest.csv", "firstlevel_manifest_candidates.csv"):
        with open(os.path.join(REPO, "FirstLevel", fn)) as f:
            rows += list(csv.DictReader(f))
    keep = [r for r in rows if r["motion_gate_fail"] == "FALSE" and r["subject"] != "sub-056"]
    subs = {r["subject"] for r in keep}
    if len(keep) != 90 or len(subs) != 47:
        sys.exit(f"run list is {len(keep)} runs / {len(subs)} participants, expected 90 / 47")
    return keep


def grid_mm(shape, affine):
    ijk = np.indices(shape).reshape(3, -1)
    return (affine @ np.vstack([ijk, np.ones(ijk.shape[1])]))[:3]


def resample_nn(path, shape, affine):
    """Sample a source image at the target voxel centres (sform, nearest neighbour)."""
    img = nib.load(path)
    src = np.asanyarray(img.dataobj)
    if src.ndim == 4:
        src = src[..., 0]
    vox = np.linalg.inv(img.affine) @ np.vstack([grid_mm(shape, affine), np.ones(np.prod(shape))])
    idx = np.rint(vox[:3]).astype(int)
    ok = np.all((idx >= 0) & (idx < np.array(src.shape)[:, None]), axis=0)
    out = np.zeros(np.prod(shape), dtype=src.dtype)
    out[ok] = src[idx[0, ok], idx[1, ok], idx[2, ok]]
    return out.reshape(shape)


def build(kind, spec, shape, affine, mm):
    if kind == "caren":
        m = np.zeros(shape, bool)
        for p in spec:
            m |= resample_nn(p, shape, affine) > 0
        return m
    if kind == "sphere":
        (cx, cy, cz), r = spec
        d2 = (mm[0] - cx) ** 2 + (mm[1] - cy) ** 2 + (mm[2] - cz) ** 2
        return (d2 <= r * r).reshape(shape)
    if kind == "ho":
        labs, hemi, part = spec
        labs = labs if isinstance(labs, tuple) else (labs,)
        ho = resample_nn(HO, shape, affine)
        m = np.isin(ho, [HO_LAB[l] for l in labs])
        x, y = mm[0].reshape(shape), mm[1].reshape(shape)
        if hemi == "L":
            m &= x < 0
        elif hemi == "R":
            m &= x > 0
        if part == "ant":
            m &= y > 0
        elif part == "post":
            m &= y <= 0
        return m
    raise ValueError(kind)


def write_mask(data, ref, path):
    hdr = ref.header.copy()
    hdr.set_data_dtype(np.uint8)
    img = nib.Nifti1Image(data.astype(np.uint8), ref.affine, hdr)
    img.set_sform(ref.affine, code=int(ref.header["sform_code"]))
    img.set_qform(ref.affine, code=int(ref.header["qform_code"]))
    nib.save(img, path)


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = run_list()
    paths = [os.path.join(FMRIPREP, r["subject"], "func",
                          f"{r['subject']}_{r['bids_entity']}_{SPACE}_desc-brain_mask.nii.gz")
             for r in runs]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit("missing brain masks:\n" + "\n".join(missing))
    ref = nib.load(paths[0])
    shape, affine = ref.shape[:3], ref.affine
    group = np.ones(shape, bool)
    for p in paths:
        im = nib.load(p)
        if im.shape[:3] != shape or not np.allclose(im.affine, affine):
            sys.exit(f"grid mismatch: {p}")
        group &= np.asanyarray(im.dataobj) > 0
    write_mask(group, ref, os.path.join(OUT, "group_brain_mask.nii"))
    mm = grid_mm(shape, affine)

    sources = {HO}
    manifest, tested = [], {}
    for s, roi, label, hyp, kind, spec in ROIS + NETWORKS:
        if kind == "caren":
            sources.update(spec)
        raw = build(kind, spec, shape, affine, mm)
        m = raw & group
        name = f"{s}_{roi}"
        path = os.path.join(OUT, f"{name}.nii")
        write_mask(m, ref, path)
        tested[name] = m
        c = mm[:, m.ravel()].mean(axis=1) if m.any() else [np.nan] * 3
        data_hash = hashlib.sha256(m.astype(np.uint8).tobytes() + affine.tobytes()).hexdigest()
        manifest.append({
            "set": s, "roi": roi, "file": f"ROIs/{name}.nii", "label": label,
            "hypotheses": hyp, "n_vox_atlas": int(raw.sum()), "n_vox_tested": int(m.sum()),
            "coverage": round(m.sum() / raw.sum(), 4) if raw.sum() else "",
            "centroid_x": round(c[0], 1), "centroid_y": round(c[1], 1), "centroid_z": round(c[2], 1),
            "sha256_file": sha256_file(path), "sha256_data": data_hash})

    with open(os.path.join(OUT, "roi_manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0]))
        w.writeheader()
        w.writerows(manifest)

    names = list(tested)
    with open(os.path.join(OUT, "roi_overlap.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["roi_a", "roi_b", "shared_vox"])
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                if a.split("_")[0] == b.split("_")[0] or "ieat" in (a.split("_")[0], b.split("_")[0]):
                    w.writerow([a, b, int((tested[a] & tested[b]).sum())])

    with open(os.path.join(OUT, "roi_sources.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "sha256"])
        for p in sorted(sources):
            rel = os.path.relpath(p, REPO) if p.startswith(REPO) else p.replace(FSLDIR, "$FSLDIR")
            w.writerow([rel, sha256_file(p)])

    print(f"group mask: {int(group.sum())} voxels from {len(paths)} runs")
    for r in manifest:
        print(f"{r['set']:6s} {r['roi']:10s} atlas {r['n_vox_atlas']:6d} tested {r['n_vox_tested']:6d} "
              f"cov {r['coverage']}  ({r['centroid_x']}, {r['centroid_y']}, {r['centroid_z']})")


if __name__ == "__main__":
    main()
