#!/usr/bin/env python
"""
intero2023_explore_rating.py -- EXPLORATORY E8 (deviations.md, "After the update", 2026-11-26):
does the post-epoch rating / report period carry the decision signal?

  python Scripts/intero2023_explore_rating.py qc               # design only, no BOLD
  python Scripts/intero2023_explore_rating.py subject sub-057  # joint two-run first level (+ per-run primary)
  python Scripts/intero2023_explore_rating.py group            # logged group tests

The rating window. Each paced and each free epoch is followed, with no gap, by one fixed-length
rating window (v2: 12.24 s, the same whatever the answer; v1: about 16.3 s, a few 24.3 s). In it
the task runs speedrating() (Faster / Slower / No Change), then confrating() (confidence 1-7),
then arousalrating() (1-7); v1 adds valence. All responses are made on a button box under the
right hand. In v2 the sub-screens are not timestamped, so the window is modelled whole, as in
Model A (R/intero_events.R).

Models (per run, then joined as for H3/E3/E7 by build_joint_design: shared task columns,
run-specific nuisance and intercept; PSC per run, AR(1), 6 mm on the group mask, float32).
Model A's epoch regressors (paced_aware, paced_unaware, paced_wrongdir, paced_nochange, free),
subrating and every nuisance column are unchanged ('spm + derivative'). Only the rating
regressor is replaced:

  E8 report model  rating after a CHANGE epoch (paced_aware / _unaware / _wrongdir), split by
                   the report:  rating_chg_faster, rating_chg_slower, rating_chg_nochange
                   (= the misses); rating after a No Change epoch: rating_ncepoch (any report,
                   no interest); rating after free: rating_free (no interest).
  E8 conf model    (separate, as Model B is separate from Model A) rating after every paced epoch
                   as one regressor rating_paced + a confidence modulator rating_paced_confpm
                   (confidence of the preceding epoch, centred within participant across runs);
                   rating_free as above. No report split.

Contrasts:
  PRIMARY    E8chgNC = reported change - reported No Change, change epochs only
                     = (n_F * rating_chg_faster + n_S * rating_chg_slower) / (n_F + n_S)
                       - rating_chg_nochange          (n_F, n_S pooled over the participant's runs)
  secondary  E8FS    = rating_chg_faster - rating_chg_slower   (motor / key-mapping diagnostic)
  secondary  E8conf  = rating_paced_confpm                      (graded confidence, report-blind)
Inclusion (fixed before fitting): a participant enters a contrast if each of its cells has >= 2
epochs pooled across runs (E8chgNC: hit + wrong-direction >= 2 and miss >= 2; E8FS: reported
Faster >= 2 and reported Slower >= 2 on change epochs; E8conf: >= 2 paced rating windows and
confidence not constant across them) and the contrast is estimable in the joint design.

Group: per H3 mask, sign-flip permutation voxel max-t FWE, two-sided, 10,000 permutations,
seed 20261126, Bonferroni x3; mask-mean dz and JZS BF10 (r = 0.707; 0.5, 1.0 robustness).
Whole brain, two-sided cluster FWE (p < .001 per tail): E8chgNC, and E8FS (diagnostic).
Robustness (descriptive): run-1 vs run-2 sign agreement of the E8chgNC mask means, in two-run
participants whose every run has >= 1 epoch in each primary cell (per-run fits).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
import intero2023_glm_common as g
from intero2023_explore_fir import vif, rank_def
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
TAG = "E8"
CHANGE = q.CHANGE_ONLY
MIN_CELL = 2
CELLS = ("rating_chg_faster", "rating_chg_slower", "rating_chg_nochange")


def rating_events(sub, ent, n_vol):
    """Model A events (clipped as in read_run) with every rating window labelled by the epoch it
    follows (0 s gap) and by the report given in it."""
    ev, cf, _ = g.read_run_physio(sub, ent, n_vol, [])
    ev = ev.sort_values("onset").reset_index(drop=True)
    lab, conf = [], []
    for i, e in ev.iterrows():
        if e.condition != "rating":
            lab.append(e.condition); conf.append(np.nan); continue
        p = ev.loc[i - 1]
        assert abs(e.onset - (p.onset + p.duration)) < 0.01, (sub, ent, e.onset)
        assert p.condition in q.PACED_MOD + ["free"], (sub, ent, p.condition)
        assert e.response == p.response, (sub, ent, e.onset)
        if p.condition == "free":
            lab.append("rating_free")
        elif p.condition == "paced_nochange":
            lab.append("rating_ncepoch")
        else:
            assert e.response in ("Faster", "Slower", "No Change"), (sub, ent, e.response)
            lab.append({"Faster": "rating_chg_faster", "Slower": "rating_chg_slower",
                        "No Change": "rating_chg_nochange"}[e.response])
        conf.append(float(p.confidence) if p.condition != "free" else np.nan)
    ev["label"] = lab; ev["prev_conf"] = conf
    prev = ev.condition.shift(1)
    assert (prev[ev.label == "rating_chg_nochange"] == "paced_unaware").all()          # misses
    assert prev[ev.label.isin(["rating_chg_faster", "rating_chg_slower"])].isin(["paced_aware", "paced_wrongdir"]).all()
    return ev, cf


def report_design(ev, cf, n_vol):
    e = ev[["onset", "duration", "label"]].rename(columns={"label": "trial_type"}).assign(modulation=1.0)
    return q.design(e, cf, n_vol)


def conf_design(ev, cf, n_vol, conf_mean):
    lab = ev.label.where(~ev.label.isin(CELLS + ("rating_ncepoch",)), "rating_paced")
    e = ev[["onset", "duration"]].assign(trial_type=lab, modulation=1.0)
    pm = ev[lab == "rating_paced"]
    pm = pm[["onset", "duration"]].assign(trial_type="rating_paced_confpm",
                                          modulation=pm.prev_conf.values - conf_mean)
    if len(pm) and np.abs(pm.modulation).max() > 0:
        e = pd.concat([e, pm], ignore_index=True)
    return q.design(e, cf, n_vol)


def load_subject(sub, runs):
    evs = []
    for r in runs.itertuples():
        ev, cf = rating_events(sub, r.bids_entity, int(r.n_vol)); evs.append((r.bids_entity, int(r.n_vol), ev, cf))
    allev = pd.concat([e for _, _, e, _ in evs])
    pc = allev[allev.label.isin(CELLS + ("rating_ncepoch",))]
    counts = {c: int((allev.label == c).sum()) for c in CELLS + ("rating_ncepoch", "rating_free")}
    conf = pc.prev_conf.dropna()
    return evs, counts, float(conf.mean()) if len(conf) else np.nan, conf


def primary_spec(counts):
    nF, nS, nN = counts["rating_chg_faster"], counts["rating_chg_slower"], counts["rating_chg_nochange"]
    if nF + nS == 0 or nN == 0:
        return None
    s = {"rating_chg_nochange": -1.0}
    if nF: s["rating_chg_faster"] = nF / (nF + nS)
    if nS: s["rating_chg_slower"] = nS / (nF + nS)
    return s


def specs_and_inclusion(counts, conf):
    """The logged contrasts and the logged >= 2-per-cell inclusion rule."""
    out = {}
    if counts["rating_chg_faster"] + counts["rating_chg_slower"] >= MIN_CELL and counts["rating_chg_nochange"] >= MIN_CELL:
        out["E8chgNC"] = primary_spec(counts)
    if counts["rating_chg_faster"] >= MIN_CELL and counts["rating_chg_slower"] >= MIN_CELL:
        out["E8FS"] = {"rating_chg_faster": 1.0, "rating_chg_slower": -1.0}
    if len(conf) >= MIN_CELL and conf.std() > 0:
        out["E8conf"] = {"rating_paced_confpm": 1.0}
    return out


def designs(sub, runs):
    evs, counts, cmean, conf = load_subject(sub, runs)
    R = [report_design(ev, cf, n) for _, n, ev, cf in evs]
    C = [conf_design(ev, cf, n, cmean) for _, n, ev, cf in evs]
    return evs, counts, conf, R, C, hj.build_joint_design(R), hj.build_joint_design(C)


def censored_fraction(dm, name):
    cens = [c for c in dm.columns if "spike" in c or "non_steady_state" in c]
    if name not in dm.columns:
        return np.nan
    x = np.abs(dm[name].values)
    hit = (dm[cens].values != 0).any(1) if cens else np.zeros(len(dm), bool)
    return float(x[hit].sum() / x.sum()) if x.sum() > 0 else np.nan


def corr_of(dm, a, b):
    """Correlation of two contrast estimates under the design: c1'Vc2 / sqrt(c1'Vc1 c2'Vc2)."""
    X = dm.values.astype(float); sd = X.std(0); sd[sd < 1e-12] = 1.0; Xs = X / sd
    V = np.linalg.pinv(Xs.T @ Xs, rcond=1e-12)
    ca = np.array([a.get(k, 0.0) for k in dm.columns]) / sd
    cb = np.array([b.get(k, 0.0) for k in dm.columns]) / sd
    return float(ca @ V @ cb / np.sqrt((ca @ V @ ca) * (cb @ V @ cb)))


def qc():
    os.chdir(g.REPO)
    man = g.run_list(); rows, jrows = [], []
    names = list(CELLS) + ["rating_ncepoch", "rating_free", "paced_aware", "paced_unaware", "paced_wrongdir",
                           "paced_nochange", "free"]
    for sub, runs in man.groupby("subject"):
        evs, counts, conf, R, C, JR, JC = designs(sub, runs)
        for (ent, n, ev, cf), dr, dc in zip(evs, R, C):
            dA = q.model_a(ev[ev.condition.isin(q.CONDITIONS)], cf, n)
            v = vif(dr, names); vc = vif(dc, ["rating_paced", "rating_paced_confpm"]); va = vif(dA, ["rating"])
            rec = dict(subject=sub, entity=ent, n_reg=dr.shape[1], rank_deficient_report=int(rank_def(dr)),
                       rank_deficient_conf=int(rank_def(dc)),
                       **{f"n_{c}": int((ev.label == c).sum()) for c in CELLS + ("rating_ncepoch", "rating_free")},
                       vif_modelA_rating=va.get("rating"),
                       **{f"vif_{c}": v.get(c) for c in names},
                       vif_rating_paced=vc.get("rating_paced"), vif_confpm=vc.get("rating_paced_confpm"),
                       **{f"censfrac_{c}": censored_fraction(dr, c) for c in CELLS})
            s = primary_spec({c: rec[f"n_{c}"] for c in CELLS})
            rec["var_E8chgNC_run"] = hj.cvar(dr, s) if s else np.nan
            rows.append(rec)
        sp = specs_and_inclusion(counts, conf)
        prim = primary_spec(counts)
        epoch = None
        nA, nW, nU = [int(sum((ev.condition == c).sum() for _, _, ev, _ in evs)) for c in CHANGE]
        if nU and nA + nW:
            epoch = {"paced_unaware": -1.0}
            if nA: epoch["paced_aware"] = nA / (nA + nW)
            if nW: epoch["paced_wrongdir"] = nW / (nA + nW)
        v = vif(JR, names); vc = vif(JC, ["rating_paced", "rating_paced_confpm"])
        jr = dict(subject=sub, n_runs=len(evs), n_reg_report=JR.shape[1], n_reg_conf=JC.shape[1],
                  rank_deficient_report=int(rank_def(JR)), rank_deficient_conf=int(rank_def(JC)),
                  **{f"n_{c}": counts[c] for c in counts},
                  n_paced_conf=len(conf), conf_sd=float(conf.std()) if len(conf) > 1 else np.nan,
                  **{f"vif_{c}": v.get(c) for c in names[:5]},
                  vif_rating_paced=vc.get("rating_paced"), vif_confpm=vc.get("rating_paced_confpm"),
                  max_censfrac_cells=np.nanmax([censored_fraction(JR, c) for c in CELLS]) if any(c in JR for c in CELLS) else np.nan,
                  var_E8chgNC=hj.cvar(JR, prim) if prim else np.nan,
                  var_E8FS=hj.cvar(JR, {"rating_chg_faster": 1.0, "rating_chg_slower": -1.0}),
                  var_E8conf=hj.cvar(JC, {"rating_paced_confpm": 1.0}),
                  var_epoch_chgNC=hj.cvar(JR, epoch) if epoch else np.nan,
                  corr_E8chgNC_epoch=corr_of(JR, prim, epoch) if prim and epoch else np.nan,
                  **{f"in_{k}": k in sp and np.isfinite(hj.cvar(JR if k != "E8conf" else JC, sp[k])) for k in
                     ("E8chgNC", "E8FS", "E8conf")})
        # rejected alternative: all paced epochs, split by report only (any change report vs any No Change)
        nch = counts["rating_chg_faster"] + counts["rating_chg_slower"]
        nfa = int(sum(((ev.label == "rating_ncepoch") & (ev.response != "No Change")).sum() for _, _, ev, _ in evs))
        ncr = counts["rating_ncepoch"] - nfa
        jr.update(n_false_alarm=nfa, n_correct_rejection=ncr,
                  alt_allpaced_eligible=(nch + nfa >= MIN_CELL) and (counts["rating_chg_nochange"] + ncr >= MIN_CELL))
        jrows.append(jr)
    d = pd.DataFrame(rows); j = pd.DataFrame(jrows)
    d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_runs.csv"), index=False)
    j.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_joint.csv"), index=False)
    print(f"=== per-run report designs: {len(d)} runs, {d.subject.nunique()} participants ===")
    print(f"rank-deficient: report {d.rank_deficient_report.sum()}, conf {d.rank_deficient_conf.sum()}")
    print("cells per run (median [max]):", {c: f"{d[f'n_{c}'].median():.0f} [{d[f'n_{c}'].max()}]"
                                            for c in CELLS + ("rating_ncepoch", "rating_free")})
    print("Model A rating VIF median %.2f [%.2f, %.2f]" % (d.vif_modelA_rating.median(), d.vif_modelA_rating.min(),
                                                          d.vif_modelA_rating.max()))
    for c in names:
        x = d[f"vif_{c}"].dropna()
        print(f"  VIF per run {c:22s} median {x.median():6.2f}  95th {x.quantile(.95):6.2f}  max {x.max():6.2f}  (runs {len(x)})")
    x = d.vif_confpm.dropna(); print(f"  VIF per run confpm median {x.median():.2f} max {x.max():.2f}")
    print("max censored fraction of a cell regressor (per run):", {c: round(d[f'censfrac_{c}'].max(), 3) for c in CELLS})
    print(f"\n=== joint designs: {len(j)} participants ===")
    print(f"rank-deficient: report {j.rank_deficient_report.sum()}, conf {j.rank_deficient_conf.sum()}")
    for c in list(counts):
        print(f"  pooled {c:22s} median {j[f'n_{c}'].median():.0f}  [{j[f'n_{c}'].min()}, {j[f'n_{c}'].max()}]  "
              f"participants with >= 2: {(j[f'n_{c}'] >= 2).sum()}")
    print(f"  false alarms {j.n_false_alarm.sum()} / correct rejections {j.n_correct_rejection.sum()} on No Change epochs")
    for c in names[:5] + ["rating_paced", "confpm"]:
        x = j[f"vif_{c}"].dropna()
        print(f"  VIF joint {c:22s} median {x.median():6.2f}  max {x.max():6.2f}")
    for k in ("E8chgNC", "E8FS", "E8conf"):
        inn = j[j[f"in_{k}"]]
        print(f"{k}: N = {len(inn)}; variance multiplier median {inn[f'var_{k}'].median():.3f} "
              f"[{inn[f'var_{k}'].min():.3f}, {inn[f'var_{k}'].max():.3f}]; out: {sorted(j.subject[~j[f'in_{k}']])}")
    inn = j[j.in_E8chgNC]
    print(f"epoch-level (aware+wrongdir) - unaware in the same model: variance median {inn.var_epoch_chgNC.median():.3f}; "
          f"corr(E8chgNC, epoch contrast) median {inn.corr_E8chgNC_epoch.median():.2f} "
          f"[{inn.corr_E8chgNC_epoch.min():.2f}, {inn.corr_E8chgNC_epoch.max():.2f}]")
    print(f"rejected alternative (all paced epochs, report split only): eligible N = {j.alt_allpaced_eligible.sum()}")
    tw = j[(j.n_runs == 2) & j.in_E8chgNC].subject
    ok = [s for s in tw if d[(d.subject == s)].var_E8chgNC_run.notna().all()]
    print(f"robustness: two-run participants in E8chgNC {len(tw)}; with the primary estimable in both runs: {len(ok)}")


def subject(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    man = g.run_list(); runs = man[man.subject == sub]; assert len(runs), sub
    od = os.path.join(g.OUT, "firstlevel", sub); os.makedirs(od, exist_ok=True)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    evs, counts, conf, R, C, JR, JC = designs(sub, runs)
    sp = specs_and_inclusion(counts, conf)
    sp = {k: s for k, s in sp.items() if np.isfinite(hj.cvar(JR if k != "E8conf" else JC, s))}
    log = dict(subject=sub, counts=counts, n_conf=len(conf), specs=sp, n_reg_report=JR.shape[1], n_reg_conf=JC.shape[1])
    if not sp:
        log["contrasts"] = []
        json.dump(log, open(os.path.join(od, f"{sub}_e8.json"), "w")); print(sub, "nothing estimable", log); return
    Ys = []
    for r in runs.itertuples():
        Y = masker.transform(g.bold_path(sub, r.bids_entity)); assert Y.shape[0] == int(r.n_vol)
        Ys.append(mean_scaling(Y, axis=0)[0])
    def save(arr, name):
        if "variance" in name and not (arr > 0).all():
            raise ValueError(f"{sub} {name}: {(arr <= 0).sum()} non-positive variances")
        g.save_float32(masker.inverse_transform(arr), os.path.join(od, name))
    done, dofs = [], {}
    rep = {k: s for k, s in sp.items() if k != "E8conf"}
    if rep:
        out, dofs["report"] = fit(np.vstack(Ys), JR, rep)
        for k, (e, v) in out.items():
            save(e, f"{sub}_var-base_con-{k}joint_effect.nii.gz"); save(v, f"{sub}_var-base_con-{k}joint_variance.nii.gz"); done.append(k)
    if "E8conf" in sp:
        out, dofs["conf"] = fit(np.vstack(Ys), JC, {"E8conf": sp["E8conf"]})
        e, v = out["E8conf"]
        save(e, f"{sub}_var-base_con-E8confjoint_effect.nii.gz"); save(v, f"{sub}_var-base_con-E8confjoint_variance.nii.gz"); done.append("E8conf")
    # robustness: the primary per run, two-run participants, primary cells present in each run
    perrun = []
    if "E8chgNC" in sp and len(runs) == 2:
        specs_r = []
        for (ent, n, ev, cf), dr in zip(evs, R):
            s = primary_spec({c: int((ev.label == c).sum()) for c in CELLS})
            specs_r.append(s if s and np.isfinite(hj.cvar(dr, s)) else None)
        if all(specs_r):
            for k, ((ent, n, ev, cf), dr, s, Y) in enumerate(zip(evs, R, specs_r, Ys)):
                out, _ = fit(Y, dr, {"E8chgNC": s})
                e, v = out["E8chgNC"]
                save(e, f"{sub}_var-base_con-E8chgNCrun{k + 1}_effect.nii.gz"); perrun.append(ent)
    log.update(dof=dofs, contrasts=done, perrun_E8chgNC=perrun)
    json.dump(log, open(os.path.join(od, f"{sub}_e8.json"), "w"))
    print(sub, "done", done, "per-run", perrun, "dof", dofs)


def group():
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.second_level import non_parametric_inference
    import intero2023_group as G
    from intero2023_explore_awareness import jzs_bf10
    man = g.run_list(); allsubs = sorted(man.subject.unique()); rows = []
    for name, role, wb in (("E8chgNC", "primary", True), ("E8FS", "secondary (motor / key-mapping diagnostic)", True),
                           ("E8conf", "secondary", False)):
        subs = [s for s in allsubs if os.path.exists(G.fl(s, "base", f"{name}joint", "effect"))]
        files = [G.fl(s, "base", f"{name}joint", "effect") for s in subs]
        imgs = [nib.load(f) for f in files]
        dm = pd.DataFrame({"intercept": np.ones(len(files))})
        json.dump(dict(N=len(subs), subjects=subs), open(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_subjects.json"), "w"))
        for roi in G.SETS["H3"]:
            rmask = nib.load(G.ROI(roi)).get_fdata() > 0
            lp = non_parametric_inference(imgs, design_matrix=dm, second_level_contrast="intercept", mask=G.ROI(roi),
                                          n_perm=G.N_PERM, two_sided_test=True, random_state=G.SEED, n_jobs=G.N_JOBS)
            tmap = G.ols_t(imgs, dm, "intercept", G.ROI(roi))
            lpv = lp.get_fdata()[rmask]; tv = tmap.get_fdata()[rmask]
            k = int(np.nanargmax(np.abs(tv))); ijk = np.argwhere(rmask)[k]
            pmin = float(10 ** -np.nanmax(lpv))
            mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
            y = np.array([mk.transform(f).mean() for f in files])
            row = dict(analysis=f"EXPLORATORY {TAG}", contrast=name, role=role, roi=roi, N=len(files),
                       peak_abs_t_signed=round(float(tv[k]), 3),
                       peak_mni=[round(float(c), 1) for c in nib.affines.apply_affine(tmap.affine, ijk)],
                       p_fwe=pmin, p_fwe_bonf3=min(1.0, 3 * pmin), n_vox_p_lt_05=int((lpv > -np.log10(0.05)).sum()),
                       mask_mean=float(y.mean()), mask_mean_dz=float(y.mean() / y.std(ddof=1)),
                       mask_mean_t=float(stats.ttest_1samp(y, 0).statistic))
            for r in (0.5, 0.707, 1.0):
                row[f"BF10_r{r}"] = jzs_bf10(row["mask_mean_t"], len(y), r)
            row["BF01_r0.707"] = 1 / row["BF10_r0.707"]
            rows.append(row); print({k_: v_ for k_, v_ in row.items()})
        if wb:
            od = os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}_{name}"); os.makedirs(od, exist_ok=True)
            res = non_parametric_inference(files, design_matrix=dm, second_level_contrast="intercept", mask=g.GROUP_MASK,
                                           n_perm=G.N_PERM, two_sided_test=True, threshold=0.002,
                                           random_state=G.SEED, n_jobs=G.N_JOBS)
            for kk in ("t", "logp_max_size"):
                g.save_float32(res[kk], os.path.join(od, f"wholebrain_{kk}.nii.gz"))
            cl = G.cluster_table(res["t"], res["logp_max_size"], len(files) - 1)
            if len(cl):
                cl = cl.assign(labels=[atlas_labels(m) for m in cl.peak_mni])
            cl.insert(0, "contrast", name)
            cl.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_wholebrain.csv"), index=False)
            print(name, "whole brain:", cl.head(15).to_string() if len(cl) else "no cluster")
    pd.DataFrame(rows).to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_rating_h3masks.csv"), index=False)
    robustness()


def atlas_labels(mni):
    """FSL atlasquery at the peak: Harvard-Oxford cortical / subcortical and the cerebellar atlas
    (descriptive)."""
    import intero2023_group as G
    out = []
    for atlas in ("Harvard-Oxford Cortical Structural Atlas", "Harvard-Oxford Subcortical Structural Atlas",
                  "Cerebellar Atlas in MNI152 space after normalization with FLIRT"):
        try:
            txt = G.sh(f'atlasquery -a "{atlas}" -c {mni[0]},{mni[1]},{mni[2]}').strip()
        except Exception as ex:
            txt = f"<error {ex}>"
        lab = txt.split("<br>", 1)[-1] if "<br>" in txt else txt
        if lab and "No label found" not in lab:
            out.append(lab.strip())
    return "; ".join(out)


def robustness():
    """Descriptive: run-1 vs run-2 sign agreement of the E8chgNC mask means (two-run participants)."""
    from nilearn.maskers import NiftiMasker
    import intero2023_group as G
    man = g.run_list(); subs = []
    for s in sorted(man.subject.unique()):
        f = [os.path.join(g.OUT, "firstlevel", s, f"{s}_var-base_con-E8chgNCrun{k}_effect.nii.gz") for k in (1, 2)]
        if all(os.path.exists(x) for x in f):
            subs.append((s, f))
    rows = []
    for roi in G.SETS["H3"]:
        mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
        y = np.array([[mk.transform(x).mean() for x in f] for _, f in subs])
        agree = np.sign(y[:, 0]) == np.sign(y[:, 1])
        rows.append(dict(analysis=f"EXPLORATORY {TAG} robustness (descriptive)", contrast="E8chgNC", roi=roi, N=len(subs),
                         run1_mean=float(y[:, 0].mean()), run1_dz=float(y[:, 0].mean() / y[:, 0].std(ddof=1)),
                         run2_mean=float(y[:, 1].mean()), run2_dz=float(y[:, 1].mean() / y[:, 1].std(ddof=1)),
                         sign_agreement=float(agree.mean()), binom_p_two_sided=float(stats.binomtest(int(agree.sum()), len(agree)).pvalue),
                         r_run1_run2=float(np.corrcoef(y[:, 0], y[:, 1])[0, 1])))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_robustness_runs.csv"), index=False)
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    a = sys.argv[1]
    {"qc": qc, "subject": lambda: subject(sys.argv[2]), "group": group}[a]()
