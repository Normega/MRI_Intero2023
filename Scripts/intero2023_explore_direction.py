#!/usr/bin/env python
"""
intero2023_explore_direction.py -- EXPLORATORY E10 (deviations.md, "After the update", 2026-11-27):
is E8's rating-period Slower > Faster effect about the REPORT, or about the STIMULUS / the
preceding paced epoch? Secondary: an epoch-level stimulus direction x noticed interaction.

  python Scripts/intero2023_explore_direction.py qc               # design only, no BOLD
  python Scripts/intero2023_explore_direction.py subject sub-057  # joint two-run first level
  python Scripts/intero2023_explore_direction.py group            # logged group tests

Background (E8, seen). In E8 the paced epoch regressors were split by outcome (aware / unaware /
wrong-direction / No Change), not by stimulus direction. The rating window that follows each
epoch with no gap was split by the report. Reported Faster - reported Slower (E8FS, N 32) was
Slower > Faster in both anterior insulae (p_FWE x3 .004, .017). About 80% of each report cell
shares the stimulus direction, and Slower epochs are longer (23 s vs 17.5 s) and breathed more
slowly, so the effect may be carried over from the epoch.

Models (per run, joined as for H3/E3/E7/E8 by build_joint_design: shared task columns,
run-specific nuisance and intercept; PSC per run, AR(1), 6 mm on the group mask, float32,
'spm + derivative'). Free, subrating and every nuisance column are Model A's, unchanged.

  Epoch regressors (both models): every paced CHANGE epoch is split by stimulus direction x
  outcome: paced_Fhit, paced_Fmiss, paced_Fwrong, paced_Shit, paced_Smiss, paced_Swrong
  (F/S = stimulus Faster/Slower; hit = correct direction reported, miss = No Change reported,
  wrong = opposite direction reported); paced_nochange; free. A cell with no epoch in a run is
  simply absent from that run.
  Model P (primary): the rating windows are split as in E8: after change epochs by the report,
  rating_chg_faster, rating_chg_slower, rating_chg_nochange (= the misses); rating_ncepoch;
  rating_free.
  Model W (dissociation, descriptive): the rating windows after change epochs are split by
  stimulus x report: rating_FF, rating_FS, rating_SF, rating_SS, rating_FN, rating_SN (first
  letter = stimulus, second = report); rating_ncepoch; rating_free.

Contrasts (fixed after the design-only check, before any BOLD was read):
  PRIMARY   E10rep = rating_chg_faster - rating_chg_slower in Model P (reported Faster - reported
            Slower, change epochs, with the epochs split by stimulus direction). Inclusion as E8FS:
            reported Faster >= 2 and reported Slower >= 2 on change epochs, pooled over runs.
  secondary E10int = (paced_Fhit - paced_Shit) - (paced_Fmiss - paced_Smiss), Model P: the
            stimulus-direction effect in the epoch when noticed (hit) minus when missed.
            Inclusion: >= 1 epoch in each of the four cells (E8's >= 2 rule leaves N 6).
  secondary E10dir = Slower - Faster stimulus main effect in the epoch, Model P: trial-weighted
            mean of the S cells minus that of the F cells (diagnostic: where a carry-over would
            come from). Inclusion: >= 2 change epochs of each direction.
  secondary E10wd  (descriptive only; not a voxelwise test): Model W rating cells FF, FS, SF, SS.
            (a) In each H3 mask, a pooled across-participant WLS of the participants' mask-mean
            cell betas on report (+1/2 F, -1/2 S), stimulus (+1/2 F, -1/2 S) and their product,
            participant fixed effects, weights = epochs in the cell, cluster-robust (CR1) SE by
            participant. b_report and b_stimulus with 95% CI. (b) In participants with >= 1 SF
            and >= 1 FS epoch, the mask-mean of rating_SF - rating_FS (report F - S with the
            stimulus opposite): a report effect predicts the sign of E8FS (negative), a stimulus
            carry-over predicts the opposite sign. Mean, dz and sign count only.

Group: per H3 mask (caren_AInsL, caren_AInsR, caren_dACC), sign-flip permutation voxel max-t FWE,
two-sided, 10,000 permutations, seed 20261126, Bonferroni x3; mask-mean dz and JZS BF10
(r = 0.707; 0.5, 1.0 robustness), equal-weighted OLS. Whole brain, two-sided cluster FWE
(p < .001 per tail): E10rep and E10dir. Peaks labelled with FSL atlasquery (Harvard-Oxford
cortical / subcortical, cerebellar FNIRT).
Descriptive: E8FS vs E10rep mask means in the same participants (dz of each, ratio of means, r).

Design QC (no BOLD; `qc`): all 47 joint designs full rank in P and W. E10rep N 32, variance
multiplier median 0.143 [0.087, 0.311], 1.21x E8FS's in the E8 model for the same participants
(efficiency 0.83); rating report-cell VIF median 6.8 (max 31, as E8); corr with the epoch-level
report contrast .31, with the epoch stimulus contrast .18. E10int N 22 (6 under >= 2), variance
0.404 [0.255, 0.881]. E10dir N 47, 0.066. E10wd: 79 wrong-direction epochs in 36 participants;
only 10 have >= 1 SF and >= 1 FS (4 with >= 2), so it is not tested voxelwise.
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats
import intero2023_glm_common as g
from intero2023_explore_fir import vif, rank_def
import intero2023_explore_rating as e8
q, hj = g.q, g.hj
warnings.simplefilter("ignore")
RD = os.path.join(g.REPO, "Analysis", "results")
TAG = "E10"
EPOCH = {("Faster", "paced_aware"): "paced_Fhit", ("Faster", "paced_unaware"): "paced_Fmiss",
         ("Faster", "paced_wrongdir"): "paced_Fwrong", ("Slower", "paced_aware"): "paced_Shit",
         ("Slower", "paced_unaware"): "paced_Smiss", ("Slower", "paced_wrongdir"): "paced_Swrong"}
ECELLS = tuple(EPOCH.values())
F_CELLS, S_CELLS = ECELLS[:3], ECELLS[3:]
P_RATE = ("rating_chg_faster", "rating_chg_slower", "rating_chg_nochange")
W_RATE = ("rating_FF", "rating_FS", "rating_SF", "rating_SS", "rating_FN", "rating_SN")
LET = {"Faster": "F", "Slower": "S", "No Change": "N"}
INT_CELLS = ("paced_Fhit", "paced_Shit", "paced_Fmiss", "paced_Smiss")
E10INT = {"paced_Fhit": 1.0, "paced_Shit": -1.0, "paced_Fmiss": -1.0, "paced_Smiss": 1.0}
E10REP = {"rating_chg_faster": 1.0, "rating_chg_slower": -1.0}
CONS = ("E10rep", "E10int", "E10dir")


def label_events(sub, ent, n_vol):
    """E8's labelled events plus the stimulus-direction x outcome epoch label and the Model W
    rating label. E8's own asserts (0 s gap, report carried over) run inside rating_events."""
    ev, cf = e8.rating_events(sub, ent, n_vol)
    ep, wl = [], []
    for i, e in ev.iterrows():
        if e.condition in q.CHANGE_ONLY:
            lab = EPOCH[(e.change_direction, e.condition)]
            want = {"paced_aware": e.change_direction, "paced_unaware": "No Change",
                    "paced_wrongdir": {"Faster": "Slower", "Slower": "Faster"}[e.change_direction]}[e.condition]
            assert e.response == want, (sub, ent, e.onset, e.condition, e.response)
            ep.append(lab); wl.append(lab)
        elif e.condition == "rating" and e.label in P_RATE:
            p = ev.loc[i - 1]
            ep.append(e.label); wl.append(f"rating_{LET[p.change_direction]}{LET[e.response]}")
        else:
            ep.append(e.label); wl.append(e.label)
    ev["labP"] = ep; ev["labW"] = wl
    return ev, cf


def _design(ev, cf, n_vol, col):
    e = ev[["onset", "duration", col]].rename(columns={col: "trial_type"}).assign(modulation=1.0)
    return q.design(e, cf, n_vol)


def load(sub, runs):
    evs = []
    for r in runs.itertuples():
        ev, cf = label_events(sub, r.bids_entity, int(r.n_vol)); evs.append((r.bids_entity, int(r.n_vol), ev, cf))
    allev = pd.concat([e for _, _, e, _ in evs])
    counts = {c: int((allev.labP == c).sum()) for c in ECELLS + P_RATE}
    counts.update({c: int((allev.labW == c).sum()) for c in W_RATE})
    P = [_design(ev, cf, n, "labP") for _, n, ev, cf in evs]
    W = [_design(ev, cf, n, "labW") for _, n, ev, cf in evs]
    return evs, counts, P, W, hj.build_joint_design(P), hj.build_joint_design(W)


def dir_spec(counts):
    nF = sum(counts[c] for c in F_CELLS); nS = sum(counts[c] for c in S_CELLS)
    s = {c: counts[c] / nS for c in S_CELLS if counts[c]}
    s.update({c: -counts[c] / nF for c in F_CELLS if counts[c]})
    return s, nF, nS


def specs_and_inclusion(counts):
    """The logged contrasts and the logged inclusion rules."""
    out = {}
    if counts["rating_chg_faster"] >= 2 and counts["rating_chg_slower"] >= 2:
        out["E10rep"] = dict(E10REP)
    if all(counts[c] >= 1 for c in INT_CELLS):
        out["E10int"] = dict(E10INT)
    s, nF, nS = dir_spec(counts)
    if nF >= 2 and nS >= 2:
        out["E10dir"] = s
    return out


def qc():
    os.chdir(g.REPO)
    man = g.run_list(); rows, jrows = [], []
    namesP = list(ECELLS) + ["paced_nochange", "free"] + list(P_RATE) + ["rating_ncepoch", "rating_free"]
    namesW = list(W_RATE)
    beh = []
    for sub, runs in man.groupby("subject"):
        evs, counts, P, W, JP, JW = load(sub, runs)
        # E8 report design (as fitted for E8FS), for comparison
        J8 = hj.build_joint_design([e8.report_design(*e8.rating_events(sub, r.bids_entity, int(r.n_vol)), int(r.n_vol))
                                    for r in runs.itertuples()])
        for (ent, n, ev, cf), dp, dw in zip(evs, P, W):
            vp = vif(dp, namesP); vw = vif(dw, namesW)
            rows.append(dict(subject=sub, entity=ent, n_reg_P=dp.shape[1], n_reg_W=dw.shape[1],
                             rank_deficient_P=int(rank_def(dp)), rank_deficient_W=int(rank_def(dw)),
                             **{f"n_{c}": int((ev.labP == c).sum()) for c in ECELLS + P_RATE},
                             **{f"vif_{c}": vp.get(c) for c in namesP},
                             **{f"vifW_{c}": vw.get(c) for c in namesW}))
            for e in ev[ev.condition.isin(q.PACED_MOD)].itertuples():
                beh.append(dict(subject=sub, entity=ent, epoch_index=int(e.epoch_index), cell=e.labP,
                                duration=float(e.duration)))
        sp = specs_and_inclusion(counts)
        dsp, nF, nS = dir_spec(counts)
        vP = vif(JP, namesP); vW = vif(JW, namesW); v8 = vif(J8, ["rating_chg_faster", "rating_chg_slower"])
        # epoch-level report F - S (hit + wrong cells, as the rating report split) in Model P
        ep_rep = {}
        nrF = counts["paced_Fhit"] + counts["paced_Swrong"]; nrS = counts["paced_Shit"] + counts["paced_Fwrong"]
        if nrF and nrS:
            for c in ("paced_Fhit", "paced_Swrong"):
                if counts[c]: ep_rep[c] = counts[c] / nrF
            for c in ("paced_Shit", "paced_Fwrong"):
                if counts[c]: ep_rep[c] = -counts[c] / nrS
        neg_dir = {k: -v for k, v in dsp.items()}          # Faster - Slower stimulus, epoch
        ok_rep = "E10rep" in sp
        wd = {"rating_SF": 1.0, "rating_FS": -1.0}
        jr = dict(subject=sub, n_runs=len(evs), n_reg_P=JP.shape[1], n_reg_W=JW.shape[1],
                  rank_deficient_P=int(rank_def(JP)), rank_deficient_W=int(rank_def(JW)),
                  **{f"n_{c}": counts[c] for c in counts},
                  **{f"vif_{c}": vP.get(c) for c in namesP},
                  **{f"vifW_{c}": vW.get(c) for c in namesW},
                  vif_E8_rating_chg_faster=v8.get("rating_chg_faster"), vif_E8_rating_chg_slower=v8.get("rating_chg_slower"),
                  var_E10rep=hj.cvar(JP, E10REP), var_E8FS_E8model=hj.cvar(J8, E10REP),
                  var_E10int=hj.cvar(JP, E10INT) if "E10int" in sp else np.nan,
                  var_E10dir=hj.cvar(JP, dsp) if "E10dir" in sp else np.nan,
                  var_E10wd_SFminusFS=hj.cvar(JW, wd),
                  corr_E10rep_epochreport=e8.corr_of(JP, E10REP, ep_rep) if ok_rep and ep_rep else np.nan,
                  corr_E10rep_epochstimFS=e8.corr_of(JP, E10REP, neg_dir) if ok_rep and nF and nS else np.nan,
                  **{f"in_{k}": k in sp and np.isfinite(hj.cvar(JP, sp[k])) for k in CONS},
                  in_E10wd_pooled=any(counts[c] for c in ("rating_FS", "rating_SF")),
                  in_E10wd_pair1=counts["rating_SF"] >= 1 and counts["rating_FS"] >= 1,
                  in_E10wd_pair2=counts["rating_SF"] >= 2 and counts["rating_FS"] >= 2,
                  in_E10int_rule2=all(counts[c] >= 2 for c in INT_CELLS))
        jrows.append(jr)
    d = pd.DataFrame(rows); j = pd.DataFrame(jrows); b = pd.DataFrame(beh)
    d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_runs.csv"), index=False)
    j.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_design_qc_joint.csv"), index=False)
    # behaviour by cell (no BOLD): epoch duration and the measured breathing-period ratio
    bt = pd.read_csv("Data/breathing_trials.csv", na_values=["NA", "n/a"])
    bt = bt[bt.trial_type == "paced_breathing"][["subject", "bids_entity", "epoch_index", "actual_period_ratio", "sdt_class"]]
    b = b.merge(bt.rename(columns={"bids_entity": "entity"}), on=["subject", "entity", "epoch_index"], how="left")
    assert b.sdt_class.notna().all()
    cb = b.groupby("cell").agg(n=("duration", "size"), duration_mean=("duration", "mean"),
                               period_ratio_mean=("actual_period_ratio", "mean"),
                               period_ratio_n=("actual_period_ratio", "count")).reset_index()
    cb.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_cell_behaviour.csv"), index=False)
    print(cb.round(3).to_string(index=False))
    m = cb.set_index("cell")
    for k in ("duration_mean", "period_ratio_mean"):
        dd = (m.at["paced_Fhit", k] - m.at["paced_Shit", k]) - (m.at["paced_Fmiss", k] - m.at["paced_Smiss", k])
        print(f"double difference (Fhit - Shit) - (Fmiss - Smiss) of {k}: {dd:+.3f}")
    print(f"\n=== per-run designs: {len(d)} runs; rank-deficient P {d.rank_deficient_P.sum()}, W {d.rank_deficient_W.sum()}")
    for c in namesP:
        x = d[f"vif_{c}"].dropna()
        print(f"  VIF per run P {c:20s} median {x.median():6.2f}  95th {x.quantile(.95):6.2f}  max {x.max():7.2f}  (runs {len(x)})")
    print(f"=== joint designs: {len(j)} participants; rank-deficient P {j.rank_deficient_P.sum()}, W {j.rank_deficient_W.sum()}")
    for c in list(ECELLS) + list(P_RATE) + list(W_RATE):
        print(f"  pooled {c:20s} total {j[f'n_{c}'].sum():4d}  median {j[f'n_{c}'].median():.0f} "
              f"[{j[f'n_{c}'].min()}, {j[f'n_{c}'].max()}]  >=1: {(j[f'n_{c}'] >= 1).sum():2d}  >=2: {(j[f'n_{c}'] >= 2).sum():2d}")
    for c in namesP:
        x = j[f"vif_{c}"].dropna(); print(f"  VIF joint P {c:20s} median {x.median():6.2f}  max {x.max():7.2f}")
    for c in namesW:
        x = j[f"vifW_{c}"].dropna(); print(f"  VIF joint W {c:20s} median {x.median():6.2f}  max {x.max():7.2f}  (n {len(x)})")
    for c in ("rating_chg_faster", "rating_chg_slower"):
        x = j[f"vif_E8_{c}"].dropna(); print(f"  VIF joint E8 model {c:20s} median {x.median():6.2f}  max {x.max():7.2f}")
    for k in CONS:
        inn = j[j[f"in_{k}"]]
        print(f"{k}: N = {len(inn)}; variance multiplier median {inn[f'var_{k}'].median():.3f} "
              f"[{inn[f'var_{k}'].min():.3f}, {inn[f'var_{k}'].max():.3f}]; out: {sorted(j.subject[~j[f'in_{k}']])}")
    inn = j[j.in_E10rep]
    r = inn.var_E10rep / inn.var_E8FS_E8model
    print(f"E10rep vs E8FS (same participants, E8 model): E8 variance median {inn.var_E8FS_E8model.median():.3f}; "
          f"ratio E10/E8 median {r.median():.2f} [{r.min():.2f}, {r.max():.2f}] -> efficiency {1 / r.median():.2f}")
    print(f"corr(E10rep, epoch-level report F-S) median {inn.corr_E10rep_epochreport.median():.2f} "
          f"[{inn.corr_E10rep_epochreport.min():.2f}, {inn.corr_E10rep_epochreport.max():.2f}]")
    print(f"corr(E10rep, epoch stimulus F-S) median {inn.corr_E10rep_epochstimFS.median():.2f} "
          f"[{inn.corr_E10rep_epochstimFS.min():.2f}, {inn.corr_E10rep_epochstimFS.max():.2f}]")
    print(f"E10int under E8's >= 2 rule: N = {j.in_E10int_rule2.sum()}")
    print(f"E10wd: participants with any wrong-direction rating {j.in_E10wd_pooled.sum()}; SF>=1 & FS>=1 {j.in_E10wd_pair1.sum()}; "
          f">=2 each {j.in_E10wd_pair2.sum()}; wrong-direction epochs total {j.n_rating_FS.sum() + j.n_rating_SF.sum()}")
    x = j.var_E10wd_SFminusFS[j.in_E10wd_pair1]
    print(f"  variance multiplier SF - FS (pair1) median {x.median():.3f} [{x.min():.3f}, {x.max():.3f}]")


def subject(sub):
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.first_level.first_level import mean_scaling
    from intero2023_firstlevel import fit
    man = g.run_list(); runs = man[man.subject == sub]; assert len(runs), sub
    od = os.path.join(g.OUT, "firstlevel", sub); os.makedirs(od, exist_ok=True)
    masker = NiftiMasker(mask_img=g.GROUP_MASK, smoothing_fwhm=6, standardize=False).fit()
    evs, counts, P, W, JP, JW = load(sub, runs)
    sp = {k: s for k, s in specs_and_inclusion(counts).items() if np.isfinite(hj.cvar(JP, s))}
    wc = {f"E10W{c.split('_')[1]}": {c: 1.0} for c in ("rating_FF", "rating_FS", "rating_SF", "rating_SS")
          if c in JW.columns and np.isfinite(hj.cvar(JW, {c: 1.0}))}
    log = dict(subject=sub, counts=counts, specs=sp, W_cells=sorted(wc), n_reg_P=JP.shape[1], n_reg_W=JW.shape[1])
    Ys = []
    for r in runs.itertuples():
        Y = masker.transform(g.bold_path(sub, r.bids_entity)); assert Y.shape[0] == int(r.n_vol)
        Ys.append(mean_scaling(Y, axis=0)[0])
    Y = np.vstack(Ys)
    def save(arr, name):
        if "variance" in name and not (arr > 0).all():
            raise ValueError(f"{sub} {name}: {(arr <= 0).sum()} non-positive variances")
        g.save_float32(masker.inverse_transform(arr), os.path.join(od, name))
    done, dofs = [], {}
    for tag, J, specs in (("P", JP, sp), ("W", JW, wc)):
        if not specs:
            continue
        out, dofs[tag] = fit(Y, J, specs)
        for k, (e, v) in out.items():
            save(e, f"{sub}_var-base_con-{k}joint_effect.nii.gz"); save(v, f"{sub}_var-base_con-{k}joint_variance.nii.gz")
            done.append(k)
    log.update(dof=dofs, contrasts=done)
    json.dump(log, open(os.path.join(od, f"{sub}_e10.json"), "w"))
    print(sub, "done", done, "dof", dofs)


def group():
    from nilearn.maskers import NiftiMasker
    from nilearn.glm.second_level import non_parametric_inference
    import intero2023_group as G
    from intero2023_explore_awareness import jzs_bf10
    man = g.run_list(); allsubs = sorted(man.subject.unique()); rows = []
    for name, role, wb in (("E10rep", "primary", True),
                           ("E10int", "secondary (epoch direction x noticed)", False),
                           ("E10dir", "secondary (diagnostic: epoch stimulus Slower - Faster)", True)):
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
                cl = cl.assign(labels=[e8.atlas_labels(m) for m in cl.peak_mni])
            cl.insert(0, "contrast", name)
            cl.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_{name}_wholebrain.csv"), index=False)
            print(name, "whole brain:", cl.head(15).to_string() if len(cl) else "no cluster")
    pd.DataFrame(rows).to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_h3masks.csv"), index=False)
    compare_e8()
    dissociation()


def compare_e8():
    """Descriptive: E8FS (E8 model) vs E10rep (direction-controlled) mask means, same participants."""
    from nilearn.maskers import NiftiMasker
    import intero2023_group as G
    man = g.run_list()
    subs = [s for s in sorted(man.subject.unique()) if os.path.exists(G.fl(s, "base", "E10repjoint", "effect"))
            and os.path.exists(G.fl(s, "base", "E8FSjoint", "effect"))]
    rows = []
    for roi in G.SETS["H3"]:
        mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
        a = np.array([mk.transform(G.fl(s, "base", "E8FSjoint", "effect")).mean() for s in subs])
        b = np.array([mk.transform(G.fl(s, "base", "E10repjoint", "effect")).mean() for s in subs])
        rows.append(dict(analysis=f"EXPLORATORY {TAG} (descriptive)", roi=roi, N=len(subs),
                         E8FS_mean=float(a.mean()), E8FS_dz=float(a.mean() / a.std(ddof=1)),
                         E10rep_mean=float(b.mean()), E10rep_dz=float(b.mean() / b.std(ddof=1)),
                         ratio_of_means=float(b.mean() / a.mean()), r_across_participants=float(np.corrcoef(a, b)[0, 1])))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_vs_E8FS.csv"), index=False)
    print(d.round(3).to_string(index=False))


def dissociation():
    """E10wd, descriptive: Model W rating cells FF, FS, SF, SS, mask means per participant."""
    from nilearn.maskers import NiftiMasker
    import intero2023_group as G
    man = g.run_list(); allsubs = sorted(man.subject.unique())
    cells = ("FF", "FS", "SF", "SS"); recs = []
    for s in allsubs:
        js = os.path.join(g.OUT, "firstlevel", s, f"{s}_e10.json")
        cnt = json.load(open(js))["counts"]
        for c in cells:
            f = G.fl(s, "base", f"E10W{c}joint", "effect")
            if os.path.exists(f):
                recs.append(dict(subject=s, cell=c, n=cnt[f"rating_{c}"], file=f))
    D = pd.DataFrame(recs)
    D["report"] = np.where(D.cell.str[1] == "F", 0.5, -0.5)
    D["stimulus"] = np.where(D.cell.str[0] == "F", 0.5, -0.5)
    D["inter"] = D.report * D.stimulus * 4 / 2      # +1/2 correct, -1/2 wrong direction
    rows = []
    for roi in G.SETS["H3"]:
        mk = NiftiMasker(mask_img=G.ROI(roi)).fit()
        D[roi] = [float(mk.transform(f).mean()) for f in D.file]
        b, se, n_sub = pooled_wls(D, roi)
        Dc = D[D.groupby("subject").cell.transform("size") >= 2]
        for k in ("report", "stimulus", "inter"):
            rows.append(dict(analysis=f"EXPLORATORY {TAG} (descriptive)", part="E10wd-a pooled WLS", roi=roi,
                             term={"report": "report F - S (stimulus held)", "stimulus": "stimulus F - S (report held)",
                                   "inter": "correct - wrong direction"}[k],
                             N_participants=n_sub, N_cells=len(Dc), N_epochs=int(Dc.n.sum()),
                             N_wrongdir_epochs=int(Dc.n[Dc.cell.isin(["FS", "SF"])].sum()),
                             estimate=b[k], se_CR1=se[k], ci95_lo=b[k] - stats.t.isf(.025, n_sub - 1) * se[k],
                             ci95_hi=b[k] + stats.t.isf(.025, n_sub - 1) * se[k]))
        w = D.pivot(index="subject", columns="cell", values=roi)
        if {"SF", "FS"} <= set(w.columns):
            x = (w.SF - w.FS).dropna()
            rows.append(dict(analysis=f"EXPLORATORY {TAG} (descriptive)", part="E10wd-b SF - FS", roi=roi,
                             term="rating SF - FS (report F - S, stimulus opposite)", N_participants=len(x),
                             estimate=float(x.mean()), dz=float(x.mean() / x.std(ddof=1)) if len(x) > 1 else np.nan,
                             n_negative=int((x < 0).sum())))
    D.drop(columns="file").to_csv(os.path.join(g.OUT, "group", f"EXPLORATORY_{TAG}_wd_cells.csv"), index=False)
    out = pd.DataFrame(rows); out.to_csv(os.path.join(RD, f"EXPLORATORY_{TAG}_wd_descriptive.csv"), index=False)
    print(out.round(4).to_string(index=False))


def pooled_wls(D, y):
    """WLS of cell mask means on report, stimulus and their product with participant fixed effects,
    weights = epochs in the cell; CR1 cluster-robust SE by participant. Only participants with >= 2
    cells contribute under the fixed effects, so only they are counted (changed after the pre-fit
    commit bda9a05: the first group run also counted 1 single-cell participant, G 47 not 46)."""
    D = D[D.groupby("subject").cell.transform("size") >= 2]
    subs = sorted(D.subject.unique())
    X = np.column_stack([D[["report", "stimulus", "inter"]].values] +
                        [(D.subject == s).values.astype(float) for s in subs])
    w = D.n.values.astype(float); yv = D[y].values
    sw = np.sqrt(w); Xw = X * sw[:, None]; yw = yv * sw
    XtXi = np.linalg.pinv(Xw.T @ Xw); beta = XtXi @ Xw.T @ yw; res = yw - Xw @ beta
    meat = np.zeros((X.shape[1], X.shape[1]))
    for s in subs:
        ix = (D.subject == s).values
        u = Xw[ix].T @ res[ix]; meat += np.outer(u, u)
    G_, n, k = len(subs), len(D), np.linalg.matrix_rank(Xw)
    V = XtXi @ meat @ XtXi * (G_ / (G_ - 1)) * ((n - 1) / (n - k))
    names = ("report", "stimulus", "inter")
    return ({nm: float(beta[i]) for i, nm in enumerate(names)},
            {nm: float(np.sqrt(V[i, i])) for i, nm in enumerate(names)}, G_)


if __name__ == "__main__":
    a = sys.argv[1]
    {"qc": qc, "subject": lambda: subject(sys.argv[2]), "group": group, "wd": dissociation}[a]()
