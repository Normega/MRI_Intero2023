#!/bin/bash
# run_group_all.sh -- every registered group analysis, in order (prereg v0.8 §3, §8, §9).
# Primary analyses get the whole-brain test; sensitivity analyses are ROI tests (GROUP_WB=0).
set -u
cd "$(dirname "$0")/.."
PY=${PY:-~/miniforge3/envs/neuro/bin/python}
LOG=${INTERO_GLM_OUT:-/media/normega/FB_Cluster_Data/intero2024/derivatives/intero2023_glm}/logs/group.log
run() { echo "[$(date +%T)] $*" >> "$LOG"; env "$@" >> "$LOG" 2>&1 || echo "FAILED: $*" >> "$LOG"; }
G="$PY Scripts/intero2023_group.py"

# primary
for h in H3 H1 H2 H4 H6; do run GROUP_WB=1 $G $h base all; done
# S1 / S2 / S2b / +both (ROI tests, base and corrected in the same subset)
for h in H3 H2 H1 H4 H6; do
  run GROUP_WB=0 $G $h base S1;  run GROUP_WB=0 $G $h co2 S1
  run GROUP_WB=0 $G $h base S2;  run GROUP_WB=0 $G $h rvt S2
  run GROUP_WB=0 $G $h base S2b; run GROUP_WB=0 $G $h rvt S2b
  run GROUP_WB=0 $G $h base both; run GROUP_WB=0 $G $h both both
done
# S3 cohort check, S4 adherence, S5 protocol covariate, S6 two-run, S7 3-trial rule
run GROUP_WB=0 $G S3 base all
for h in H3 H2; do run GROUP_WB=0 $G $h base S4; run GROUP_WB=0 GROUP_COV=protocol $G $h base all; done
for h in H3 H1 H2 H4 H6; do run GROUP_WB=0 $G $h base S6; done
run GROUP_WB=1 $G H3 base S7
run GROUP_WB=0 $G S7ols base S7
# H3 per-run robustness form (v0.6)
run GROUP_WB=1 $G H3perrun base all
# Harvard-Oxford atlas sensitivity (reported, never inferential)
for h in H3 H1 H2 H4; do run GROUP_WB=0 GROUP_ATLAS=ho $G $h base all; done
# H5 TOST and the reliability diagnostic
run $PY Scripts/intero2023_h5_reliability.py h5
run $PY Scripts/intero2023_h5_reliability.py reliability
echo "[$(date +%T)] ALL GROUP ANALYSES DONE" >> "$LOG"
