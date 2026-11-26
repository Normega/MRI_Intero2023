#!/bin/bash
# run_after_firstlevel.sh -- everything after the registered first levels, in the background.
#  (a) Model A-g first levels (exploratory) and H7 per-participant coupling, 8 jobs in parallel;
#  (b) at the same time, every registered group analysis (run_group_all.sh, 8 permutation jobs);
#  (c) then the Model A-g group tests and the H7 group tests.
set -u
cd "$(dirname "$0")/.."
PY=~/miniforge3/envs/neuro/bin/python
OUT=/media/normega/FB_Cluster_Data/intero2024/derivatives/intero2023_glm
SUBS=$OUT/logs/subjects.txt
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
( sed 's/^/ag /' $SUBS; sed 's/^/h7 /' $SUBS ) | xargs -P 8 -L 1 sh -c '
  if [ "$0" = ag ]; then '$PY' Scripts/intero2023_model_ag.py "$1" > '$OUT'/logs/ag_"$1".log 2>&1 || echo "FAILED ag $1" >> '$OUT'/logs/after_failed.txt
  else '$PY' Scripts/intero2023_h7_connectivity.py subject "$1" > '$OUT'/logs/h7_"$1".log 2>&1 || echo "FAILED h7 $1" >> '$OUT'/logs/after_failed.txt; fi' &
PA=$!
GROUP_NJOBS=8 bash Scripts/run_group_all.sh
wait $PA
LOG=$OUT/logs/group.log
for h in AG AGrank AGaw; do echo "[$(date +%T)] $h" >> $LOG; GROUP_NJOBS=12 GROUP_WB=1 $PY Scripts/intero2023_group.py $h base all >> $LOG 2>&1 || echo "FAILED: $h" >> $LOG; done
echo "[$(date +%T)] H7 group" >> $LOG; $PY Scripts/intero2023_h7_connectivity.py group >> $LOG 2>&1 || echo "FAILED: H7 group" >> $LOG
echo "[$(date +%T)] EVERYTHING DONE" >> $LOG
