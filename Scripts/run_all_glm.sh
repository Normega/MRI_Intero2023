#!/bin/bash
# run_all_glm.sh -- registered first levels for all participants, then everything downstream.
set -u
cd "$(dirname "$0")/.."
PY=~/miniforge3/envs/neuro/bin/python
OUT=/media/normega/FB_Cluster_Data/intero2024/derivatives/intero2023_glm
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
date > $OUT/logs/firstlevel_start.txt
xargs -P 11 -I{} sh -c "$PY Scripts/intero2023_firstlevel.py {} > $OUT/logs/fl_{}.log 2>&1 || echo 'FAILED {}' >> $OUT/logs/fl_failed.txt" < $OUT/logs/subjects.txt
date > $OUT/logs/firstlevel_end.txt
n=$(ls $OUT/firstlevel/*/*_firstlevel.json | wc -l)
echo "first levels complete: $n" > $OUT/logs/group.log
if [ "$n" -ne 47 ] || [ -s $OUT/logs/fl_failed.txt ]; then echo "FAILED: first level incomplete" >> $OUT/logs/group.log; exit 1; fi
bash Scripts/run_after_firstlevel.sh
