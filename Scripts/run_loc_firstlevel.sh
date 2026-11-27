#!/bin/bash
# run_loc_firstlevel.sh -- localizer first levels (OSF p5ztc) for all 54 participants.
set -u
REPO=/media/normega/FB_Cluster_Data/intero2024/MRI_Intero2023
PY=/home/normega/miniforge3/envs/neuro/bin/python
OUT=/media/normega/FB_Cluster_Data/intero2024/derivatives/intero2023_glm
LOGS=$OUT/logs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p $LOGS; rm -f $LOGS/loc_fl_failed.txt
cut -d, -f1,16 $REPO/FirstLevel/localizer/localizer_manifest.csv | grep ',True$' | cut -d, -f1 | sort -u > $LOGS/loc_subjects.txt
date > $LOGS/loc_firstlevel_start.txt
xargs -P 8 -I{} sh -c "$PY $REPO/Scripts/intero2023_loc_firstlevel.py {} > $LOGS/loc_fl_{}.log 2>&1 || echo 'FAILED {}' >> $LOGS/loc_fl_failed.txt" < $LOGS/loc_subjects.txt
date > $LOGS/loc_firstlevel_end.txt
n=$(ls $OUT/firstlevel/*/*_task-localizer_firstlevel.json | wc -l)
echo "localizer first levels complete: $n of $(wc -l < $LOGS/loc_subjects.txt)"
[ -s $LOGS/loc_fl_failed.txt ] && cat $LOGS/loc_fl_failed.txt
exit 0
