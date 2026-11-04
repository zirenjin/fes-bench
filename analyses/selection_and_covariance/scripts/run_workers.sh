#!/bin/bash
# Launch N parallel MLIP workers on thu-GenSi; each worker loads the model once and processes its share of entries.
# Usage: run_workers.sh <head> <n_workers>
set -e
HEAD=$1; N=$2
ROOT=/GenSIvePFS/users/zirenj/fes-alexandria-polymorph-audit-20261030
CK=/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt
PY=/root/miniconda3/envs/mat-agent/bin/python
OUT=$ROOT/mlip/$HEAD; mkdir -p $OUT $ROOT/logs
cd $ROOT
ls data/raw/*.yaml.bz2 | sort > $OUT/entries.txt
for w in $(seq 0 $((N-1))); do
  awk -v n=$N -v w=$w 'NR % n == w' $OUT/entries.txt > $OUT/worker_$w.txt
  OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 nohup $PY scripts/mlip_vs_dft_phonons.py $CK $HEAD $OUT $(cat $OUT/worker_$w.txt) > logs/${HEAD}_worker_$w.log 2>&1 &
done
echo "launched $N workers for $HEAD; entries $(wc -l < $OUT/entries.txt)"
