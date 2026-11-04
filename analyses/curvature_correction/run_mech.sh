#!/bin/bash
# PROTOCOL_RE1_MECH.md: all 521 primary entries, 6 workers (runs alongside RE1).
R=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; mkdir -p data_mech logs_mech
for w in 0 1 2 3 4 5; do
  awk -v n=6 -v w=$w 'NR % n == w' entries.txt > logs_mech/worker_$w.txt
  OMP_NUM_THREADS=4 $PY scripts/re1_mech.py data_mech $(cat logs_mech/worker_$w.txt) > logs_mech/worker_$w.log 2>&1 &
done
wait
echo "DONE entries $(ls data_mech/*.npz | wc -l) failed $(ls data_mech/*FAILED* 2>/dev/null | wc -l)"
