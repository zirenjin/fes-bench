#!/bin/bash
# PROTOCOL_1c_A1.md: all 521 primary entries, 6 workers.
R=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; mkdir -p data1c logs1c
for w in 0 1 2 3 4 5; do
  awk -v n=6 -v w=$w 'NR % n == w' entries.txt > logs1c/worker_$w.txt
  OMP_NUM_THREADS=8 $PY scripts/gate1c_a1.py data1c $(cat logs1c/worker_$w.txt) > logs1c/worker_$w.log 2>&1 &
done
wait
echo "DONE entries $(ls data1c/*.npz | wc -l) failed $(ls data1c/*FAILED* 2>/dev/null | wc -l)"
