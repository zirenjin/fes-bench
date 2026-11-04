#!/bin/bash
# PROTOCOL.md gate 1: all 521 primary entries, 6 workers on the A100.
R=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; mkdir -p data logs
cp /GenSIvePFS/users/zirenj/fes-mattersim-replication-20261102/entries.txt entries.txt
for w in 0 1 2 3 4 5; do
  awk -v n=6 -v w=$w 'NR % n == w' entries.txt > logs/worker_$w.txt
  OMP_NUM_THREADS=8 $PY scripts/gate1.py data $(cat logs/worker_$w.txt) > logs/worker_$w.log 2>&1 &
done
wait
echo "DONE entries $(ls data/*.npz | wc -l) failed $(ls data/*FAILED* 2>/dev/null | wc -l)"
