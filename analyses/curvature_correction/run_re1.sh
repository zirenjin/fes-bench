#!/bin/bash
# PROTOCOL_RE1.md: all 521 primary entries, 16 workers.
R=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; mkdir -p data_re1 logs_re1
for w in $(seq 0 15); do
  awk -v n=16 -v w=$w 'NR % n == w' entries.txt > logs_re1/worker_$w.txt
  OMP_NUM_THREADS=4 $PY scripts/re1.py data_re1 $(cat logs_re1/worker_$w.txt) > logs_re1/worker_$w.log 2>&1 &
done
wait
echo "DONE entries $(ls data_re1/*.npz | wc -l) failed $(ls data_re1/*FAILED* 2>/dev/null | wc -l)"
