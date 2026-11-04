#!/bin/bash
# PROTOCOL.md data step: band-matched spectra for the 521 primary entries, 6 workers on the shared A100.
R=/GenSIvePFS/users/zirenj/fes-correction-routes-20261102; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; mkdir -p data/spectra logs
cp /GenSIvePFS/users/zirenj/fes-mattersim-replication-20261102/entries.txt entries.txt
for w in 0 1 2 3 4 5; do
  awk -v n=6 -v w=$w 'NR % n == w' entries.txt > logs/worker_$w.txt
  OMP_NUM_THREADS=8 $PY scripts/spectra.py data/spectra $(cat logs/worker_$w.txt) > logs/worker_$w.log 2>&1 &
done
wait
echo "DONE entries $(ls data/spectra/*.npz | wc -l) failed $(ls data/spectra/*FAILED* 2>/dev/null | wc -l)"
