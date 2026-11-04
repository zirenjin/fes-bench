#!/bin/bash
# DEVIATIONS_RE1 D2: rerun GPU-OOM failures with 4 workers.
R=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; RAW=/GenSIvePFS/users/zirenj/fes-alexandria-polymorph-audit-20261030/data/raw; PY=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python
cd $R; ls data_re1 | grep FAILED | sed 's/.FAILED.txt//' > logs_re1/rerun_ids.txt; rm -f data_re1/*FAILED*
for w in 0 1 2 3; do
  awk -v n=4 -v w=$w -v raw=$RAW 'NR % n == w {print raw "/" $1 ".yaml.bz2"}' logs_re1/rerun_ids.txt > logs_re1/rerun_$w.txt
  OMP_NUM_THREADS=8 $PY scripts/re1.py data_re1 $(cat logs_re1/rerun_$w.txt) > logs_re1/rerun_$w.log 2>&1 &
done
wait
echo "RERUN_DONE entries $(ls data_re1/*.npz | wc -l) failed $(ls data_re1/*FAILED* 2>/dev/null | wc -l)"
