#!/bin/bash
# PROTOCOL.md: verify checkpoint hash, then run the 521 primary entries with 6 workers, then the unchanged scorers.
set -e
R=/GenSIvePFS/users/zirenj/fes-mattersim-replication-20261102; A=/GenSIvePFS/users/zirenj/fes-alexandria-polymorph-audit-20261030
AR=/GenSIvePFS/users/zirenj/fes-anchored-response-20261102; PYM=/GenSIvePFS/users/zirenj/envs/mattersim/bin/python; PYD=/root/miniconda3/envs/mat-agent/bin/python
CK=/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth
cd $R; mkdir -p mlip/MatterSim-5M logs analysis
echo "e3df9fa708725e3d453140646c7d1838324b347a3d1214cf1440522146f872b5  $CK" | sha256sum -c -
cp /GenSIvePFS/users/zirenj/fes-relaxation-channel-20261101/entries.txt entries.txt; wc -l < entries.txt
for w in 0 1 2 3 4 5; do
  awk -v n=6 -v w=$w 'NR % n == w' entries.txt > logs/worker_$w.txt
  OMP_NUM_THREADS=8 $PYM scripts/mattersim_vs_dft_phonons.py $CK - mlip/MatterSim-5M $(cat logs/worker_$w.txt) > logs/worker_$w.log 2>&1 &
done
wait
echo "entries $(ls mlip/MatterSim-5M/*.npz | wc -l) failed $(ls mlip/MatterSim-5M/*FAILED* 2>/dev/null | wc -l)"
G=$A/metadata/selection/FINAL_GROUPS.csv
$PYD $A/scripts/analyze_covariance.py $G mlip/MatterSim-5M analysis/covariance > /dev/null
$PYD $A/scripts/abstention_check.py $G mlip/MatterSim-5M analysis/ABSTENTION.json > /dev/null
$PYD $AR/scripts/anchored_alexandria.py $G mlip/MatterSim-5M analysis/ANCHORED.json > /dev/null
$PYD $AR/scripts/supp_same_subset.py $G mlip/MatterSim-5M analysis/SUPP_SAME_SUBSET.json > /dev/null
echo SCORED
