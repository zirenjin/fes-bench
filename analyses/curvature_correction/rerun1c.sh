#!/bin/bash
# DEVIATIONS_1c D1: rerun the short-supercell entries with the guarded script.
G=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; R=/GenSIvePFS/users/zirenj/fes-alexandria-polymorph-audit-20261030/data/raw
cd $G; rm -f data1c/*FAILED*; args=""
for m in mp-18377 mp-553880 mp-556576 mp-561286 mp-569510 mp-627401 mp-680205 mp-9947; do test -f data1c/$m.npz || args="$args $R/$m.yaml.bz2"; done
OMP_NUM_THREADS=8 /GenSIvePFS/users/zirenj/envs/mattersim/bin/python scripts/gate1c_a1_guard.py data1c $args
echo RERUN_DONE
