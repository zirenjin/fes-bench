#!/bin/bash
G=/GenSIvePFS/users/zirenj/fes-local-harmonic-gate1-20261102; R=/GenSIvePFS/users/zirenj/fes-alexandria-polymorph-audit-20261030/data/raw
cd $G; mkdir -p smoke_re1; rm -f smoke_re1/*
OMP_NUM_THREADS=8 /GenSIvePFS/users/zirenj/envs/mattersim/bin/python scripts/re1.py smoke_re1 $R/mp-22865.yaml.bz2 $R/mp-570935.yaml.bz2 $R/mp-4604.yaml.bz2
echo SMOKE_DONE
