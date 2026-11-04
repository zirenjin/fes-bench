# FES-Bench

Reproduce every metric, table and figure from the stored curves (no model inference, no DFT; about 20 s):

```bash
pip install -r requirements.txt
bash reproduce.sh
```

Outputs are written to `database_pairs/{tables,figures,audit}`, `finite_temperature/results` and `analyses/curvature_correction/figures`.
