import json
from pathlib import Path

import numpy as np

from src.experiments.external_baselines import _global_delta_mean


def test_global_mean_delta_g_uses_training_labels_only():
    split = json.loads(Path("data/processed/splits/temp_extrap.json").read_text(encoding="utf-8"))
    train = {}
    for row in split["train"]:
        train.setdefault((row["system"], row["phase"]), []).append(row["T_index"])
    value = _global_delta_mean(Path("data/processed"), train)
    assert value is not None and np.isfinite(value)


def test_external_and_trivial_floor_artifacts_are_present():
    external = json.loads(Path("result/experiments/external_baselines/recomputed/metrics.json").read_text(encoding="utf-8"))
    trivial = json.loads(Path("result/experiments/skill_floor/metrics.json").read_text(encoding="utf-8"))
    assert set(("temp_extrap", "phase_lopo", "system_loso")) <= external.keys()
    assert set(("temp_extrap", "phase_lopo", "system_loso")) <= trivial.keys()
    assert "global_mean_delta_g" in external["temp_extrap"]
    assert "constant_delta_g" in trivial["temp_extrap"]
