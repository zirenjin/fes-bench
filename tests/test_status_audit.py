from __future__ import annotations

import json
from pathlib import Path

from fes_bench.status.audit import _check_metric


def test_metric_audit_accepts_reference_and_noise_shapes(tmp_path: Path) -> None:
    reference = tmp_path / "reference.json"
    reference.write_text(
        json.dumps(
            {
                "predictor": "reference",
                "G_MAE_eV_per_atom": 0.0,
                "pairs": {
                    "x": {
                        "Tc_error_K": [0.0],
                        "false_crossings": 0,
                        "missed_crossings": 0,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    noise = tmp_path / "noise.json"
    noise.write_text(
        json.dumps(
            {
                "predictor": "reference_noise:0.005",
                "pairs": {
                    "x": {
                        "delta_G_MAE_at_crossing_eV_per_atom": [0.005],
                        "false_crossings": 0,
                        "missed_crossings": 0,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    checks: dict[str, str] = {}
    _check_metric(reference, "reference", checks)
    _check_metric(noise, "reference_noise:0.005", checks)
    assert len(checks) == 2
